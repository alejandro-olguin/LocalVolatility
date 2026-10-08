// [[Rcpp::interfaces(r, cpp)]]

#include <Rcpp.h>
#include "utils.h"
#include "sigma_time.h"

using namespace Rcpp;

namespace LocalVolatility {

double american_option_lv_2d(double s_0,
                             double x_0,
                             double k,
                             double tau,
                             double r_d,
                             double r_f,
                             double q,
                             NumericMatrix sigma_s,
                             NumericMatrix sigma_x,
                             double rho,
                             String type,
                             double s_min,
                             double s_max,
                             double x_min,
                             double x_max,
                             int n_s,
                             int n_x,
                             int n_t,
                             double alpha,
                             double lambda,
                             double tolerance,
                             double sigma_t_max) {

  // Input guards
  if (n_s < 3 || n_x < 3) stop("n_s and n_x must be >= 3");
  if (n_t < 1)             stop("n_t must be >= 1");
  if (s_max <= s_min || x_max <= x_min) stop("grid bounds must be increasing");
  if (tau <= 0.0)          stop("tau must be > 0");
  if (rho < -1.0 || rho > 1.0) stop("rho must be within [-1, 1]");
  if (lambda <= 0.0)       stop("lambda must be > 0");
  if (tolerance <= 0.0)    stop("tolerance must be > 0");
  if (sigma_s.nrow() != n_s + 1 || sigma_s.ncol() < 1)
    stop("sigma_s must have n_s+1 rows (columns = calendar-time nodes on [0, sigma_t_max])");
  if (sigma_x.nrow() != n_x + 1 || sigma_x.ncol() < 1)
    stop("sigma_x must have n_x+1 rows (columns = calendar-time nodes on [0, sigma_t_max])");

  // Grids (Tavella-Randall non-uniform)
  NumericVector s = tavella_randall(s_0, alpha, s_min, s_max, n_s);
  NumericVector x = tavella_randall(x_0, alpha, x_min, x_max, n_x);

  // Non-uniform steps
  NumericVector hs(n_s), hx(n_x);
  for (int i = 0; i < n_s; ++i) hs[i] = s[i + 1] - s[i];
  for (int j = 0; j < n_x; ++j) hx[j] = x[j + 1] - x[j];

  const double dt = tau / n_t;
  // Penalty scaled relative to the identity (the 2D system is divided by dt),
  // so `lambda` has the same meaning as in the 1D solvers.
  const double half_lambda = 0.5 * lambda / dt;
  const double t_sig = (sigma_t_max > 0.0) ? sigma_t_max : tau;
  const int max_iter = 200;

  // Payoff on S*X
  NumericMatrix payoff(n_s + 1, n_x + 1);
  if (type == "call") {
    for (int i = 0; i <= n_s; ++i)
      for (int j = 0; j <= n_x; ++j)
        payoff(i, j) = std::max(s[i] * x[j] - k, 0.0);
  } else if (type == "put") {
    for (int i = 0; i <= n_s; ++i)
      for (int j = 0; j <= n_x; ++j)
        payoff(i, j) = std::max(k - s[i] * x[j], 0.0);
  } else {
    stop("type must be \"call\" or \"put\"");
  }

  NumericMatrix u(n_s + 1, n_x + 1);
  for (int i = 0; i <= n_s; ++i)
    for (int j = 0; j <= n_x; ++j)
      u(i, j) = payoff(i, j);

  // Analytic far-field Dirichlet boundaries
  auto set_S_boundaries = [&](double tau_rem,
                              NumericVector& leftS, NumericVector& rightS) {
    leftS  = NumericVector(n_x + 1);
    rightS = NumericVector(n_x + 1);
    if (type == "call") {
      for (int j = 0; j <= n_x; ++j) {
        leftS[j]  = 0.0;
        // Far field: max(European asymptote, immediate exercise). Under the
        // domestic measure S*X drifts at r_d - q, so the asymptote is
        // S*X*exp(-q*tau) - K*exp(-r_d*tau).
        rightS[j] = std::max(std::max(s[n_s] * x[j] * std::exp(-q * tau_rem)
                                - k * std::exp(-r_d * tau_rem), 0.0), payoff(n_s, j));
      }
    } else {
      for (int j = 0; j <= n_x; ++j) {
        leftS[j]  = std::max(k - s[0] * x[j], 0.0);
        rightS[j] = 0.0;
      }
    }
  };

  auto set_X_boundaries = [&](double tau_rem,
                              NumericVector& leftX, NumericVector& rightX) {
    leftX  = NumericVector(n_s + 1);
    rightX = NumericVector(n_s + 1);
    if (type == "call") {
      for (int i = 0; i <= n_s; ++i) {
        leftX[i]  = 0.0;
        rightX[i] = std::max(std::max(s[i] * x[n_x] * std::exp(-q * tau_rem)
                                - k * std::exp(-r_d * tau_rem), 0.0), payoff(i, n_x));
      }
    } else {
      for (int i = 0; i <= n_s; ++i) {
        leftX[i]  = std::max(k - s[i] * x[0], 0.0);
        rightX[i] = 0.0;
      }
    }
  };

  // Tridiagonal buffers
  NumericVector as(n_s - 1), bs(n_s - 1), cs(n_s - 1), fs(n_s - 1);
  NumericVector ax(n_x - 1), bx(n_x - 1), cx(n_x - 1), fx(n_x - 1);

  // Time stepping (backward)
  for (int n = 0; n < n_t; ++n) {
    const double tau_rem = (n_t - (n + 1)) * dt;
    // Calendar time of the new (implicit) level; sigma columns span [0, t_sig].
    const double t_new = tau - tau_rem;
    const lvtime::ColLerp cl_s = lvtime::col_lerp(t_new, t_sig, sigma_s.ncol());
    const lvtime::ColLerp cl_x = lvtime::col_lerp(t_new, t_sig, sigma_x.ncol());

    NumericVector leftS, rightS, leftX, rightX;
    set_S_boundaries(tau_rem, leftS, rightS);
    set_X_boundaries(tau_rem, leftX, rightX);

    // Penalty fixed-point iteration
    NumericMatrix u_i = clone(u);
    double rel_err = 1e6;
    int iter = 0;

    while (rel_err > tolerance && iter < max_iter) {

      // Penalty indicator (lambda/2 per substep)
      NumericMatrix p_u(n_s + 1, n_x + 1);
      for (int i = 0; i <= n_s; ++i)
        for (int j = 0; j <= n_x; ++j)
          p_u(i, j) = (u_i(i, j) < payoff(i, j)) ? half_lambda : 0.0;

      // === Pass 1: implicit in S, explicit cross ===
      NumericMatrix v(n_s + 1, n_x + 1);
      for (int j = 1; j < n_x; ++j) {
        const double L = leftS[j], R = rightS[j];

        for (int i = 1; i < n_s; ++i) {
          const double h_im1 = hs[i - 1], h_i = hs[i];
          const double Si = s[i], Xj = x[j];

          const double sigS = lvtime::at(sigma_s, i, cl_s);
          const double sigX = lvtime::at(sigma_x, j, cl_x);
          const double a2S = sigS * sigS * Si * Si;

          const double muS = r_f - q - rho * sigS * sigX;

          as[i - 1] = (muS * Si * h_i - a2S) / (h_im1 * (h_im1 + h_i));
          bs[i - 1] = 1.0 / dt
            + (a2S - muS * Si * (h_i - h_im1)) / (h_im1 * h_i)
            + 0.5 * r_d
            + 0.5 * p_u(i, j);
          cs[i - 1] = (-muS * Si * h_im1 - a2S) / (h_i * (h_im1 + h_i));

          const double denom = (h_im1 + h_i) * (hx[j - 1] + hx[j]);
          double cross = 0.0;
          if (denom != 0.0) {
            cross = 0.5 * rho * sigS * sigX * Si * Xj *
              (u(i+1,j+1) + u(i-1,j-1) - u(i-1,j+1) - u(i+1,j-1)) / denom;
          }

          fs[i - 1] = u(i, j) / dt + cross + 0.5 * p_u(i, j) * payoff(i, j);
        }

        fs[0]       -= as[0]      * L;
        fs[n_s - 2] -= cs[n_s - 2] * R;
        as[0] = 0.0; cs[n_s - 2] = 0.0;

        NumericVector v_col = thomas_algorithm(as, bs, cs, fs);
        v(0, j) = L; v(n_s, j) = R;
        for (int i = 1; i < n_s; ++i) v(i, j) = v_col[i - 1];
      }

      for (int i = 0; i <= n_s; ++i) { v(i, 0) = leftX[i]; v(i, n_x) = rightX[i]; }
      for (int j = 0; j <= n_x; ++j) { v(0, j) = leftS[j]; v(n_s, j) = rightS[j]; }

      // === Pass 2: implicit in X, explicit cross ===
      NumericMatrix u_new(n_s + 1, n_x + 1);

      NumericMatrix p_v(n_s + 1, n_x + 1);
      for (int i = 0; i <= n_s; ++i)
        for (int j = 0; j <= n_x; ++j)
          p_v(i, j) = (v(i, j) < payoff(i, j)) ? half_lambda : 0.0;

      for (int i = 1; i < n_s; ++i) {
        const double L = leftX[i], R = rightX[i];

        for (int j = 1; j < n_x; ++j) {
          const double h_jm1 = hx[j - 1], h_j = hx[j];
          const double Si = s[i], Xj = x[j];

          const double sigS = lvtime::at(sigma_s, i, cl_s);
          const double sigX = lvtime::at(sigma_x, j, cl_x);
          const double a2X = sigX * sigX * Xj * Xj;

          ax[j - 1] = ((r_d - r_f) * Xj * h_j - a2X) / (h_jm1 * (h_jm1 + h_j));
          bx[j - 1] = 1.0 / dt
            + (a2X - (r_d - r_f) * Xj * (h_j - h_jm1)) / (h_jm1 * h_j)
            + 0.5 * r_d
            + 0.5 * p_v(i, j);
          cx[j - 1] = (-(r_d - r_f) * Xj * h_jm1 - a2X) / (h_j * (h_jm1 + h_j));

          const double denom = (hs[i - 1] + hs[i]) * (h_jm1 + h_j);
          double cross = 0.0;
          if (denom != 0.0) {
            cross = 0.5 * rho * sigS * sigX * Si * Xj *
              (v(i+1,j+1) + v(i-1,j-1) - v(i-1,j+1) - v(i+1,j-1)) / denom;
          }

          fx[j - 1] = v(i, j) / dt + cross + 0.5 * p_v(i, j) * payoff(i, j);
        }

        fx[0]       -= ax[0]      * L;
        fx[n_x - 2] -= cx[n_x - 2] * R;
        ax[0] = 0.0; cx[n_x - 2] = 0.0;

        NumericVector row = thomas_algorithm(ax, bx, cx, fx);
        u_new(i, 0) = L; u_new(i, n_x) = R;
        for (int j = 1; j < n_x; ++j) u_new(i, j) = row[j - 1];
      }

      for (int i = 0; i <= n_s; ++i) { u_new(i, 0) = leftX[i]; u_new(i, n_x) = rightX[i]; }
      for (int j = 0; j <= n_x; ++j) { u_new(0, j) = leftS[j]; u_new(n_s, j) = rightS[j]; }

      // Project to obstacle and compute relative error
      double num = 0.0, den = 0.0;
      for (int i = 0; i <= n_s; ++i) {
        for (int j = 0; j <= n_x; ++j) {
          const double vproj = std::max(u_new(i, j), payoff(i, j));
          const double diff  = vproj - u_i(i, j);
          num += diff * diff;
          den += u_i(i, j) * u_i(i, j);
          u_i(i, j) = vproj;
        }
      }
      rel_err = std::sqrt(num / (den + 1e-16));
      ++iter;
    }

    if (iter >= max_iter)
      Rcpp::warning("Penalty iteration did not converge at time step %d (error = %.2e)", n, rel_err);

    u = clone(u_i);
  }

  return bilinear_interpolation(s, x, u, s_0, x_0);
}

}

//' American Option 2D (local volatility, penalty method)
//'
//' Prices an American option on an equity in a foreign currency with local
//' volatilities for both S and X, using Yanenko operator splitting and a
//' penalty-projection scheme for the early exercise constraint on a
//' non-uniform Tavella-Randall grid.
//'
//' @param s_0 Stock spot price.
//' @param x_0 FX spot price (domestic per foreign).
//' @param k Strike price (domestic currency).
//' @param tau Time to expiry (in years).
//' @param r_d Risk-free rate (domestic).
//' @param r_f Risk-free rate (foreign).
//' @param q Dividend yield.
//' @param sigma_s Local volatility matrix for the stock with \code{n_s + 1} rows.
//'   Columns are uniformly spaced in calendar time on \code{[0, sigma_t_max]}:
//'   column 1 is the valuation date, the last column is \code{sigma_t_max}.
//'   Any number of columns is accepted; values are interpolated in time.
//' @param sigma_x Local volatility matrix for the FX with \code{n_x + 1} rows,
//'   same time convention as \code{sigma_s}.
//' @param rho Correlation between the stock and the FX in \code{[-1, 1]}.
//' @param type Either \code{"call"} or \code{"put"}.
//' @param s_min,s_max Min/Max of the stock grid.
//' @param x_min,x_max Min/Max of the FX grid.
//' @param n_s Number of intervals in stock grid (\code{n_s + 1} nodes).
//' @param n_x Number of intervals in FX grid (\code{n_x + 1} nodes).
//' @param n_t Number of time steps.
//' @param alpha Grid clustering parameter for Tavella-Randall grids (> 0).
//' @param lambda Penalty parameter (> 0). Scaled internally by \code{1/dt} so it
//'   has the same meaning as in the 1D solvers (relative to the identity).
//' @param tolerance Convergence tolerance for penalty iterations.
//' @param sigma_t_max Calendar horizon (years) spanned by the sigma columns.
//'   Defaults to \code{tau} when not positive, i.e. the matrices describe
//'   exactly this option's life. Pass the common horizon when one matrix is
//'   shared across options of different maturities.
//'
//' @return Option price as a numeric scalar.
//'
//' @details
//' Domestic-measure dynamics with pointwise volatilities:
//' \eqn{\mu_S = r_f - q - \rho \sigma_S \sigma_X} in S, \eqn{r_d - r_f} in X,
//' discounting at \eqn{r_d}. Far-field Dirichlet values are the maximum of the
//' European asymptote \eqn{SX e^{-q\tau} - K e^{-r_d\tau}} and the payoff. The
//' penalty parameter is split as \eqn{\lambda/2} across the two
//' operator-splitting sub-steps. Maximum 200 penalty iterations
//' per time step with a warning on non-convergence.
//'
//' @examples
//' SigmaS <- matrix(0.2, nrow = 21, ncol = 25)
//' SigmaX <- matrix(0.15, nrow = 21, ncol = 25)
//' american_option_lv_2d(100, 20, 2000, 1, 0.05, 0.02, 0.01,
//'                       SigmaS, SigmaX, 0.3, "put",
//'                       10, 300, 1, 60, 20, 20, 25, 3, 1e4, 1e-8)
//'
//' @export
// [[Rcpp::export]]
double american_option_lv_2d(double s_0,
                             double x_0,
                             double k,
                             double tau,
                             double r_d,
                             double r_f,
                             double q,
                             NumericMatrix sigma_s,
                             NumericMatrix sigma_x,
                             double rho,
                             String type,
                             double s_min,
                             double s_max,
                             double x_min,
                             double x_max,
                             int n_s,
                             int n_x,
                             int n_t,
                             double alpha,
                             double lambda,
                             double tolerance,
                             double sigma_t_max = -1.0) {

  return LocalVolatility::american_option_lv_2d(s_0, x_0, k, tau, r_d, r_f, q,
                                                sigma_s, sigma_x, rho, type,
                                                s_min, s_max, x_min, x_max,
                                                n_s, n_x, n_t, alpha,
                                                lambda, tolerance, sigma_t_max);
}
