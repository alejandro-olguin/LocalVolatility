// [[Rcpp::interfaces(r, cpp)]]

#include <Rcpp.h>
#include "utils.h"

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
                             double tolerance) {

  // --- sanity checks ---
  if (n_s < 3 || n_x < 3) stop("n_s and n_x must be >= 3");
  if (n_t < 1)            stop("n_t must be >= 1");
  if (s_max <= s_min || x_max <= x_min) stop("grid bounds must be increasing");
  if (tau <= 0.0)         stop("tau must be > 0");

  // Grids (Tavella–Randall non-uniform)
  NumericVector s = tavella_randall(s_0, alpha, s_min, s_max, n_s);
  NumericVector x = tavella_randall(x_0, alpha, x_min, x_max, n_x);

  // Non-uniform steps (forward differences)
  NumericVector hs(n_s), hx(n_x);
  for (int i = 0; i < n_s; ++i) hs[i] = s[i + 1] - s[i];
  for (int j = 0; j < n_x; ++j) hx[j] = x[j + 1] - x[j];

  const double dt = tau / n_t;
  const double half_lambda = 0.5 * lambda; // split penalty over the two sub-steps

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

  // Start from payoff (terminal condition at t = tau)
  NumericMatrix u(n_s + 1, n_x + 1);
  for (int i = 0; i <= n_s; ++i)
    for (int j = 0; j <= n_x; ++j)
      u(i, j) = payoff(i, j);

  // Tridiagonal buffers (S- and X- passes)
  NumericVector as(n_s - 1), bs(n_s - 1), cs(n_s - 1), fs(n_s - 1);
  NumericVector ax(n_x - 1), bx(n_x - 1), cx(n_x - 1), fx(n_x - 1);

  // Helpers: analytic Dirichlet boundaries for product options
  auto set_S_boundaries = [&](double tau_rem,
                              NumericVector& leftS, NumericVector& rightS) {
    leftS = NumericVector(n_x + 1);
    rightS = NumericVector(n_x + 1);
    if (type == "call") {
      // S -> 0 or X -> 0 => value ~ 0; S -> ∞: SX e^{-(q+r_f) tau} - K e^{-r_d tau}
      for (int j = 0; j <= n_x; ++j) {
        leftS[j]  = 0.0;
        double up = std::max(s[n_s] * x[j] * std::exp(-q * tau_rem)
                               - k * std::exp(-r_d * tau_rem), 0.0);
        rightS[j] = up;
      }
    } else { // put
      // S -> 0: ~ K e^{-r_d tau}, S -> ∞: 0
      for (int j = 0; j <= n_x; ++j) {
        leftS[j]  = k * std::exp(-r_d * tau_rem);
        rightS[j] = 0.0;
      }
    }
  };

  auto set_X_boundaries = [&](double tau_rem,
                              NumericVector& leftX, NumericVector& rightX) {
    leftX = NumericVector(n_s + 1);
    rightX = NumericVector(n_s + 1);
    if (type == "call") {
      // X -> 0 => 0; X -> ∞: SX e^{-(q+r_f) tau} - K e^{-r_d tau}
      for (int i = 0; i <= n_s; ++i) {
        leftX[i]  = 0.0;
        double up = std::max(s[i] * x[n_x] * std::exp(-q * tau_rem)
                               - k * std::exp(-r_d * tau_rem), 0.0);
        rightX[i] = up;
      }
    } else { // put
      // X -> 0: ~ K e^{-r_d tau}, X -> ∞: 0
      for (int i = 0; i <= n_s; ++i) {
        leftX[i]  = k * std::exp(-r_d * tau_rem);
        rightX[i] = 0.0;
      }
    }
  };

  // Time stepping (backward): n = 0..n_t-1
  for (int n = 0; n < n_t; ++n) {

    // Time remaining AFTER completing this full step (for boundaries)
    const double tau_rem = (n_t - (n + 1)) * dt;

    // Prepare boundary vectors (Dirichlet) for both passes
    NumericVector leftS, rightS, leftX, rightX;
    set_S_boundaries(tau_rem, leftS, rightS);
    set_X_boundaries(tau_rem, leftX, rightX);

    // Iterative penalty-projection per time level
    NumericMatrix u_i = clone(u);
    double rel_err = 1e6;
    int iters = 0;

    while (rel_err > tolerance) {
      // Penalty from current guess
      NumericMatrix p_u(n_s + 1, n_x + 1);
      for (int i = 0; i <= n_s; ++i)
        for (int j = 0; j <= n_x; ++j)
          p_u(i, j) = (u_i(i, j) < payoff(i, j)) ? half_lambda : 0.0; // λ/2 here

      // ===== Pass 1: implicit in S (solve per fixed j), explicit in X and cross =====
      NumericMatrix v(n_s + 1, n_x + 1);
      for (int j = 1; j < n_x; ++j) {
        // Known boundary values in S for this column j
        double L = leftS[j];
        double R = rightS[j];

        for (int i = 1; i < n_s; ++i) {
          double h_im1 = hs[i - 1], h_i = hs[i];
          double Si = s[i], Xj = x[j];

          double sigS = sigma_s(i, n);
          double sigX = sigma_x(j, n);

          // Non-uniform S-coefficients (sub/diag/super) with domestic-measure drift
          double muS = r_f - q - rho * sigS * sigX;
          as[i - 1] = ( (muS * Si * h_i) - std::pow(sigS * Si, 2.0) ) / (h_im1 * (h_im1 + h_i));
          bs[i - 1] = 1.0 / dt
          + (std::pow(sigS * Si, 2.0) - muS * Si * (h_i - h_im1)) / (h_im1 * h_i)
            + 0.5 * r_d
            + 0.5 * p_u(i, j);
            cs[i - 1] = ( -(muS * Si * h_im1) - std::pow(sigS * Si, 2.0) ) / (h_i * (h_im1 + h_i));

            // Explicit cross term (centered mixed derivative)
            double denom_cross = (h_im1 * hx[j - 1] + h_i * hx[j] + h_i * hx[j - 1] + h_im1 * hx[j]);
            double cross = 0.0;
            if (denom_cross != 0.0) {
              cross = 0.5 * rho * sigS * sigX * Si * Xj *
                (u(i + 1, j + 1) + u(i - 1, j - 1) - u(i - 1, j + 1) - u(i + 1, j - 1)) / denom_cross;
            }

            fs[i - 1] = u(i, j) / dt + cross + 0.5 * p_u(i, j) * payoff(i, j);
        }

        // Boundary contributions to RHS, then zero those coefficients
        fs[0]         -= as[0]       * L;
        fs[n_s - 2]   -= cs[n_s - 2] * R;
        as[0]          = 0.0;
        cs[n_s - 2]    = 0.0;

        // Solve in S
        NumericVector v_col = thomas_algorithm(as, bs, cs, fs);
        v(0, j)    = L;
        v(n_s, j)  = R;
        for (int i = 1; i < n_s; ++i) v(i, j) = v_col[i - 1];
      }

      // Enforce X-boundaries (Dirichlet) on v’s outer columns
      for (int i = 0; i <= n_s; ++i) { v(i, 0) = leftX[i];  v(i, n_x) = rightX[i]; }
      // Enforce S-boundaries on v’s outer rows
      for (int j = 0; j <= n_x; ++j) { v(0, j) = leftS[j];  v(n_s, j) = rightS[j]; }

      // ===== Pass 2: implicit in X (solve per fixed i), explicit in S and cross =====
      NumericMatrix u_new(n_s + 1, n_x + 1);

      // Penalty for this pass from v (fresh)
      NumericMatrix p_v(n_s + 1, n_x + 1);
      for (int i = 0; i <= n_s; ++i)
        for (int j = 0; j <= n_x; ++j)
          p_v(i, j) = (v(i, j) < payoff(i, j)) ? half_lambda : 0.0; // λ/2 here

      for (int i = 1; i < n_s; ++i) {
        // Known boundary values in X for this row i
        double L = leftX[i];
        double R = rightX[i];

        for (int j = 1; j < n_x; ++j) {
          double h_jm1 = hx[j - 1], h_j = hx[j];
          double Si = s[i], Xj = x[j];

          double sigS = sigma_s(i, n);
          double sigX = sigma_x(j, n);

          // Non-uniform X-coefficients (sub/diag/super)
          ax[j - 1] = ((r_d - r_f) * Xj * h_j - std::pow(sigX * Xj, 2.0)) / (h_jm1 * (h_jm1 + h_j));
          bx[j - 1] = 1.0 / dt
          + (std::pow(sigX * Xj, 2.0) - (r_d - r_f) * Xj * (h_j - h_jm1)) / (h_jm1 * h_j)
            + 0.5 * r_d
            + 0.5 * p_v(i, j);
            cx[j - 1] = (-(r_d - r_f) * Xj * h_jm1 - std::pow(sigX * Xj, 2.0)) / (h_j * (h_jm1 + h_j));

            // Explicit cross term using v
            double denom_cross = (h_jm1 * hs[i - 1] + h_j * hs[i] + h_jm1 * hs[i - 1] + h_j * hs[i]);
            double cross = 0.0;
            if (denom_cross != 0.0) {
              cross = 0.5 * rho * sigS * sigX * Si * Xj *
                (v(i + 1, j + 1) + v(i - 1, j - 1) - v(i - 1, j + 1) - v(i + 1, j - 1)) / denom_cross;
            }

            fx[j - 1] = v(i, j) / dt + cross + 0.5 * p_v(i, j) * payoff(i, j);
        }

        // Boundary contributions to RHS in X, then zero
        fx[0]         -= ax[0]       * L;
        fx[n_x - 2]   -= cx[n_x - 2] * R;
        ax[0]          = 0.0;
        cx[n_x - 2]    = 0.0;

        // Solve in X
        NumericVector row = thomas_algorithm(ax, bx, cx, fx);
        u_new(i, 0)   = L;
        u_new(i, n_x) = R;
        for (int j = 1; j < n_x; ++j) u_new(i, j) = row[j - 1];
      }

      // Enforce boundaries on the outer frame of u_new
      for (int i = 0; i <= n_s; ++i) { u_new(i, 0) = leftX[i];  u_new(i, n_x) = rightX[i]; }
      for (int j = 0; j <= n_x; ++j) { u_new(0, j) = leftS[j];  u_new(n_s, j) = rightS[j]; }

      // Project to obstacle (American) and compute relative error
      double num = 0.0, den = 0.0;
      for (int i = 0; i <= n_s; ++i) {
        for (int j = 0; j <= n_x; ++j) {
          double vproj = std::max(u_new(i, j), payoff(i, j));
          double diff  = vproj - u_i(i, j);
          num += diff * diff;
          den += u_i(i, j) * u_i(i, j);
          u_i(i, j) = vproj;
        }
      }
      rel_err = std::sqrt(num / (den + 1e-16));
      ++iters;
    }

    u = clone(u_i);
  }

  // Price at (s_0, x_0) by bilinear interpolation on the (s,x) grid
  double result = bilinear_interpolation(s, x, u, s_0, x_0);
  return result;
}
}


//' American Option 2D (local volatility, penalty method)
//'
//' Prices an American option on an equity in a foreign currency with local volatilities for both S and X.
//' Uses operator splitting and a penalty-projection scheme for the obstacle.
//'
//' @param s_0 Stock spot price.
//' @param x_0 FX spot price (domestic per foreign).
//' @param k Strike price (domestic currency).
//' @param tau Time to expiry (in years).
//' @param r_d Risk-free rate (domestic).
//' @param r_f Risk-free rate (foreign).
//' @param q Dividend yield.
//' @param sigma_s Local volatility matrix for the stock of size (n_s + 1) x n_t.
//' @param sigma_x Local volatility matrix for the FX of size (n_x + 1) x n_t.
//' @param rho Correlation between the stock and the FX in `[-1, 1]`.
//' @param type Either "call" or "put".
//' @param s_min,s_max Min/Max of the stock grid.
//' @param x_min,x_max Min/Max of the FX grid.
//' @param n_s Number of intervals in stock grid (n_s + 1 nodes).
//' @param n_x Number of intervals in FX grid (n_x + 1 nodes).
//' @param n_t Number of time steps.
//' @param alpha Grid clustering parameter for Tavella–Randall grids (> 0).
//' @param lambda Penalty parameter (> 1).
//' @param tolerance Relative error tolerance for the penalty fixed-point iterations.
//'
//' @details
//' Domestic-measure drift and discounting as in the constant-vol case, using local pointwise vols:
//' mu_S = r_f - q - rho * sigma_S * sigma_X in S, r_d - r_f in X, discounting at r_d. Far-field
//' boundaries use exp(-q * tau) on S·X and exp(-r_d * tau) on K.
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
                             double tolerance) {


  return LocalVolatility::american_option_lv_2d(s_0,
                                                x_0,
                                                k,
                                                tau,
                                                r_d,
                                                r_f,
                                                q,
                                                sigma_s,
                                                sigma_x,
                                                rho,
                                                type,
                                                s_min,
                                                s_max,
                                                x_min,
                                                x_max,
                                                n_s,
                                                n_x,
                                                n_t,
                                                alpha,
                                                lambda,
                                                tolerance);

}
