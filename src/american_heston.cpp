// [[Rcpp::interfaces(r, cpp)]]

#include <Rcpp.h>
#include "utils.h"

using namespace Rcpp;

namespace LocalVolatility {

double american_option_heston(double s_0,
                              double k,
                              double tau,
                              double r_d,
                              double q,
                              double kappa,
                              double theta,
                              double xi,
                              double rho,
                              double v_0,
                              String type,
                              double s_min,
                              double s_max,
                              double v_min,
                              double v_max,
                              int n_s,
                              int n_v,
                              int n_t,
                              double alpha,
                              double lambda,
                              double tolerance) {

  // Input guards
  if (n_s < 3 || n_v < 3)       stop("n_s and n_v must be >= 3");
  if (n_t < 1)                   stop("n_t must be >= 1");
  if (s_max <= s_min)            stop("s_max must be > s_min");
  if (v_max <= v_min)            stop("v_max must be > v_min");
  if (v_min <= 0.0)              stop("v_min must be > 0");
  if (tau <= 0.0)                stop("tau must be > 0");
  if (rho < -1.0 || rho > 1.0)  stop("rho must be within [-1, 1]");
  if (kappa <= 0.0)              stop("kappa must be > 0");
  if (theta <= 0.0)              stop("theta must be > 0");
  if (xi <= 0.0)                 stop("xi must be > 0");
  if (v_0 <= 0.0)               stop("v_0 must be > 0");
  if (s_0 <= 0.0 || k <= 0.0)   stop("s_0 and k must be > 0");
  if (lambda <= 0.0)             stop("lambda must be > 0");
  if (tolerance <= 0.0)          stop("tolerance must be > 0");

  const bool is_call = (type == "call");
  if (!is_call && type != "put") stop("type must be \"call\" or \"put\"");

  // Warn if Feller condition violated
  if (2.0 * kappa * theta <= xi * xi)
    Rcpp::warning("Feller condition 2*kappa*theta > xi^2 is violated; "
                  "variance may reach zero.");

  // Grids (Tavella-Randall non-uniform)
  NumericVector s = tavella_randall(s_0, alpha, s_min, s_max, n_s);
  NumericVector v = tavella_randall(v_0, alpha, v_min, v_max, n_v);

  // Non-uniform steps
  NumericVector hs(n_s), hv(n_v);
  for (int i = 0; i < n_s; ++i) hs[i] = s[i + 1] - s[i];
  for (int j = 0; j < n_v; ++j) hv[j] = v[j + 1] - v[j];

  const double dt = tau / n_t;
  const double half_lambda = 0.5 * lambda;
  const int max_iter = 200;

  // Payoff (depends only on S, constant across v)
  NumericMatrix payoff(n_s + 1, n_v + 1);
  for (int i = 0; i <= n_s; ++i) {
    const double p = is_call ? std::max(s[i] - k, 0.0)
                             : std::max(k - s[i], 0.0);
    for (int j = 0; j <= n_v; ++j) payoff(i, j) = p;
  }

  NumericMatrix u = clone(payoff);

  // Black-Scholes for v_max boundary
  auto bs_price = [&](double Si, double tau_rem, double sigma_bs) -> double {
    if (tau_rem <= 0.0)
      return is_call ? std::max(Si - k, 0.0) : std::max(k - Si, 0.0);
    if (sigma_bs <= 0.0) {
      const double F = Si * std::exp((r_d - q) * tau_rem);
      const double disc = std::exp(-r_d * tau_rem);
      return is_call ? disc * std::max(F - k, 0.0)
                     : disc * std::max(k - F, 0.0);
    }
    const double sqrt_t = std::sqrt(tau_rem);
    const double d1 = (std::log(Si / k) + (r_d - q + 0.5 * sigma_bs * sigma_bs) * tau_rem)
                     / (sigma_bs * sqrt_t);
    const double d2 = d1 - sigma_bs * sqrt_t;
    if (is_call)
      return Si * std::exp(-q * tau_rem) * R::pnorm5(d1, 0, 1, 1, 0)
           - k * std::exp(-r_d * tau_rem) * R::pnorm5(d2, 0, 1, 1, 0);
    else
      return k * std::exp(-r_d * tau_rem) * R::pnorm5(-d2, 0, 1, 1, 0)
           - Si * std::exp(-q * tau_rem) * R::pnorm5(-d1, 0, 1, 1, 0);
  };

  // S-direction boundary conditions (independent of v)
  auto set_S_boundaries = [&](double tau_rem,
                              NumericVector& leftS, NumericVector& rightS) {
    leftS  = NumericVector(n_v + 1);
    rightS = NumericVector(n_v + 1);
    if (is_call) {
      for (int j = 0; j <= n_v; ++j) {
        leftS[j]  = 0.0;
        rightS[j] = std::max(s[n_s] * std::exp(-q * tau_rem)
                            - k * std::exp(-r_d * tau_rem), 0.0);
      }
    } else {
      for (int j = 0; j <= n_v; ++j) {
        leftS[j]  = std::max(k - s[0], 0.0);
        rightS[j] = 0.0;
      }
    }
  };

  // v-direction boundary conditions
  auto set_V_boundaries = [&](double tau_rem,
                              NumericVector& leftV, NumericVector& rightV) {
    leftV  = NumericVector(n_s + 1);
    rightV = NumericVector(n_s + 1);
    const double sigma_vmax = std::sqrt(v[n_v]);
    for (int i = 0; i <= n_s; ++i) {
      leftV[i]  = u(i, 0); // placeholder for v_min
      rightV[i] = bs_price(s[i], tau_rem, sigma_vmax);
    }
  };

  const double muS_coeff = r_d - q;

  // Tridiagonal buffers
  NumericVector aS(n_s - 1), bS(n_s - 1), cS(n_s - 1), fS(n_s - 1);
  NumericVector aV(n_v - 1), bV(n_v - 1), cV(n_v - 1), fV(n_v - 1);

  // Time stepping (backward)
  for (int n = 0; n < n_t; ++n) {
    const double tau_rem = (n_t - (n + 1)) * dt;

    NumericVector leftS, rightS, leftV, rightV;
    set_S_boundaries(tau_rem, leftS, rightS);
    set_V_boundaries(tau_rem, leftV, rightV);

    // Penalty fixed-point iteration
    NumericMatrix u_i = clone(u);
    double rel_err = 1e6;
    int iter = 0;

    while (rel_err > tolerance && iter < max_iter) {

      // Penalty indicator from current iterate (lambda/2 per substep)
      NumericMatrix p_u(n_s + 1, n_v + 1);
      for (int i = 0; i <= n_s; ++i)
        for (int j = 0; j <= n_v; ++j)
          p_u(i, j) = (u_i(i, j) < payoff(i, j)) ? half_lambda : 0.0;

      // === Pass 1: implicit in S, explicit cross ===
      NumericMatrix w(n_s + 1, n_v + 1);
      for (int j = 1; j < n_v; ++j) {
        const double L = leftS[j], R = rightS[j];
        const double vj = v[j];

        for (int i = 1; i < n_s; ++i) {
          const double h_im1 = hs[i - 1], h_i = hs[i];
          const double Si = s[i];
          const double a2S = vj * Si * Si;
          const double muS = muS_coeff * Si;

          aS[i - 1] = (muS * h_i - a2S) / (h_im1 * (h_im1 + h_i));
          bS[i - 1] = 1.0 / dt
            + (a2S - muS * (h_i - h_im1)) / (h_im1 * h_i)
            + 0.5 * r_d
            + 0.5 * p_u(i, j);
          cS[i - 1] = (-muS * h_im1 - a2S) / (h_i * (h_im1 + h_i));

          const double denom = (h_im1 + h_i) * (hv[j - 1] + hv[j]);
          double cross = 0.0;
          if (denom != 0.0) {
            cross = 0.5 * rho * xi * vj * Si *
              (u(i+1,j+1) + u(i-1,j-1) - u(i-1,j+1) - u(i+1,j-1)) / denom;
          }

          fS[i - 1] = u(i, j) / dt + cross + 0.5 * p_u(i, j) * payoff(i, j);
        }

        fS[0]       -= aS[0]       * L;
        fS[n_s - 2] -= cS[n_s - 2] * R;
        aS[0] = 0.0; cS[n_s - 2] = 0.0;

        NumericVector col = thomas_algorithm(aS, bS, cS, fS);
        w(0, j) = L; w(n_s, j) = R;
        for (int i = 1; i < n_s; ++i) w(i, j) = col[i - 1];
      }

      for (int i = 0; i <= n_s; ++i) { w(i, 0) = leftV[i]; w(i, n_v) = rightV[i]; }
      for (int j = 0; j <= n_v; ++j) { w(0, j) = leftS[j]; w(n_s, j) = rightS[j]; }

      // === Pass 2: implicit in v, explicit cross ===
      NumericMatrix u_new(n_s + 1, n_v + 1);

      NumericMatrix p_w(n_s + 1, n_v + 1);
      for (int i = 0; i <= n_s; ++i)
        for (int j = 0; j <= n_v; ++j)
          p_w(i, j) = (w(i, j) < payoff(i, j)) ? half_lambda : 0.0;

      for (int i = 1; i < n_s; ++i) {
        const double Si = s[i];

        for (int j = 1; j < n_v; ++j) {
          const double h_jm1 = hv[j - 1], h_j = hv[j];
          const double vj = v[j];
          const double a2V = xi * xi * vj;
          const double muV = kappa * (theta - vj);

          aV[j - 1] = (muV * h_j - a2V) / (h_jm1 * (h_jm1 + h_j));
          bV[j - 1] = 1.0 / dt
            + (a2V - muV * (h_j - h_jm1)) / (h_jm1 * h_j)
            + 0.5 * r_d
            + 0.5 * p_w(i, j);
          cV[j - 1] = (-muV * h_jm1 - a2V) / (h_j * (h_jm1 + h_j));

          const double denom = (hs[i - 1] + hs[i]) * (h_jm1 + h_j);
          double cross = 0.0;
          if (denom != 0.0) {
            cross = 0.5 * rho * xi * vj * Si *
              (w(i+1,j+1) + w(i-1,j-1) - w(i-1,j+1) - w(i+1,j-1)) / denom;
          }

          fV[j - 1] = w(i, j) / dt + cross + 0.5 * p_w(i, j) * payoff(i, j);
        }

        const double LV = leftV[i];
        const double RV = rightV[i];

        fV[0]       -= aV[0]       * LV;
        fV[n_v - 2] -= cV[n_v - 2] * RV;
        aV[0] = 0.0; cV[n_v - 2] = 0.0;

        NumericVector row = thomas_algorithm(aV, bV, cV, fV);
        u_new(i, 0) = LV; u_new(i, n_v) = RV;
        for (int j = 1; j < n_v; ++j) u_new(i, j) = row[j - 1];
      }

      // Set boundaries
      for (int j = 0; j <= n_v; ++j) { u_new(0, j) = leftS[j]; u_new(n_s, j) = rightS[j]; }
      for (int i = 0; i <= n_s; ++i) { u_new(i, n_v) = rightV[i]; }

      // v_min boundary: linear extrapolation from j=1, j=2, clamped to non-negative
      for (int i = 0; i <= n_s; ++i) {
        u_new(i, 0) = std::max(0.0, u_new(i, 1) - hv[0] * (u_new(i, 2) - u_new(i, 1)) / hv[1]);
      }

      // Project to obstacle and compute relative error
      double num = 0.0, den = 0.0;
      for (int i = 0; i <= n_s; ++i) {
        for (int j = 0; j <= n_v; ++j) {
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

    std::swap(u, u_i);
  }

  return bilinear_interpolation(s, v, u, s_0, v_0);
}

}

//' American Option (Heston stochastic volatility, penalty method)
//'
//' Prices an American option under the Heston stochastic volatility model
//' using a 2D finite-difference PDE with Yanenko operator splitting and a
//' penalty-projection scheme for the early exercise constraint.
//'
//' @param s_0 Stock spot price (> 0).
//' @param k Strike price (> 0).
//' @param tau Time to expiry in years (> 0).
//' @param r_d Risk-free rate.
//' @param q Dividend yield.
//' @param kappa Mean reversion speed (> 0).
//' @param theta Long-run variance (> 0).
//' @param xi Volatility of variance (vol-of-vol, > 0).
//' @param rho Correlation between stock and variance in \code{[-1, 1]}.
//' @param v_0 Initial variance (> 0).
//' @param type Either \code{"call"} or \code{"put"}.
//' @param s_min,s_max Min/Max of the stock grid.
//' @param v_min Min of the variance grid (must be > 0, e.g. 0.001).
//' @param v_max Max of the variance grid.
//' @param n_s Number of intervals in stock grid (\code{n_s + 1} nodes, >= 3).
//' @param n_v Number of intervals in variance grid (\code{n_v + 1} nodes, >= 3).
//' @param n_t Number of time steps (>= 1).
//' @param alpha Grid clustering parameter for Tavella-Randall grids (> 0).
//' @param lambda Penalty parameter (> 0, typically 1e4 to 1e6).
//' @param tolerance Convergence tolerance for penalty iterations.
//'
//' @return Option price as a numeric scalar.
//'
//' @details
//' Solves the Heston LCP \eqn{\min(-\mathcal{L}V, V - \phi) = 0} where
//' the differential operator is:
//' \deqn{\mathcal{L}V = V_t + \frac{1}{2}vS^2 V_{SS} + \rho\xi v S V_{Sv}
//'   + \frac{1}{2}\xi^2 v V_{vv} + (r_d - q)S V_S
//'   + \kappa(\theta - v) V_v - r_d V.}
//' The penalty parameter is split as \eqn{\lambda/2} across the two
//' operator-splitting sub-steps. Maximum 200 penalty iterations per time
//' step with a warning on non-convergence.
//'
//' @examples
//' american_option_heston(100, 100, 1, 0.05, 0, 2, 0.04, 0.5, -0.7, 0.04,
//'                        "put", 20, 300, 0.001, 1.0, 80, 40, 100, 3,
//'                        1e4, 1e-8)
//'
//' @export
// [[Rcpp::export]]
double american_option_heston(double s_0,
                              double k,
                              double tau,
                              double r_d,
                              double q,
                              double kappa,
                              double theta,
                              double xi,
                              double rho,
                              double v_0,
                              String type,
                              double s_min,
                              double s_max,
                              double v_min,
                              double v_max,
                              int n_s,
                              int n_v,
                              int n_t,
                              double alpha,
                              double lambda,
                              double tolerance) {

  return LocalVolatility::american_option_heston(
    s_0, k, tau, r_d, q, kappa, theta, xi, rho, v_0, type,
    s_min, s_max, v_min, v_max, n_s, n_v, n_t, alpha,
    lambda, tolerance);
}
