// [[Rcpp::interfaces(r, cpp)]]

#include <Rcpp.h>
#include "utils.h"

using namespace Rcpp;

namespace LocalVolatility {

double american_option_lv(double s_0,
                          double k,
                          double tau,
                          double r_d,
                          double q,
                          NumericMatrix sigma,
                          String type,
                          double s_min,
                          double s_max,
                          int n_s,
                          int n_t,
                          double lambda,
                          double tolerance) {

  // Input guards
  if (n_s < 2 || n_t < 1) stop("n_s >= 2 and n_t >= 1 required");
  if (s_max <= s_min)     stop("s_max must be > s_min");
  if (tau <= 0.0)         stop("tau must be > 0");
  if (lambda <= 0.0)      stop("lambda must be > 0");
  if (tolerance <= 0.0)   stop("tolerance must be > 0");
  if (sigma.nrow() != n_s + 1 || sigma.ncol() != n_t + 1)
    stop("sigma must have dimensions (n_s+1) x (n_t+1)");

  const double ds = (s_max - s_min) / n_s;
  const double dt = tau / n_t;
  const int max_iter = 200;

  // Storage
  NumericVector old_prices(n_s + 1), new_prices(n_s + 1),
    payoff(n_s + 1), u_bound(n_t + 1), l_bound(n_t + 1);

  // Payoff and boundary conditions
  if (type == "call") {
    for (int i = 0; i <= n_s; ++i)
      payoff[i] = std::max(s_min + i * ds - k, 0.0);

    for (int n = 0; n <= n_t; ++n) {
      const double t_rem = (n_t - n) * dt;
      u_bound[n] = 0.0;  // S -> 0
      // Far-field European asymptote (penalty handles early exercise)
      l_bound[n] = std::max(s_max * std::exp(-q * t_rem)
                              - k * std::exp(-r_d * t_rem), 0.0);
    }

  } else if (type == "put") {
    for (int i = 0; i <= n_s; ++i)
      payoff[i] = std::max(k - (s_min + i * ds), 0.0);

    for (int n = 0; n <= n_t; ++n) {
      u_bound[n] = k;      // S -> 0: immediate exercise value
      l_bound[n] = 0.0;    // S -> +inf
    }

  } else {
    stop("type must be \"call\" or \"put\"");
  }

  // Terminal condition
  old_prices = payoff;

  // Tridiagonal system buffers
  NumericVector a(n_s - 1), b(n_s - 1), c(n_s - 1), d(n_s - 1);

  for (int n = n_t - 1; n >= 0; --n) {

    new_prices[0]   = u_bound[n];
    new_prices[n_s] = l_bound[n];

    // Build base tridiagonal system (without penalty)
    for (int i = 1; i < n_s; ++i) {
      const double S_i = s_min + i * ds;

      const double sig_n   = sigma(i, n);
      const double sig_np1 = sigma(i, n + 1);

      const double alpha_n  = 0.5 * sig_n * sig_n * S_i * S_i;
      const double alpha_np = 0.5 * sig_np1 * sig_np1 * S_i * S_i;

      // Implicit part (t_n)
      a[i - 1] = -0.5 * dt * (alpha_n / (ds * ds) - (r_d - q) * S_i / (2.0 * ds));
      b[i - 1] =  1.0      + 0.5 * dt * (2.0 * alpha_n / (ds * ds) + r_d);
      c[i - 1] = -0.5 * dt * (alpha_n / (ds * ds) + (r_d - q) * S_i / (2.0 * ds));

      // RHS (t_{n+1})
      d[i - 1] =  0.5 * dt * (alpha_np / (ds * ds) - (r_d - q) * S_i / (2.0 * ds)) * old_prices[i - 1]
        + (1.0 - 0.5 * dt * (2.0 * alpha_np / (ds * ds) + r_d))                     * old_prices[i]
        + 0.5 * dt * (alpha_np / (ds * ds) + (r_d - q) * S_i / (2.0 * ds))          * old_prices[i + 1];
    }

    // Boundary contributions to RHS and zero touching off-diagonals
    d[0]       -= a[0]      * u_bound[n];
    d[n_s - 2] -= c[n_s - 2] * l_bound[n];
    a[0]        = 0.0;
    c[n_s - 2]  = 0.0;

    // Penalty iteration initialization
    NumericVector temp_prices = old_prices[Rcpp::Range(1, n_s - 1)];

    // Base copies (avoid accumulating penalty on b/d)
    const NumericVector a0 = clone(a);
    const NumericVector b0 = clone(b);
    const NumericVector c0 = clone(c);
    const NumericVector d0 = clone(d);

    double error = std::numeric_limits<double>::infinity();
    int iter = 0;

    while (error > tolerance && iter < max_iter) {
      // Obstacle violation indicator (interior)
      NumericVector p_v(n_s - 1);
      for (int i = 0; i < n_s - 1; ++i) {
        p_v[i] = (temp_prices[i] < payoff[i + 1]) ? lambda : 0.0;
      }

      // Rebuild system for this iteration (do NOT accumulate)
      NumericVector bb = clone(b0);
      NumericVector dd = clone(d0);
      for (int i = 0; i < n_s - 1; ++i) {
        bb[i] += p_v[i];
        dd[i] += p_v[i] * payoff[i + 1];
      }

      // Solve
      NumericVector solution = thomas_algorithm(a0, bb, c0, dd);

      // L-infinity error (economically meaningful)
      double maxdiff = 0.0;
      for (int i = 0; i < n_s - 1; ++i) {
        double diff = std::fabs(solution[i] - temp_prices[i]);
        if (diff > maxdiff) maxdiff = diff;
      }
      error = maxdiff;
      temp_prices = solution;
      ++iter;
    }

    if (iter >= max_iter)
      Rcpp::warning("Penalty iteration did not converge at time step %d (error = %.2e)", n, error);

    // Write interior solution and project to obstacle
    for (int i = 1; i < n_s; ++i)
      new_prices[i] = std::max(temp_prices[i - 1], payoff[i]);

    new_prices[0]   = u_bound[n];
    new_prices[n_s] = l_bound[n];

    old_prices = new_prices;
  }

  // Interpolation at spot
  if (s_0 <= s_min) return old_prices[0];
  if (s_0 >= s_max) return old_prices[n_s];
  return linear_interpolation(s_0, s_min, ds, old_prices);
}

}

//' American Option 1D (local volatility, penalty method)
//'
//' Prices an American option on a single equity with local volatility using
//' Crank-Nicolson finite differences and a penalty-projection method for the
//' early exercise constraint.
//'
//' @param s_0 Stock spot price.
//' @param k Strike price.
//' @param tau Time to expiry (in years).
//' @param r_d Risk-free rate (domestic).
//' @param q Dividend yield.
//' @param sigma Local volatility matrix of size \code{(n_s + 1) x (n_t + 1)},
//'   sampled on the spatial grid (rows) and time grid including both endpoints (columns).
//' @param type Either \code{"call"} or \code{"put"}.
//' @param s_min,s_max Min/Max of the underlying prices grid.
//' @param n_s Number of intervals in the asset grid (\code{n_s + 1} nodes).
//' @param n_t Number of time steps.
//' @param lambda Penalty parameter (> 0, typically 1e4 to 1e6).
//' @param tolerance Convergence tolerance for the penalty iterations.
//'
//' @return Option price as a numeric scalar.
//'
//' @details
//' The solver uses the standard Black-Scholes PDE with time- and state-dependent
//' volatility \eqn{\sigma(S, t)}, discretised via Crank-Nicolson. The early exercise
//' constraint is enforced through a penalty-projection fixed-point iteration at each
//' time step, with a maximum of 200 iterations and a warning on non-convergence.
//'
//' @export
// [[Rcpp::export]]
double american_option_lv(double s_0,
                          double k,
                          double tau,
                          double r_d,
                          double q,
                          NumericMatrix sigma,
                          String type,
                          double s_min,
                          double s_max,
                          int n_s,
                          int n_t,
                          double lambda,
                          double tolerance) {

  return LocalVolatility::american_option_lv(s_0, k, tau, r_d, q, sigma, type,
                                             s_min, s_max, n_s, n_t,
                                             lambda, tolerance);
}
