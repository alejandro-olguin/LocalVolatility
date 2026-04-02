// [[Rcpp::interfaces(r, cpp)]]

#include <Rcpp.h>
#include "utils.h"

using namespace Rcpp;

namespace LocalVolatility {

double european_option_lv(double s_0,
                          double k,
                          double tau,
                          double r_d,
                          double q,
                          NumericMatrix sigma,
                          String type,
                          double s_min,
                          double s_max,
                          int n_s,
                          int n_t) {

  // Input guards
  if (n_s < 2 || n_t < 1) stop("n_s >= 2 and n_t >= 1 required");
  if (s_max <= s_min)     stop("s_max must be > s_min");
  if (tau <= 0.0)         stop("tau must be > 0");
  if (sigma.nrow() != n_s + 1 || sigma.ncol() != n_t + 1)
    stop("sigma must have dimensions (n_s+1) x (n_t+1)");

  const double ds = (s_max - s_min) / n_s;
  const double dt = tau / n_t;

  // Storage
  NumericVector old_prices(n_s + 1), new_prices(n_s + 1),
    payoff(n_s + 1), u_bound(n_t + 1), l_bound(n_t + 1);

  // Payoff and time-dependent boundaries
  if (type == "call") {
    for (int i = 0; i <= n_s; ++i)
      payoff[i] = std::max(s_min + i * ds - k, 0.0);

    for (int n = 0; n <= n_t; ++n) {
      const double t_rem = (n_t - n) * dt;
      u_bound[n] = 0.0;  // S -> 0
      l_bound[n] = std::max(s_max * std::exp(-q * t_rem)
                              - k * std::exp(-r_d * t_rem), 0.0);
    }

  } else if (type == "put") {
    for (int i = 0; i <= n_s; ++i)
      payoff[i] = std::max(k - (s_min + i * ds), 0.0);

    for (int n = 0; n <= n_t; ++n) {
      const double t_rem = (n_t - n) * dt;
      u_bound[n] = k * std::exp(-r_d * t_rem);  // S -> 0
      l_bound[n] = 0.0;                         // S -> +inf
    }

  } else {
    stop("type must be \"call\" or \"put\"");
  }

  // Terminal condition
  old_prices = payoff;

  // Tridiagonal system
  NumericVector a(n_s - 1), b(n_s - 1), c(n_s - 1), d(n_s - 1);

  for (int n = n_t - 1; n >= 0; --n) {
    new_prices[0]   = u_bound[n];
    new_prices[n_s] = l_bound[n];

    for (int i = 1; i < n_s; ++i) {
      const double S_i = s_min + i * ds;

      // Local vols at time levels n and n+1
      const double sig_n   = sigma(i, n);
      const double sig_np1 = sigma(i, n + 1);

      // alpha := 0.5 * sigma^2 * S^2
      const double alpha_n  = 0.5 * sig_n * sig_n * S_i * S_i;
      const double alpha_np = 0.5 * sig_np1 * sig_np1 * S_i * S_i;

      // CN coefficients (implicit part at t_n)
      a[i - 1] = -0.5 * dt * (alpha_n / (ds * ds) - (r_d - q) * S_i / (2.0 * ds));
      b[i - 1] =  1.0      + 0.5 * dt * (2.0 * alpha_n / (ds * ds) + r_d);
      c[i - 1] = -0.5 * dt * (alpha_n / (ds * ds) + (r_d - q) * S_i / (2.0 * ds));

      // RHS from explicit part at t_{n+1}
      d[i - 1] =  0.5 * dt * (alpha_np / (ds * ds) - (r_d - q) * S_i / (2.0 * ds)) * old_prices[i - 1]
        + (1.0 - 0.5 * dt * (2.0 * alpha_np / (ds * ds) + r_d))                     * old_prices[i]
        + 0.5 * dt * (alpha_np / (ds * ds) + (r_d - q) * S_i / (2.0 * ds))          * old_prices[i + 1];
    }

    // Boundary contributions to RHS
    d[0]       -= a[0]      * u_bound[n];
    d[n_s - 2] -= c[n_s - 2] * l_bound[n];

    // Zero the touching off-diagonals before solving
    a[0]       = 0.0;
    c[n_s - 2] = 0.0;

    NumericVector solution = thomas_algorithm(a, b, c, d);

    for (int i = 1; i < n_s; ++i)
      new_prices[i] = solution[i - 1];

    old_prices = new_prices;
  }

  // Interpolation at spot
  if (s_0 <= s_min) return old_prices[0];
  if (s_0 >= s_max) return old_prices[n_s];
  return linear_interpolation(s_0, s_min, ds, old_prices);
}

}

//' European Option 1D (local volatility)
//'
//' Prices a European option on a single equity with local volatility using a
//' 1D Crank-Nicolson finite-difference scheme, time-dependent Dirichlet
//' boundaries, and linear interpolation at the spot.
//'
//' @param s_0 Stock spot price.
//' @param k Strike price.
//' @param tau Time to expiry (in years).
//' @param r_d Risk-free rate (domestic).
//' @param q Dividend yield.
//' @param sigma Local volatility matrix of size \code{(n_s + 1) x (n_t + 1)},
//'   sampled on the spatial grid (rows) and time grid including both endpoints (columns).
//' @param type Either \code{"call"} or \code{"put"}.
//' @param s_min,s_max Min/Max of the asset grid.
//' @param n_s Number of intervals in the asset grid (\code{n_s + 1} nodes).
//' @param n_t Number of time steps.
//'
//' @return Option price as a numeric scalar.
//'
//' @export
// [[Rcpp::export]]
double european_option_lv(double s_0,
                          double k,
                          double tau,
                          double r_d,
                          double q,
                          NumericMatrix sigma,
                          String type,
                          double s_min,
                          double s_max,
                          int n_s,
                          int n_t) {

  return LocalVolatility::european_option_lv(s_0, k, tau, r_d, q, sigma, type,
                                             s_min, s_max, n_s, n_t);
}
