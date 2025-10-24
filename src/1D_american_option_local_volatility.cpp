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

  // pequeñas salvaguardas
  if (n_s < 2 || n_t < 1) stop("n_s >= 2 y n_t >= 1 requeridos");
  if (s_max <= s_min)     stop("s_max debe ser > s_min");
  if (tau <= 0.0)         stop("tau debe ser > 0");

  const double ds = (s_max - s_min)/n_s;
  const double dt = tau/n_t;

  // almacenamiento
  NumericVector old_prices(n_s+1), new_prices(n_s+1),
  payoff(n_s+1), u_bound(n_t+1), l_bound(n_t+1);

  // payoff + fronteras (American)
  if (type == "call") {
    for (int i = 0; i <= n_s; ++i)
      payoff[i] = std::max(s_min + i*ds - k, 0.0);

    // American call: S->0 -> 0 ; S->∞ -> S - K (no descontado)
    for (int n = 0; n <= n_t; ++n) {
      u_bound[n] = 0.0;
      l_bound[n] = std::max(s_max - k, 0.0);
    }

  } else if (type == "put") {
    for (int i = 0; i <= n_s; ++i)
      payoff[i] = std::max(k - (s_min + i*ds), 0.0); // fix

    // American put: S->0 -> K e^{-r_d τ_rem}; S->∞ -> 0
    for (int n = 0; n <= n_t; ++n) {
      const double t_rem = (n_t - n) * dt;
      u_bound[n] = k * std::exp(-r_d * t_rem); // borde inferior (S≈0)
      l_bound[n] = 0.0;                        // borde superior (S≈∞)
    }

  } else {
    stop("type debe ser \"call\" o \"put\"");
  }

  // condición terminal
  old_prices = payoff;

  // tridiagonales
  NumericVector a(n_s-1), b(n_s-1), c(n_s-1), d(n_s-1), temp_prices(n_s-1);

  for (int n = n_t - 1; n >= 0; --n) {

    new_prices[0]   = u_bound[n];
    new_prices[n_s] = l_bound[n];

    // construir sistema base (sin penalización)
    for (int i = 1; i < n_s; ++i) {
      const double S_i    = s_min + i*ds;

      const double sig_n   = sigma(i, n);
      const double sig_np1 = sigma(i, n+1);

      const double alpha_n  = 0.5 * sig_n*sig_n    * S_i*S_i;
      const double alpha_np = 0.5 * sig_np1*sig_np1 * S_i*S_i;

      // parte implícita (t_n)
      a[i-1] = -0.5*dt*( alpha_n/(ds*ds) - (r_d - q)*S_i/(2.0*ds) );
      b[i-1] =  1.0   + 0.5*dt*( 2.0*alpha_n/(ds*ds) + r_d );
      c[i-1] = -0.5*dt*( alpha_n/(ds*ds) + (r_d - q)*S_i/(2.0*ds) );

      // RHS (t_{n+1})
      d[i-1] =  0.5*dt*( alpha_np/(ds*ds) - (r_d - q)*S_i/(2.0*ds) ) * old_prices[i-1]
      + (1.0 - 0.5*dt*( 2.0*alpha_np/(ds*ds) + r_d ))        * old_prices[i]
      + 0.5*dt*( alpha_np/(ds*ds) + (r_d - q)*S_i/(2.0*ds) ) * old_prices[i+1];
    }

    // aportar fronteras al RHS (clave en CN) y anular coef pegados a borde
    d[0]       -= a[0]     * u_bound[n];
    d[n_s - 2] -= c[n_s-2] * l_bound[n];
    a[0]        = 0.0;
    c[n_s - 2]  = 0.0;

    // inicialización de la iteración de penalización
    temp_prices = old_prices[Rcpp::Range(1, n_s-1)];

    // Copias base (no acumular penalización sobre b/d)
    const NumericVector a0 = clone(a);
    const NumericVector b0 = clone(b);
    const NumericVector c0 = clone(c);
    const NumericVector d0 = clone(d);

    double error = std::numeric_limits<double>::infinity();

    while (error > tolerance) {
      // indicador de violación del obstáculo (interior)
      NumericVector p_v(n_s-1);
      for (int i = 0; i < n_s-1; ++i) {
        const double pay_i = payoff[i+1];
        p_v[i] = (temp_prices[i] < pay_i) ? lambda : 0.0;
      }

      // reconstruir sistema para esta iteración (NO acumular)
      NumericVector bb = clone(b0);
      NumericVector dd = clone(d0);
      for (int i = 0; i < n_s-1; ++i) {
        const double pay_i = payoff[i+1];
        bb[i] += p_v[i];
        dd[i] += p_v[i] * pay_i;
      }

      // resolver
      NumericVector solution = thomas_algorithm(a0, bb, c0, dd);

      // error como norma infinito (significativa económicamente)
      double maxdiff = 0.0;
      for (int i = 0; i < n_s-1; ++i) {
        double diff = std::fabs(solution[i] - temp_prices[i]);
        if (diff > maxdiff) maxdiff = diff;
      }
      error = maxdiff;

      temp_prices = solution;
    }

    // escribir solución interior + proyectar al obstáculo
    for (int i = 1; i < n_s; ++i)
      new_prices[i] = std::max(temp_prices[i-1], payoff[i]);

    // bordes
    new_prices[0]   = u_bound[n];
    new_prices[n_s] = l_bound[n];

    // avanzar
    old_prices = new_prices;
  }

  // interpolación segura
  if (s_0 <= s_min) return old_prices[0];
  if (s_0 >= s_max) return old_prices[n_s];
  return linear_interpolation(s_0, ds, old_prices);
}

}

//' American Option 1D (local volatility, penalty method)
//'
//' Prices an American option on a single equity with local volatility using operator splitting
//' and a penalty-projection method for the early exercise constraint.
//'
//' @param s_0 Stock spot price.
//' @param k Strike price.
//' @param tau Time to expiry (in years).
//' @param r_d Risk-free rate (domestic).
//' @param q Dividend yield.
//' @param sigma Local volatility matrix sampled on the grid.
//' @param type Either "call" or "put".
//' @param s_min,s_max Min/Max of the underlying prices grid.
//' @param n_s Number of intervals in the asset grid (n_s + 1 nodes).
//' @param n_t Number of time steps.
//' @param lambda Penalty parameter (> 1).
//' @param tolerance Relative error tolerance for the penalty iterations.
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

  return LocalVolatility::american_option_lv(s_0,
                                             k,
                                             tau,
                                             r_d,
                                             q,
                                             sigma,
                                             type,
                                             s_min,
                                             s_max,
                                             n_s,
                                             n_t,
                                             lambda,
                                             tolerance);

}
