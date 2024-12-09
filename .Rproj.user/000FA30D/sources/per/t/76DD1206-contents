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

    double ds = (s_max - s_min)/n_s;
    double dt = tau/n_t;

    // create storage for the option prices
    NumericVector old_prices(n_s+1), new_prices(n_s+1), payoff(n_s+1), u_bound(n_t+1), l_bound(n_t+1);

    if(type == "call"){

      // Payoff boundary condition
      for(int i = 0; i <= n_s; i++){
        payoff[i] = std::max(s_min + i*ds - k, 0.0); // Call payoff
      }

      // Upper and lower boundary condition
      for(int n = n_t; n >= 0; n--){
        u_bound[n] = 0;
        l_bound[n] = std::max(s_max - k, 0.0);
      }

    }else if(type == "put"){

      // Payoff boundary condition
      for(int i = 0; i <= n_s; i++){
        payoff[i] = std::max(k - s_min - i*ds, 0.0); // Put payoff
      }

      // Upper and lower boundary condition
      for(int n = n_t; n >= 0; n--){
        u_bound[n] = k * std::exp(-r_d * (n_t-n) * dt);
        l_bound[n] = 0;
      }
    }

    old_prices = payoff;

    // declare vectors for matrix equations
    NumericVector a(n_s-1), b(n_s-1), c(n_s-1), d(n_s-1), temp_prices(n_s-1);

    for(int n = n_t-1; n >= 0; n--){

      new_prices[0] = u_bound[n];
      new_prices[n_s] = l_bound[n];

      temp_prices = old_prices[Rcpp::Range(1,n_s-1)];

      for(int i = 1; i < n_s; i++){

        a[i-1] = -(dt/4 * (std::pow(sigma(i,n) * (s_min + i*ds)/ds, 2.) - (r_d-q)*(s_min + i*ds)/ds));
        b[i-1] = 1 + (dt/2 * (std::pow(sigma(i,n) * (s_min + i*ds)/ds, 2.) + r_d));
        c[i-1] = -(dt/4 * (std::pow(sigma(i,n) * (s_min + i*ds)/ds, 2.) + (r_d-q)*(s_min + i*ds)/ds));

        d[i-1] = (dt/4 * (std::pow(sigma(i,n+1) * (s_min + i*ds)/ds, 2.) - (r_d-q)*(s_min + i*ds)/ds)) * old_prices[i-1] +
          (1 - (dt/2 * (std::pow(sigma(i,n+1) * (s_min + i*ds)/ds, 2.) + r_d))) * old_prices[i] +
          (dt/4 * (std::pow(sigma(i,n+1) * (s_min + i*ds)/ds, 2.) + (r_d-q)*(s_min + i*ds)/ds)) * old_prices[i+1];
      }

      // Adding boundary conditions
      d[0]   = d[0]   - a[0]   * u_bound[n];
      d[n_s-2] = d[n_s-2] - c[n_s-2] * l_bound[n];

      // Resetting already used diagonal values
      a[0]   = 0;
      c[n_s-2] = 0;

      // To enter the loop
      double error  = 100;

      while(error > tolerance){

        NumericVector p_v = ifelse(temp_prices < payoff[Rcpp::Range(1, n_s-1)], lambda, 0.0);

        b += p_v;
        d += p_v * payoff[Rcpp::Range(1, n_s-1)];

        // Solve Equation With Thomas Algorithm
        NumericVector solution = thomas_algorithm(a, b, c, d);

        // Compute error
        error = sum(pow(temp_prices - solution, 2.0));

        // Update New Prices
        temp_prices = solution;
      }

      // Update New Prices
      new_prices[Range(1,n_s-1)] = temp_prices;

      // Set old prices as the new prices for next iteration
      old_prices = new_prices;
    }

    double price = linear_interpolation(s_0, ds, old_prices);

    return(price);
  }
}

//' American Option Local Volatility 1D
//'
//' This function evaluates an American-style option on a common stock with local volatility using finite differences.
//'
//' @param s_0 Stock spot price
//' @param k Strike price
//' @param tau Time to expiry
//' @param r_d Risk-free rate (domestic)
//' @param q Dividend yield
//' @param sigma Local volatility matrix
//' @param type Either "call" or "put"
//' @param s_min max value of the underlying prices grid
//' @param s_max min value of the underlying prices grid
//' @param n_s Size of asset grid for finite difference grid
//' @param n_t Size of time grid for finite difference grid
//' @param lambda Penalty parameter greater than 1
//' @param tolerance Error tolerance
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
