// [[Rcpp::interfaces(r, cpp)]]

#include <Rcpp.h>
#include "utils.h"

using namespace Rcpp;

namespace LocalVolatility {

  double european_option_lv_2d(double s_0,
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
                               double alpha) {

    // Initialize grid and parameters
    NumericVector s = tavella_randall(s_0, alpha, s_min, s_max, n_s);
    NumericVector x = tavella_randall(x_0, alpha, x_min, x_max, n_x);

    NumericMatrix u(n_s + 1, n_x + 1);

    if(type == "call"){
      for (int i = 0; i < n_s; ++i) {
        for (int j = 0; j < n_x; ++j) {
          u(i, j) = std::max(s[i] * x[j] - k, 0.0);
        }
      }

      // Boundary conditions payoff
      // for (int i = 1; i < n_s + 1; i++) u(i, n_x) = u(i, n_x - 1);
      // for (int j = 1; j < n_x + 1; j++) u(n_s, j) = u(n_s - 1, j);

    }else if(type == "put"){
      for (int i = 0; i < n_s; ++i) {
        for (int j = 0; j < n_x; ++j) {
          u(i, j) = std::max(k - s[i] * x[j], 0.0);
        }
      }

      // Boundary conditions payoff
      // Homogeneous Neumann boundary conditions
      // for (int i = 0; i <= n_s; i++) {
      //   u(i, 0) = u(i, 1);         // Left boundary (x_min)
      //   u(i, n_x) = u(i, n_x - 1); // Right boundary (x_max)
      // }
      // for (int j = 0; j <= n_x; j++) {
      //   u(0, j) = u(1, j);         // Bottom boundary (s_min)
      //   u(n_s, j) = u(n_s - 1, j); // Top boundary (s_max)
      // }
    }

    // Mixed Boundary Condition Example (linear extrapolation)
    for (int i = 0; i <= n_s; i++) {
      u(i, n_x) = 2 * u(i, n_x - 1) - u(i, n_x - 2);  // Right boundary extrapolated
    }
    for (int j = 0; j <= n_x; j++) {
      u(n_s, j) = 2 * u(n_s - 1, j) - u(n_s - 2, j);  // Top boundary extrapolated
    }


    NumericMatrix v = clone(u);


    NumericVector hs(n_s), hx(n_x);
    for (int i = 0; i < n_s-1; i++) hs[i] = s[i + 1] - s[i];
    hs[n_s - 1] = hs[n_s - 2];

    for (int j = 0; j < n_x-1; j++) hx[j] = x[j + 1] - x[j];
    hx[n_x - 1] = hx[n_x - 2];

    double dt = tau / n_t;

    NumericVector as(n_s - 1), bs(n_s - 1), gs(n_s), fs(n_s - 1);
    NumericVector ax(n_x - 1), bx(n_x - 1), gx(n_x), fx(n_x - 1);

    int n = 2 * n_t;

    // Main time loop
    for (int t = 0; t < n_t; ++t) {

      n -= 1;

      // Step in the s-direction
      for (int j = 1; j < n_x; ++j) {
        for (int i = 1; i < n_s; ++i) {

          // Rcpp::Rcout << "n: " << n << std::endl;
          // Rcpp::Rcout << "s: " << s[i] << std::endl;

          as[i - 1] = ((r_d - q) * s[i] * hs[i] - std::pow(sigma_s(i,n) * s[i], 2)) / (hs[i - 1] * (hs[i - 1] + hs[i]));
          bs[i - 1] = 1 / dt + (std::pow(sigma_s(i,n) * s[i], 2) - (r_d - q) * s[i] * (hs[i] - hs[i - 1])) / (hs[i - 1] * hs[i]) + 0.5 * r_d;
          gs[i - 1] = (-(r_d - q) * s[i] * hs[i - 1] - std::pow(sigma_s(i,n) * s[i], 2)) / (hs[i] * (hs[i - 1] + hs[i]));

          fs[i - 1] = u(i, j) / dt + 0.5 * rho * sigma_s(i,n) * sigma_x(j,n) * s[i] * x[j] *
            (u(i + 1, j + 1) + u(i - 1, j - 1) - u(i - 1, j + 1) - u(i + 1, j - 1)) /
              (hs[i - 1] * hx[j] + hs[i] * hx[j] + hs[i] * hx[j - 1] + hs[i - 1] * hx[j - 1]);
        }

        bs[n_s - 2] += gs[n_s - 2];
        as[0] = 0;
        gs[n_s - 1] = 0;

        NumericVector v_slice = thomas_algorithm(as, bs, gs, fs);
        for (int i = 1; i < n_s; i++) v(i, j) = v_slice[i - 1];
      }

      // Boundary conditions for v
      // for (int i = 1; i < n_s + 1; i++) v(i, n_x) = v(i, n_x - 1);
      // for (int j = 1; j < n_x + 1; j++) v(n_s, j) = v(n_s - 1, j);

      // Mixed Boundary Condition Example (linear extrapolation)
      for (int i = 0; i <= n_s; i++) {
        v(i, n_x) = 2 * v(i, n_x - 1) - v(i, n_x - 2);  // Right boundary extrapolated
      }
      for (int j = 0; j <= n_x; j++) {
        v(n_s, j) = 2 * v(n_s - 1, j) - v(n_s - 2, j);  // Top boundary extrapolated
      }

      n -= 1;

      // Rcpp::Rcout << "n: " << n << std::endl;

      // Step in the x-direction
      for (int i = 1; i < n_s; ++i) {
        for (int j = 1; j < n_x; ++j) {
          ax[j - 1] = ((r_d - r_f) * x[j] * hx[j] - std::pow(sigma_x(j,n) * x[j], 2)) / (hx[j - 1] * (hx[j - 1] + hx[j]));
          bx[j - 1] = 1 / dt + (std::pow(sigma_x(j,n) * x[j], 2) - (r_d - r_f) * x[j] * (hx[j] - hx[j - 1])) / (hx[j - 1] * hx[j]) + 0.5 * r_d;
          gx[j - 1] = (-(r_d - r_f) * x[j] * hx[j - 1] - std::pow(sigma_x(j,n) * x[j], 2)) / (hx[j] * (hx[j - 1] + hx[j]));

          fx[j - 1] = v(i, j) / dt + 0.5 * rho * sigma_s(i,n) * sigma_x(j,n) * s[i] * x[j] *
            (v(i + 1, j + 1) + v(i - 1, j - 1) - v(i - 1, j + 1) - v(i + 1, j - 1)) /
              (hx[j - 1] * hs[i] + hx[j] * hs[i] + hx[j] * hs[i - 1] + hx[j - 1] * hs[i - 1]);
        }

        bx[n_x - 2] += gx[n_x - 2];
        ax[0] = 0;
        gx[n_x - 1] = 0;

        NumericVector u_slice = thomas_algorithm(ax, bx, gx, fx);
        for (int j = 1; j < n_x; j++) u(i, j) = u_slice[j - 1];
      }

      // Boundary conditions for u_new
      // for (int i = 1; i < n_s + 1; i++) u(i, n_x) = u(i, n_x - 1);
      // for (int j = 1; j < n_x + 1; j++) u(n_s, j) = u(n_s - 1, j);

      // Mixed Boundary Condition Example (linear extrapolation)
      for (int i = 0; i <= n_s; i++) {
        u(i, n_x) = 2 * u(i, n_x - 1) - u(i, n_x - 2);  // Right boundary extrapolated
      }
      for (int j = 0; j <= n_x; j++) {
        u(n_s, j) = 2 * u(n_s - 1, j) - u(n_s - 2, j);  // Top boundary extrapolated
      }

    }
    // Interpolation step, needs C++ interpolation function equivalent to akima::bilinear
    double result = bilinear_interpolation(s, x, u, s_0, x_0);  // Placeholder for interpolation function
    return result;
  }

}

//' European Option Local Volatility 2D
//'
//' This function evaluates an European-style option on a common stock with local volatility in a foreing currency using finite differences.
//'
//' @param s_0 Stock spot price
//' @param x_0 FX spot price
//' @param k Strike price
//' @param tau Time to expiry
//' @param r_d Risk-free rate (domestic)
//' @param r_f Risk-free rate (foreign)
//' @param q Dividend yield
//' @param sigma_s Local volatility matrix of the stock
//' @param sigma_x Local volatility matrix of the FX
//' @param rho Correlatin between the stock and the FX
//' @param type Either "call" or "put"
//' @param s_min Min value of the stock prices grid
//' @param s_max Man value of the stock prices grid
//' @param x_min Min value of the FX prices grid
//' @param x_max Max value of the FX prices grid
//' @param n_s Size of stock grid for finite difference grid
//' @param n_x Size of FX grid for finite difference grid
//' @param n_t Size of time grid for finite difference grid
//' @param alpha Parameter that defines the uniformity of the grid
//'
//' @export
// [[Rcpp::export]]
double european_option_lv_2d(double s_0,
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
                                    double alpha) {

  return LocalVolatility::european_option_lv_2d(s_0,
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
                                                alpha);
}
