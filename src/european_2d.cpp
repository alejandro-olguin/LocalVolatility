// [[Rcpp::interfaces(r, cpp)]]

#include <Rcpp.h>
#include "utils.h"

using namespace Rcpp;

namespace LocalVolatility {

double european_option_2d(double s_0,
                          double x_0,
                          double k,
                          double tau,
                          double r_d,
                          double r_f,
                          double q,
                          double sigma_s,
                          double sigma_x,
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

  // Input guards
  if (n_s < 3 || n_x < 3) stop("n_s and n_x must be >= 3");
  if (n_t < 1)             stop("n_t must be >= 1");
  if (s_max <= s_min || x_max <= x_min) stop("grid bounds must be increasing");
  if (tau <= 0.0)          stop("tau must be > 0");
  if (rho < -1.0 || rho > 1.0) stop("rho must be within [-1, 1]");
  if (sigma_s < 0.0 || sigma_x < 0.0) stop("volatilities must be >= 0");

  // Grids (Tavella-Randall non-uniform)
  NumericVector s = tavella_randall(s_0, alpha, s_min, s_max, n_s);
  NumericVector x = tavella_randall(x_0, alpha, x_min, x_max, n_x);

  // Non-uniform steps
  NumericVector hs(n_s), hx(n_x);
  for (int i = 0; i < n_s; ++i) hs[i] = s[i + 1] - s[i];
  for (int j = 0; j < n_x; ++j) hx[j] = x[j + 1] - x[j];

  const double dt = tau / n_t;

  // Terminal payoff
  NumericMatrix U(n_s + 1, n_x + 1);
  if (type == "call") {
    for (int i = 0; i <= n_s; ++i)
      for (int j = 0; j <= n_x; ++j)
        U(i, j) = std::max(s[i] * x[j] - k, 0.0);
  } else if (type == "put") {
    for (int i = 0; i <= n_s; ++i)
      for (int j = 0; j <= n_x; ++j)
        U(i, j) = std::max(k - s[i] * x[j], 0.0);
  } else {
    stop("type must be \"call\" or \"put\"");
  }

  // Tridiagonal buffers (S-pass, X-pass)
  NumericVector aS(n_s - 1), bS(n_s - 1), cS(n_s - 1), fS(n_s - 1);
  NumericVector aX(n_x - 1), bX(n_x - 1), cX(n_x - 1), fX(n_x - 1);

  // Analytic far-field Dirichlet boundaries for product options
  // Under domestic measure: forward of S*X is S*X*exp((r_d - r_f - q)*tau)
  auto set_S_boundaries = [&](double tau_rem,
                              NumericVector& leftS, NumericVector& rightS) {
    leftS  = NumericVector(n_x + 1);
    rightS = NumericVector(n_x + 1);
    if (type == "call") {
      for (int j = 0; j <= n_x; ++j) {
        leftS[j]  = 0.0;
        rightS[j] = std::max(s[n_s] * x[j] * std::exp(-(r_f + q) * tau_rem)
                                - k * std::exp(-r_d * tau_rem), 0.0);
      }
    } else {
      for (int j = 0; j <= n_x; ++j) {
        leftS[j]  = k * std::exp(-r_d * tau_rem);
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
        rightX[i] = std::max(s[i] * x[n_x] * std::exp(-(r_f + q) * tau_rem)
                                - k * std::exp(-r_d * tau_rem), 0.0);
      }
    } else {
      for (int i = 0; i <= n_s; ++i) {
        leftX[i]  = k * std::exp(-r_d * tau_rem);
        rightX[i] = 0.0;
      }
    }
  };

  // Domestic-measure drift for S
  const double muS = r_f - q - rho * sigma_s * sigma_x;

  // --- Yanenko operator splitting (backward in time) ---
  for (int n = 0; n < n_t; ++n) {
    const double tau_rem = (n_t - (n + 1)) * dt;

    NumericVector leftS, rightS, leftX, rightX;
    set_S_boundaries(tau_rem, leftS, rightS);
    set_X_boundaries(tau_rem, leftX, rightX);

    // === Pass 1: implicit in S, explicit cross ===
    NumericMatrix V(n_s + 1, n_x + 1);
    for (int j = 1; j < n_x; ++j) {
      const double L = leftS[j];
      const double R = rightS[j];

      for (int i = 1; i < n_s; ++i) {
        const double h_im1 = hs[i - 1], h_i = hs[i];
        const double Si = s[i], Xj = x[j];
        const double a2S = sigma_s * sigma_s * Si * Si;

        aS[i - 1] = (muS * Si * h_i - a2S) / (h_im1 * (h_im1 + h_i));
        bS[i - 1] = 1.0 / dt
          + (a2S - muS * Si * (h_i - h_im1)) / (h_im1 * h_i)
          + 0.5 * r_d;
        cS[i - 1] = (-muS * Si * h_im1 - a2S) / (h_i * (h_im1 + h_i));

        // Explicit cross term (centered mixed derivative)
        const double denom = (h_im1 + h_i) * (hx[j - 1] + hx[j]);
        double cross = 0.0;
        if (denom != 0.0) {
          cross = 0.5 * rho * sigma_s * sigma_x * Si * Xj *
            (U(i+1,j+1) + U(i-1,j-1) - U(i-1,j+1) - U(i+1,j-1)) / denom;
        }

        fS[i - 1] = U(i, j) / dt + cross;
      }

      fS[0]       -= aS[0]      * L;
      fS[n_s - 2] -= cS[n_s - 2] * R;
      aS[0] = 0.0; cS[n_s - 2] = 0.0;

      NumericVector col = thomas_algorithm(aS, bS, cS, fS);
      V(0, j) = L; V(n_s, j) = R;
      for (int i = 1; i < n_s; ++i) V(i, j) = col[i - 1];
    }

    for (int i = 0; i <= n_s; ++i) { V(i, 0) = leftX[i]; V(i, n_x) = rightX[i]; }
    for (int j = 0; j <= n_x; ++j) { V(0, j) = leftS[j]; V(n_s, j) = rightS[j]; }

    // === Pass 2: implicit in X, explicit cross ===
    NumericMatrix Unew(n_s + 1, n_x + 1);
    for (int i = 1; i < n_s; ++i) {
      const double L = leftX[i];
      const double R = rightX[i];

      for (int j = 1; j < n_x; ++j) {
        const double h_jm1 = hx[j - 1], h_j = hx[j];
        const double Si = s[i], Xj = x[j];
        const double a2X = sigma_x * sigma_x * Xj * Xj;

        aX[j - 1] = ((r_d - r_f) * Xj * h_j - a2X) / (h_jm1 * (h_jm1 + h_j));
        bX[j - 1] = 1.0 / dt
          + (a2X - (r_d - r_f) * Xj * (h_j - h_jm1)) / (h_jm1 * h_j)
          + 0.5 * r_d;
        cX[j - 1] = (-(r_d - r_f) * Xj * h_jm1 - a2X) / (h_j * (h_jm1 + h_j));

        // Explicit cross term using V
        const double denom = (hs[i - 1] + hs[i]) * (h_jm1 + h_j);
        double cross = 0.0;
        if (denom != 0.0) {
          cross = 0.5 * rho * sigma_s * sigma_x * Si * Xj *
            (V(i+1,j+1) + V(i-1,j-1) - V(i-1,j+1) - V(i+1,j-1)) / denom;
        }

        fX[j - 1] = V(i, j) / dt + cross;
      }

      fX[0]       -= aX[0]      * L;
      fX[n_x - 2] -= cX[n_x - 2] * R;
      aX[0] = 0.0; cX[n_x - 2] = 0.0;

      NumericVector row = thomas_algorithm(aX, bX, cX, fX);
      Unew(i, 0) = L; Unew(i, n_x) = R;
      for (int j = 1; j < n_x; ++j) Unew(i, j) = row[j - 1];
    }

    for (int i = 0; i <= n_s; ++i) { Unew(i, 0) = leftX[i]; Unew(i, n_x) = rightX[i]; }
    for (int j = 0; j <= n_x; ++j) { Unew(0, j) = leftS[j]; Unew(n_s, j) = rightS[j]; }

    std::swap(U, Unew);
  }

  return bilinear_interpolation(s, x, U, s_0, x_0);
}

}

//' European Option 2D (constant volatility)
//'
//' Prices a European option on an equity quoted in a foreign currency using a
//' 2D finite-difference PDE with Yanenko operator splitting, constant
//' volatilities, and a centered mixed-derivative stencil on a non-uniform
//' Tavella-Randall grid.
//'
//' @param s_0 Stock spot price.
//' @param x_0 FX spot price (domestic per foreign).
//' @param k Strike price (domestic currency).
//' @param tau Time to expiry (in years).
//' @param r_d Risk-free rate (domestic).
//' @param r_f Risk-free rate (foreign).
//' @param q Dividend yield.
//' @param sigma_s Constant volatility of the stock (>= 0).
//' @param sigma_x Constant volatility of the FX (>= 0).
//' @param rho Correlation between the stock and the FX in \code{[-1, 1]}.
//' @param type Either \code{"call"} or \code{"put"}.
//' @param s_min,s_max Min/Max of the stock grid.
//' @param x_min,x_max Min/Max of the FX grid.
//' @param n_s Number of intervals in stock grid (\code{n_s + 1} nodes).
//' @param n_x Number of intervals in FX grid (\code{n_x + 1} nodes).
//' @param n_t Number of time steps.
//' @param alpha Grid clustering parameter for Tavella-Randall grids (> 0).
//'
//' @return Option price as a numeric scalar.
//'
//' @details
//' Domestic-measure dynamics: drift \eqn{\mu_S = r_f - q - \rho \sigma_S \sigma_X}
//' in S, \eqn{r_d - r_f} in X, discounting at \eqn{r_d}. Far-field Dirichlet
//' boundaries use \eqn{e^{-(r_f + q)\tau}} on the \eqn{S \cdot X} leg and
//' \eqn{e^{-r_d \tau}} on \eqn{K}.
//'
//' @examples
//' european_option_2d(100, 20, 2000, 1, 0.05, 0.02, 0.01,
//'                    0.2, 0.15, 0.3, "call",
//'                    10, 300, 1, 60, 20, 20, 25, 3)
//'
//' @export
// [[Rcpp::export]]
double european_option_2d(double s_0,
                          double x_0,
                          double k,
                          double tau,
                          double r_d,
                          double r_f,
                          double q,
                          double sigma_s,
                          double sigma_x,
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

  return LocalVolatility::european_option_2d(s_0, x_0, k, tau, r_d, r_f, q,
                                             sigma_s, sigma_x, rho, type,
                                             s_min, s_max, x_min, x_max,
                                             n_s, n_x, n_t, alpha);
}
