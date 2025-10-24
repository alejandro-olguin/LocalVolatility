// [[Rcpp::interfaces(r, cpp)]]

#include <Rcpp.h>
#include "utils.h"

using namespace Rcpp;

//' European Option 2D (constant volatility)
//'
//' Prices a European option on an equity quoted in a foreign currency using a 2D finite-difference PDE
//' with Yanenko operator splitting, constant volatilities (sigma_s, sigma_x), and a centered mixed-derivative.
//'
//' @param s_0 Stock spot price.
//' @param x_0 FX spot price (domestic per foreign).
//' @param k Strike price (domestic currency).
//' @param tau Time to expiry (in years).
//' @param r_d Risk-free rate (domestic).
//' @param r_f Risk-free rate (foreign).
//' @param q Dividend yield.
//' @param sigma_s Constant volatility of the stock.
//' @param sigma_x Constant volatility of the FX.
//' @param rho Correlation between the stock and the FX in `[-1, 1]`.
//' @param type Either "call" or "put".
//' @param s_min,s_max Min/Max of the stock grid.
//' @param x_min,x_max Min/Max of the FX grid.
//' @param n_s Number of intervals in stock grid (n_s + 1 nodes).
//' @param n_x Number of intervals in FX grid (n_x + 1 nodes).
//' @param n_t Number of time steps.
//' @param alpha Grid clustering parameter for Tavella–Randall grids (> 0).
//'
//' @details
//' Domestic-measure drift and discounting are applied. In the S-direction, the drift is
//' mu_S = r_f - q - rho * sigma_s * sigma_x. In the X-direction, the drift is r_d - r_f, and
//' discounting is at r_d. Far-field Dirichlet boundaries use exp(-q * tau) on the S·X term and
//' exp(-r_d * tau) on K. The spot (s_0, x_0) price is obtained via bilinear interpolation.
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
  // --- guards (don’t change API) ---
  if (n_s < 3 || n_x < 3) stop("n_s and n_x must be >= 3");
  if (n_t < 1)            stop("n_t must be >= 1");
  if (s_max <= s_min || x_max <= x_min) stop("grid bounds must be increasing");
  if (tau <= 0.0)         stop("tau must be > 0");

  // Grids (Tavella–Randall non-uniform, size n_s+1 / n_x+1)
  NumericVector s = tavella_randall(s_0, alpha, s_min, s_max, n_s);
  NumericVector x = tavella_randall(x_0, alpha, x_min, x_max, n_x);

  // Non-uniform steps
  NumericVector hs(n_s), hx(n_x);
  for (int i = 0; i < n_s; ++i) hs[i] = s[i + 1] - s[i];
  for (int j = 0; j < n_x; ++j) hx[j] = x[j + 1] - x[j];

  const double dt = tau / n_t;

  // Terminal payoff U(T)
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

  // Helpers: analytic far-field Dirichlet boundaries for product options
  auto set_S_boundaries = [&](double tau_rem,
                              NumericVector& leftS, NumericVector& rightS) {
    leftS  = NumericVector(n_x + 1);
    rightS = NumericVector(n_x + 1);
    if (type == "call") {
      for (int j = 0; j <= n_x; ++j) {
        leftS[j]  = 0.0; // S -> 0
        double up = std::max(s[n_s] * x[j] * std::exp(-q * tau_rem)
                               - k * std::exp(-r_d * tau_rem), 0.0);
        rightS[j] = up;  // S -> ∞
      }
    } else { // put
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
        leftX[i]  = 0.0; // X -> 0
        double up = std::max(s[i] * x[n_x] * std::exp(-q * tau_rem)
                               - k * std::exp(-r_d * tau_rem), 0.0);
        rightX[i] = up;  // X -> ∞
      }
    } else { // put
      for (int i = 0; i <= n_s; ++i) {
        leftX[i]  = k * std::exp(-r_d * tau_rem);
        rightX[i] = 0.0;
      }
    }
  };

  // --- Yanenko operator splitting (backward in time) ---
  for (int n = 0; n < n_t; ++n) {
    double tau_rem = (n_t - (n + 1)) * dt;

    // Boundaries for this time level
    NumericVector leftS, rightS, leftX, rightX;
    set_S_boundaries(tau_rem, leftS, rightS);
    set_X_boundaries(tau_rem, leftX, rightX);

    // ===== Pass 1: implicit in S (for each fixed j), explicit in X and cross =====
    NumericMatrix V(n_s + 1, n_x + 1);
    for (int j = 1; j < n_x; ++j) {
      double L = leftS[j];   // U at i=0, col j
      double R = rightS[j];  // U at i=n_s, col j

      for (int i = 1; i < n_s; ++i) {
        double h_im1 = hs[i - 1], h_i = hs[i];
        double Si = s[i], Xj = x[j];

        // S-direction coefficients (constant vol in S)
        double a2S = std::pow(sigma_s * Si, 2.0);

        // *** domestic-measure drift for S: muS = r_f - q - rho*sigma_s*sigma_x ***
        double muS = r_f - q - rho * sigma_s * sigma_x;

        aS[i - 1] = ( (muS * Si * h_i) - a2S ) / (h_im1 * (h_im1 + h_i));
        bS[i - 1] = 1.0 / dt
        + ( a2S - muS * Si * (h_i - h_im1) ) / (h_im1 * h_i)
          + 0.5 * r_d;
          cS[i - 1] = ( -(muS * Si * h_im1) - a2S ) / (h_i * (h_im1 + h_i));

          // Explicit cross term (centered mixed derivative using U)
          double denom = (h_im1 * hx[j - 1] + h_i * hx[j] + h_i * hx[j - 1] + h_im1 * hx[j]);
          double cross = 0.0;
          if (denom != 0.0) {
            cross = 0.5 * rho * sigma_s * sigma_x * Si * Xj *
              (U(i + 1, j + 1) + U(i - 1, j - 1) - U(i - 1, j + 1) - U(i + 1, j - 1)) / denom;
          }

          fS[i - 1] = U(i, j) / dt + cross;
      }

      // Boundary contributions to RHS, then zero the touching off-diagonals
      fS[0]         -= aS[0]       * L;
      fS[n_s - 2]   -= cS[n_s - 2] * R;
      aS[0]          = 0.0;
      cS[n_s - 2]    = 0.0;

      // Solve in S
      NumericVector col = thomas_algorithm(aS, bS, cS, fS);
      V(0, j)    = L;
      V(n_s, j)  = R;
      for (int i = 1; i < n_s; ++i) V(i, j) = col[i - 1];
    }

    // Enforce frame boundaries on V
    for (int i = 0; i <= n_s; ++i) { V(i, 0) = leftX[i];  V(i, n_x) = rightX[i]; }
    for (int j = 0; j <= n_x; ++j) { V(0, j) = leftS[j];  V(n_s, j) = rightS[j]; }

    // ===== Pass 2: implicit in X (for each fixed i), explicit in S and cross =====
    NumericMatrix Unew(n_s + 1, n_x + 1);
    for (int i = 1; i < n_s; ++i) {
      double L = leftX[i];   // V at j=0, row i
      double R = rightX[i];  // V at j=n_x, row i

      for (int j = 1; j < n_x; ++j) {
        double h_jm1 = hx[j - 1], h_j = hx[j];
        double Si = s[i], Xj = x[j];

        // X-direction coefficients (constant vol in X)
        double a2X = std::pow(sigma_x * Xj, 2.0);
        aX[j - 1] = ((r_d - r_f) * Xj * h_j - a2X) / (h_jm1 * (h_jm1 + h_j));
        bX[j - 1] = 1.0 / dt
        + (a2X - (r_d - r_f) * Xj * (h_j - h_jm1)) / (h_jm1 * h_j)
          + 0.5 * r_d;
          cX[j - 1] = (-(r_d - r_f) * Xj * h_jm1 - a2X) / (h_j * (h_jm1 + h_j));

          // Explicit cross using V
          double denom = (h_jm1 * hs[i - 1] + h_j * hs[i] + h_jm1 * hs[i - 1] + h_j * hs[i]);
          double cross = 0.0;
          if (denom != 0.0) {
            cross = 0.5 * rho * sigma_s * sigma_x * Si * Xj *
              (V(i + 1, j + 1) + V(i - 1, j - 1) - V(i - 1, j + 1) - V(i + 1, j - 1)) / denom;
          }

          fX[j - 1] = V(i, j) / dt + cross;
      }

      // Boundary contributions to RHS in X, then zero
      fX[0]         -= aX[0]       * L;
      fX[n_x - 2]   -= cX[n_x - 2] * R;
      aX[0]          = 0.0;
      cX[n_x - 2]    = 0.0;

      // Solve in X
      NumericVector row = thomas_algorithm(aX, bX, cX, fX);
      Unew(i, 0)   = L;
      Unew(i, n_x) = R;
      for (int j = 1; j < n_x; ++j) Unew(i, j) = row[j - 1];
    }

    // Enforce frame boundaries on Unew and advance
    for (int i = 0; i <= n_s; ++i) { Unew(i, 0) = leftX[i];  Unew(i, n_x) = rightX[i]; }
    for (int j = 0; j <= n_x; ++j) { Unew(0, j) = leftS[j];  Unew(n_s, j) = rightS[j]; }

    U = clone(Unew);
  }

  // Bilinear interpolation at (s_0, x_0)
  return bilinear_interpolation(s, x, U, s_0, x_0);
}
