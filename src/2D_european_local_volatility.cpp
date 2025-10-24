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

  // --- guards ---
  if (n_s < 3 || n_x < 3) stop("n_s and n_x must be >= 3");
  if (n_t < 1)            stop("n_t must be >= 1");
  if (s_max <= s_min || x_max <= x_min) stop("grid bounds must be increasing");
  if (tau <= 0.0)         stop("tau must be > 0");

  // Grids (Tavella–Randall, non-uniform)
  NumericVector s = tavella_randall(s_0, alpha, s_min, s_max, n_s);
  NumericVector x = tavella_randall(x_0, alpha, x_min, x_max, n_x);

  // Non-uniform steps
  NumericVector hs(n_s), hx(n_x);
  for (int i = 0; i < n_s; ++i) hs[i] = s[i + 1] - s[i];
  for (int j = 0; j < n_x; ++j) hx[j] = x[j + 1] - x[j];

  const double dt = tau / n_t;

  // Terminal payoff V(T) = max(±(S·X − K), 0)
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

  // Tridiagonal buffers (S- then X-pass)
  NumericVector aS(n_s - 1), bS(n_s - 1), cS(n_s - 1), fS(n_s - 1);
  NumericVector aX(n_x - 1), bX(n_x - 1), cX(n_x - 1), fX(n_x - 1);

  // Helpers: analytic Dirichlet far-field boundaries for product options
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
        leftS[j]  = k * std::exp(-r_d * tau_rem); // S -> 0
        rightS[j] = 0.0;                          // S -> ∞
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
        leftX[i]  = k * std::exp(-r_d * tau_rem); // X -> 0
        rightX[i] = 0.0;                          // X -> ∞
      }
    }
  };

  // --- Yanenko operator splitting over n_t steps (backward in time) ---
  for (int n = 0; n < n_t; ++n) {
    // time remaining AFTER this full step (for boundaries)
    double tau_rem = (n_t - (n + 1)) * dt;

    // Precompute Dirichlet boundaries for this step
    NumericVector leftS, rightS, leftX, rightX;
    set_S_boundaries(tau_rem, leftS, rightS);
    set_X_boundaries(tau_rem, leftX, rightX);

    // ===== Pass 1: implicit in S (per fixed j), explicit in X and cross =====
    NumericMatrix V(n_s + 1, n_x + 1);
    for (int j = 1; j < n_x; ++j) {
      // Known S-boundary values for this column j
      double L = leftS[j];
      double R = rightS[j];

      for (int i = 1; i < n_s; ++i) {
        double h_im1 = hs[i - 1], h_i = hs[i];
        double Si = s[i], Xj = x[j];

        double sigS = sigma_s(i, n);
        double sigX = sigma_x(j, n);

        // *** NEW: domestic-measure drift of S ***
        // mu_S = r_f - q - rho * sigma_S * sigma_X   (pointwise)
        double muS = r_f - q - rho * sigS * sigX;

        // S-direction coefficients on non-uniform mesh (replace (r_d - q) with muS)
        aS[i - 1] = ( (muS * Si * h_i) - std::pow(sigS * Si, 2.0) ) / (h_im1 * (h_im1 + h_i));

        bS[i - 1] = 1.0 / dt
          + ( std::pow(sigS * Si, 2.0) - muS * Si * (h_i - h_im1) ) / (h_im1 * h_i)
          + 0.5 * r_d;

        cS[i - 1] = ( -(muS * Si * h_im1) - std::pow(sigS * Si, 2.0) ) / (h_i * (h_im1 + h_i));

        // Explicit cross term (unchanged)
        double denom_cross = (h_im1 * hx[j - 1] + h_i * hx[j] + h_i * hx[j - 1] + h_im1 * hx[j]);
        double cross = 0.0;
        if (denom_cross != 0.0) {
          cross = 0.5 * rho * sigS * sigX * Si * Xj *
            (U(i + 1, j + 1) + U(i - 1, j - 1) - U(i - 1, j + 1) - U(i + 1, j - 1)) / denom_cross;
        }

        fS[i - 1] = U(i, j) / dt + cross;
      }

      // Boundary contributions… (unchanged)
      fS[0]       -= aS[0]       * L;
      fS[n_s - 2] -= cS[n_s - 2] * R;
      aS[0] = 0.0; cS[n_s - 2] = 0.0;

      // Solve in S (unchanged)
      NumericVector col = thomas_algorithm(aS, bS, cS, fS);
      V(0, j)   = L;
      V(n_s, j) = R;
      for (int i = 1; i < n_s; ++i) V(i, j) = col[i - 1];
    }

    // Enforce outer-frame boundaries on V
    for (int i = 0; i <= n_s; ++i) { V(i, 0) = leftX[i];  V(i, n_x) = rightX[i]; }
    for (int j = 0; j <= n_x; ++j) { V(0, j) = leftS[j];  V(n_s, j) = rightS[j]; }

    // ===== Pass 2: implicit in X (per fixed i), explicit in S and cross =====
    NumericMatrix Unew(n_s + 1, n_x + 1);

    for (int i = 1; i < n_s; ++i) {
      double L = leftX[i];
      double R = rightX[i];

      for (int j = 1; j < n_x; ++j) {
        double h_jm1 = hx[j - 1], h_j = hx[j];
        double Si = s[i], Xj = x[j];

        double sigS = sigma_s(i, n);
        double sigX = sigma_x(j, n);

        // X-direction coefficients on non-uniform mesh
        aX[j - 1] = ((r_d - r_f) * Xj * h_j - std::pow(sigX * Xj, 2.0)) / (h_jm1 * (h_jm1 + h_j));
        bX[j - 1] = 1.0 / dt
        + (std::pow(sigX * Xj, 2.0) - (r_d - r_f) * Xj * (h_j - h_jm1)) / (h_jm1 * h_j)
          + 0.5 * r_d;
          cX[j - 1] = (-(r_d - r_f) * Xj * h_jm1 - std::pow(sigX * Xj, 2.0)) / (h_j * (h_jm1 + h_j));

          // Explicit cross term using V
          double denom_cross = (h_jm1 * hs[i - 1] + h_j * hs[i] + h_jm1 * hs[i - 1] + h_j * hs[i]);
          double cross = 0.0;
          if (denom_cross != 0.0) {
            cross = 0.5 * rho * sigS * sigX * Si * Xj *
              (V(i + 1, j + 1) + V(i - 1, j - 1) - V(i - 1, j + 1) - V(i + 1, j - 1)) / denom_cross;
          }

          fX[j - 1] = V(i, j) / dt + cross;
      }

      // Boundary contributions to RHS in X, then zero those coefficients
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

    // Enforce boundaries on Unew frame
    for (int i = 0; i <= n_s; ++i) { Unew(i, 0) = leftX[i];  Unew(i, n_x) = rightX[i]; }
    for (int j = 0; j <= n_x; ++j) { Unew(0, j) = leftS[j];  Unew(n_s, j) = rightS[j]; }

    // Move to next time level
    U = clone(Unew);
  }

  // Bilinear interpolation at (s_0, x_0)
  return bilinear_interpolation(s, x, U, s_0, x_0);
}

}

//' European Option 2D (local volatility)
//'
//' Prices a European option on an equity in a foreign currency using a 2D finite-difference PDE
//' with time- and state-dependent local volatilities for both S and X.
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
//'
//' @details
//' Domestic-measure drift and discounting are applied as in the constant-vol solver, but with
//' pointwise volatilities. In the S-direction, mu_S = r_f - q - rho * sigma_S * sigma_X (local).
//' In the X-direction, drift is r_d - r_f. Far-field Dirichlet boundaries use exp(-q * tau) on S·X
//' and exp(-r_d * tau) on K. When `sigma_s` and `sigma_x` are constant-valued matrices, prices
//' should match the constant-vol solver up to numerical tolerance.
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
