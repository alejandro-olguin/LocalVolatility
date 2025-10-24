// [[Rcpp::interfaces(r, cpp)]]

#include <Rcpp.h>
#include "utils.h"

using namespace Rcpp;

namespace LocalVolatility {

double american_option_2d(double s_0,
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
                          double alpha,
                          double lambda,
                          double tolerance) {

  // --- guards ---
  if (n_s < 3 || n_x < 3) stop("n_s and n_x must be >= 3");
  if (n_t < 1)            stop("n_t must be >= 1");
  if (s_max <= s_min || x_max <= x_min) stop("grid bounds must be increasing");
  if (tau <= 0.0)         stop("tau must be > 0");

  // Initialize grids (Tavella–Randall; returns n+1 nodes)
  NumericVector s = tavella_randall(s_0, alpha, s_min, s_max, n_s);
  NumericVector x = tavella_randall(x_0, alpha, x_min, x_max, n_x);

  // Nonuniform steps
  NumericVector hs(n_s), hx(n_x);
  for (int i = 0; i < n_s; ++i) hs[i] = s[i + 1] - s[i];
  for (int j = 0; j < n_x; ++j) hx[j] = x[j + 1] - x[j];

  const double dt = tau / n_t;
  const double half_lambda = 0.5 * lambda; // split penalty over the two substeps

  // Payoff on product S*X
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

  // Start value at maturity
  NumericMatrix u = clone(payoff);

  // --- far-field Dirichlet boundaries (more stable than linear extrapolation) ---
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

  // Tridiagonal buffers (proper sub/diag/super)
  NumericVector aS(n_s - 1), bS(n_s - 1), cS(n_s - 1), fS(n_s - 1);
  NumericVector aX(n_x - 1), bX(n_x - 1), cX(n_x - 1), fX(n_x - 1);

  // Main time loop (backward)
  for (int n = 0; n < n_t; ++n) {
    const double tau_rem = (n_t - (n + 1)) * dt;

    // Precompute boundaries for this step
    NumericVector leftS, rightS, leftX, rightX;
    set_S_boundaries(tau_rem, leftS, rightS);
    set_X_boundaries(tau_rem, leftX, rightX);

    // Penalty fixed-point iteration at this time level
    NumericMatrix u_i = clone(u);
    double rel_err = 1e6;

    while (rel_err > tolerance) {
      // Penalty indicator from current iterate (λ/2 in this substep)
      NumericMatrix p_u(n_s + 1, n_x + 1);
      for (int i = 0; i <= n_s; ++i)
        for (int j = 0; j <= n_x; ++j)
          p_u(i, j) = (u_i(i, j) < payoff(i, j)) ? half_lambda : 0.0;

      // ===== Pass 1: implicit in S (per fixed j), explicit in X & cross =====
      NumericMatrix v(n_s + 1, n_x + 1);
      for (int j = 1; j < n_x; ++j) {
        const double L = leftS[j], R = rightS[j];

        for (int i = 1; i < n_s; ++i) {
          const double h_im1 = hs[i - 1], h_i = hs[i];
          const double Si = s[i], Xj = x[j];

          // constant vol in S
          const double a2S = std::pow(sigma_s * Si, 2.0);

          // domestic-measure drift for S
          const double muS = r_f - q - rho * sigma_s * sigma_x;

          aS[i - 1] = ( (muS * Si * h_i) - a2S ) / (h_im1 * (h_im1 + h_i));
          bS[i - 1] = 1.0 / dt
          + ( a2S - muS * Si * (h_i - h_im1) ) / (h_im1 * h_i)
            + 0.5 * r_d
            + 0.5 * p_u(i, j);
            cS[i - 1] = ( -(muS * Si * h_im1) - a2S ) / (h_i * (h_im1 + h_i));

            // explicit mixed term (centered)
            const double denom = (h_im1 * hx[j - 1] + h_i * hx[j]
                                    + h_i * hx[j - 1] + h_im1 * hx[j]);
            double cross = 0.0;
            if (denom != 0.0) {
              cross = 0.5 * rho * sigma_s * sigma_x * Si * Xj *
                (u(i + 1, j + 1) + u(i - 1, j - 1) - u(i - 1, j + 1) - u(i + 1, j - 1)) / denom;
            }

            fS[i - 1] = u(i, j) / dt + cross + 0.5 * p_u(i, j) * payoff(i, j);
        }

        // boundary contributions to RHS then zero touching coefficients
        fS[0]         -= aS[0]       * L;
        fS[n_s - 2]   -= cS[n_s - 2] * R;
        aS[0]          = 0.0;
        cS[n_s - 2]    = 0.0;

        // solve in S
        NumericVector col = thomas_algorithm(aS, bS, cS, fS);
        v(0, j)    = L;
        v(n_s, j)  = R;
        for (int i = 1; i < n_s; ++i) v(i, j) = col[i - 1];
      }

      // enforce frame boundaries on v
      for (int i = 0; i <= n_s; ++i) { v(i, 0) = leftX[i];  v(i, n_x) = rightX[i]; }
      for (int j = 0; j <= n_x; ++j) { v(0, j) = leftS[j];  v(n_s, j) = rightS[j]; }

      // ===== Pass 2: implicit in X (per fixed i), explicit in S & cross =====
      NumericMatrix u_new(n_s + 1, n_x + 1);

      // penalty from v for this pass (again λ/2)
      NumericMatrix p_v(n_s + 1, n_x + 1);
      for (int i = 0; i <= n_s; ++i)
        for (int j = 0; j <= n_x; ++j)
          p_v(i, j) = (v(i, j) < payoff(i, j)) ? half_lambda : 0.0;

      for (int i = 1; i < n_s; ++i) {
        const double L = leftX[i], R = rightX[i];

        for (int j = 1; j < n_x; ++j) {
          const double h_jm1 = hx[j - 1], h_j = hx[j];
          const double Si = s[i], Xj = x[j];

          const double a2X = std::pow(sigma_x * Xj, 2.0);

          aX[j - 1] = ((r_d - r_f) * Xj * h_j - a2X) / (h_jm1 * (h_jm1 + h_j));
          bX[j - 1] = 1.0 / dt
          + (a2X - (r_d - r_f) * Xj * (h_j - h_jm1)) / (h_jm1 * h_j)
            + 0.5 * r_d
            + 0.5 * p_v(i, j);
            cX[j - 1] = (-(r_d - r_f) * Xj * h_jm1 - a2X) / (h_j * (h_jm1 + h_j));

            const double denom = (h_jm1 * hs[i - 1] + h_j * hs[i]
                                    + h_jm1 * hs[i - 1] + h_j * hs[i]);
            double cross = 0.0;
            if (denom != 0.0) {
              cross = 0.5 * rho * sigma_s * sigma_x * Si * Xj *
                (v(i + 1, j + 1) + v(i - 1, j - 1) - v(i - 1, j + 1) - v(i + 1, j - 1)) / denom;
            }

            fX[j - 1] = v(i, j) / dt + cross + 0.5 * p_v(i, j) * payoff(i, j);
        }

        // boundary contributions (X) then zero touching coefficients
        fX[0]         -= aX[0]       * L;
        fX[n_x - 2]   -= cX[n_x - 2] * R;
        aX[0]          = 0.0;
        cX[n_x - 2]    = 0.0;

        // solve in X
        NumericVector row = thomas_algorithm(aX, bX, cX, fX);
        u_new(i, 0)   = L;
        u_new(i, n_x) = R;
        for (int j = 1; j < n_x; ++j) u_new(i, j) = row[j - 1];
      }

      // enforce boundaries on u_new
      for (int i = 0; i <= n_s; ++i) { u_new(i, 0) = leftX[i];  u_new(i, n_x) = rightX[i]; }
      for (int j = 0; j <= n_x; ++j) { u_new(0, j) = leftS[j];  u_new(n_s, j) = rightS[j]; }

      // project to obstacle + compute relative error
      double num = 0.0, den = 0.0;
      for (int i = 0; i <= n_s; ++i)
        for (int j = 0; j <= n_x; ++j) {
          const double vproj = std::max(u_new(i, j), payoff(i, j));
          const double diff  = vproj - u_i(i, j);
          num += diff * diff;
          den += u_i(i, j) * u_i(i, j);
          u_i(i, j) = vproj;
        }
        rel_err = std::sqrt(num / (den + 1e-16));
    }

    u = clone(u_i);
  }

  // Price by bilinear interpolation
  return bilinear_interpolation(s, x, u, s_0, x_0);
}
}

//' American Option 2D (constant volatility, penalty method)
//'
//' Prices an American option on an equity quoted in a foreign currency using a 2D finite-difference PDE
//' with Yanenko operator splitting and a penalty-projection scheme for the early exercise constraint.
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
//' @param lambda Penalty parameter (> 1).
//' @param tolerance Relative error tolerance for the penalty fixed-point iterations.
//'
//' @details
//' Domestic-measure drift and discounting are applied (mu_S = r_f - q - rho * sigma_s * sigma_x in S,
//' r_d - r_f in X, discounting at r_d). Far-field Dirichlet boundaries use exp(-q * tau) on S·X and
//' exp(-r_d * tau) on K. The penalty parameter is split across the two sub-steps.
//'
// [[Rcpp::export]]
 double american_option_2d(double s_0,
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
                           double alpha,
                           double lambda,
                           double tolerance) {
   return LocalVolatility::american_option_2d(s_0, x_0, k, tau, r_d, r_f, q,
                                              sigma_s, sigma_x, rho, type,
                                              s_min, s_max, x_min, x_max,
                                              n_s, n_x, n_t, alpha, lambda, tolerance);
 }
