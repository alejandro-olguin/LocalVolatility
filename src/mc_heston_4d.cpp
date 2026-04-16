// [[Rcpp::interfaces(r, cpp)]]

#include <Rcpp.h>
#include <thread>
#include <vector>
#include <random>
#include <cmath>
#include <algorithm>

using namespace Rcpp;

namespace {

// 4x4 Cholesky decomposition (lower triangular L such that R = L L^T)
// Returns false if matrix is not positive semi-definite.
bool cholesky4(const double R[4][4], double L[4][4]) {
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j)
      L[i][j] = 0.0;

  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j <= i; ++j) {
      double sum = R[i][j];
      for (int k = 0; k < j; ++k) sum -= L[i][k] * L[j][k];
      if (i == j) {
        if (sum <= 0.0) return false;
        L[i][j] = std::sqrt(sum);
      } else {
        L[i][j] = sum / L[j][j];
      }
    }
  }
  return true;
}

} // namespace

//' European Option 4D Monte Carlo (Double-Heston Quanto)
//'
//' Prices a European option on a foreign equity quoted in domestic currency
//' under double Heston stochastic volatility (one Heston process for the
//' equity, one for the FX rate) using Monte Carlo simulation with
//' \code{std::thread} parallelism and antithetic variates.
//'
//' @param s_0 Equity spot price (> 0).
//' @param x_0 FX spot price, domestic per foreign (> 0).
//' @param k Strike price in domestic currency (> 0).
//' @param tau Time to expiry in years (> 0).
//' @param r_d Domestic risk-free rate.
//' @param r_f Foreign risk-free rate.
//' @param q Equity dividend yield.
//' @param kappa_s Mean reversion speed for equity variance (> 0).
//' @param theta_s Long-run equity variance (> 0).
//' @param xi_s Vol-of-vol for equity (> 0).
//' @param v_s0 Initial equity variance (> 0).
//' @param kappa_x Mean reversion speed for FX variance (> 0).
//' @param theta_x Long-run FX variance (> 0).
//' @param xi_x Vol-of-vol for FX (> 0).
//' @param v_x0 Initial FX variance (> 0).
//' @param rho_sx Correlation between S and X.
//' @param rho_sv Correlation between S and v_S.
//' @param rho_xv Correlation between X and v_X.
//' @param rho_svx Correlation between S and v_X (typically ~0).
//' @param rho_xvs Correlation between X and v_S (typically ~0).
//' @param rho_vsvx Correlation between v_S and v_X (typically ~0).
//' @param type Either \code{"call"} or \code{"put"}.
//' @param n_paths Number of Monte Carlo paths (e.g. 1e6).
//' @param n_steps Number of time steps per path (e.g. 100).
//' @param seed RNG seed for reproducibility.
//'
//' @return A list with components \code{price} (option price) and
//'   \code{std_error} (Monte Carlo standard error).
//'
//' @details
//' Simulates the 4D system (S, X, v_S, v_X) under the domestic measure using
//' Euler-Maruyama with full truncation for variance processes. Correlated
//' Brownian motions are generated via Cholesky decomposition of the 6-parameter
//' correlation matrix. Antithetic variates halve the variance at no extra
//' simulation cost.
//'
//' The domestic-measure dynamics are:
//' \deqn{dS = (r_f - q - \rho_{SX}\sqrt{v_S}\sqrt{v_X}) S \, dt + \sqrt{v_S} S \, dW_S}
//' \deqn{dX = (r_d - r_f) X \, dt + \sqrt{v_X} X \, dW_X}
//' \deqn{dv_S = \kappa_S(\theta_S - v_S) \, dt + \xi_S \sqrt{v_S} \, dW_{v_S}}
//' \deqn{dv_X = \kappa_X(\theta_X - v_X) \, dt + \xi_X \sqrt{v_X} \, dW_{v_X}}
//'
//' @examples
//' mc_european_heston_4d(100, 20, 2000, 1, 0.05, 0.02, 0.01,
//'   2, 0.04, 0.5, 0.04, 1.5, 0.02, 0.3, 0.02,
//'   0.3, -0.7, -0.5, 0, 0, 0, "call", 100000, 100, 42)
//'
//' @export
// [[Rcpp::export]]
Rcpp::List mc_european_heston_4d(
    double s_0, double x_0, double k, double tau,
    double r_d, double r_f, double q,
    double kappa_s, double theta_s, double xi_s, double v_s0,
    double kappa_x, double theta_x, double xi_x, double v_x0,
    double rho_sx, double rho_sv, double rho_xv,
    double rho_svx, double rho_xvs, double rho_vsvx,
    String type,
    int n_paths, int n_steps, int seed) {

  // Input validation
  if (s_0 <= 0 || x_0 <= 0 || k <= 0) stop("s_0, x_0, and k must be > 0");
  if (tau <= 0)    stop("tau must be > 0");
  if (kappa_s <= 0 || theta_s <= 0 || xi_s <= 0 || v_s0 <= 0)
    stop("kappa_s, theta_s, xi_s, and v_s0 must be > 0");
  if (kappa_x <= 0 || theta_x <= 0 || xi_x <= 0 || v_x0 <= 0)
    stop("kappa_x, theta_x, xi_x, and v_x0 must be > 0");
  if (n_paths < 1) stop("n_paths must be >= 1");
  if (n_steps < 1) stop("n_steps must be >= 1");

  // Validate correlations
  auto check_rho = [](double r, const char* name) {
    if (r < -1.0 || r > 1.0) {
      std::string msg = std::string(name) + " must be within [-1, 1]";
      Rcpp::stop(msg.c_str());
    }
  };
  check_rho(rho_sx, "rho_sx"); check_rho(rho_sv, "rho_sv");
  check_rho(rho_xv, "rho_xv"); check_rho(rho_svx, "rho_svx");
  check_rho(rho_xvs, "rho_xvs"); check_rho(rho_vsvx, "rho_vsvx");

  const bool is_call = (type == "call");
  if (!is_call && type != "put") stop("type must be \"call\" or \"put\"");

  // Feller condition warnings
  if (2.0 * kappa_s * theta_s <= xi_s * xi_s)
    Rcpp::warning("Feller condition violated for equity variance process.");
  if (2.0 * kappa_x * theta_x <= xi_x * xi_x)
    Rcpp::warning("Feller condition violated for FX variance process.");

  // Build 4x4 correlation matrix: order (S, X, v_S, v_X)
  double R[4][4] = {
    {1.0,    rho_sx,  rho_sv,  rho_svx},
    {rho_sx, 1.0,     rho_xvs, rho_xv},
    {rho_sv, rho_xvs, 1.0,     rho_vsvx},
    {rho_svx,rho_xv,  rho_vsvx,1.0}
  };

  double L[4][4];
  if (!cholesky4(R, L))
    stop("Correlation matrix is not positive semi-definite.");

  const double dt = tau / n_steps;
  const double sqrt_dt = std::sqrt(dt);
  const double disc = std::exp(-r_d * tau);

  // Thread setup
  unsigned int n_threads = std::thread::hardware_concurrency();
  if (n_threads == 0) n_threads = 1;
  // Each antithetic pair = 1 "path unit"; n_paths is actual number of path pairs


  // Per-thread accumulators
  std::vector<double> thread_sum(n_threads, 0.0);
  std::vector<double> thread_sum_sq(n_threads, 0.0);

  // Worker lambda
  auto simulate = [&](unsigned int tid, int path_begin, int path_end) {
    std::mt19937_64 rng(static_cast<uint64_t>(seed) + tid * 1000003ULL);
    std::normal_distribution<double> randn(0.0, 1.0);

    double local_sum = 0.0;
    double local_sum_sq = 0.0;

    for (int p = path_begin; p < path_end; ++p) {
      // Simulate two paths (antithetic pair)
      double log_s[2]  = {std::log(s_0), std::log(s_0)};
      double log_x[2]  = {std::log(x_0), std::log(x_0)};
      double vs[2]     = {v_s0, v_s0};
      double vx[2]     = {v_x0, v_x0};

      for (int t = 0; t < n_steps; ++t) {
        double z1 = randn(rng), z2 = randn(rng);
        double z3 = randn(rng), z4 = randn(rng);

        // Correlated Brownians via Cholesky
        double w[4];
        w[0] = L[0][0]*z1;
        w[1] = L[1][0]*z1 + L[1][1]*z2;
        w[2] = L[2][0]*z1 + L[2][1]*z2 + L[2][2]*z3;
        w[3] = L[3][0]*z1 + L[3][1]*z2 + L[3][2]*z3 + L[3][3]*z4;

        // Antithetic Brownians
        double wa[4] = {-w[0], -w[1], -w[2], -w[3]};

        for (int a = 0; a < 2; ++a) {
          double* ww = (a == 0) ? w : wa;
          double vs_pos = std::max(vs[a], 0.0);
          double vx_pos = std::max(vx[a], 0.0);
          double sqrt_vs = std::sqrt(vs_pos);
          double sqrt_vx = std::sqrt(vx_pos);

          // Quanto drift adjustment
          double drift_s = r_f - q - rho_sx * sqrt_vs * sqrt_vx;

          log_s[a] += (drift_s - 0.5 * vs_pos) * dt + sqrt_vs * sqrt_dt * ww[0];
          log_x[a] += (r_d - r_f - 0.5 * vx_pos) * dt + sqrt_vx * sqrt_dt * ww[1];
          vs[a]    += kappa_s * (theta_s - vs_pos) * dt + xi_s * sqrt_vs * sqrt_dt * ww[2];
          vx[a]    += kappa_x * (theta_x - vx_pos) * dt + xi_x * sqrt_vx * sqrt_dt * ww[3];
        }
      }

      // Payoff for both paths
      for (int a = 0; a < 2; ++a) {
        double S_T = std::exp(log_s[a]);
        double X_T = std::exp(log_x[a]);
        double payoff = is_call ? std::max(S_T * X_T - k, 0.0)
                                : std::max(k - S_T * X_T, 0.0);
        double disc_payoff = disc * payoff;
        local_sum += disc_payoff;
        local_sum_sq += disc_payoff * disc_payoff;
      }
    }

    thread_sum[tid] = local_sum;
    thread_sum_sq[tid] = local_sum_sq;
  };

  // Launch threads
  if (n_threads > static_cast<unsigned int>(n_paths)) n_threads = n_paths;

  std::vector<std::thread> threads;
  threads.reserve(n_threads);

  int chunk = n_paths / n_threads;
  int remainder = n_paths % n_threads;
  int start = 0;

  for (unsigned int t = 0; t < n_threads; ++t) {
    int end = start + chunk + (static_cast<int>(t) < remainder ? 1 : 0);
    threads.emplace_back(simulate, t, start, end);
    start = end;
  }

  for (auto& th : threads) th.join();

  // Aggregate results
  double total_sum = 0.0, total_sum_sq = 0.0;
  for (unsigned int t = 0; t < n_threads; ++t) {
    total_sum += thread_sum[t];
    total_sum_sq += thread_sum_sq[t];
  }

  // n_paths pairs × 2 paths each = total 2*n_paths individual payoffs
  const double total_paths = 2.0 * n_paths;
  const double price = total_sum / total_paths;
  const double variance = total_sum_sq / total_paths - price * price;
  const double std_error = std::sqrt(std::max(variance, 0.0) / total_paths);

  return Rcpp::List::create(
    Rcpp::Named("price") = price,
    Rcpp::Named("std_error") = std_error
  );
}
