// [[Rcpp::interfaces(r, cpp)]]

#include <Rcpp.h>
#include <thread>
#include <vector>
#include <random>
#include <cmath>
#include <algorithm>

using namespace Rcpp;

namespace {

bool cholesky4_am(const double R[4][4], double L[4][4]) {
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

// Solve normal equations (X'X) beta = X'y for small systems
// X is n x p, y is n x 1, beta is p x 1
// Uses simple Gauss elimination for the p x p system
void solve_normal_equations(const std::vector<double>& X, // n*p column-major
                            const std::vector<double>& y, // n
                            int n, int p,
                            std::vector<double>& beta) {
  // Build X'X (p x p) and X'y (p x 1)
  std::vector<double> XtX(p * p, 0.0);
  std::vector<double> Xty(p, 0.0);

  for (int j = 0; j < p; ++j) {
    for (int k = 0; k <= j; ++k) {
      double sum = 0.0;
      for (int i = 0; i < n; ++i)
        sum += X[i + j * n] * X[i + k * n];
      XtX[j * p + k] = sum;
      XtX[k * p + j] = sum;
    }
    double sum = 0.0;
    for (int i = 0; i < n; ++i)
      sum += X[i + j * n] * y[i];
    Xty[j] = sum;
  }

  // Add small ridge for numerical stability
  for (int j = 0; j < p; ++j)
    XtX[j * p + j] += 1e-10;

  // Gauss elimination with partial pivoting
  std::vector<double> A(p * (p + 1));
  for (int i = 0; i < p; ++i) {
    for (int j = 0; j < p; ++j)
      A[i * (p + 1) + j] = XtX[i * p + j];
    A[i * (p + 1) + p] = Xty[i];
  }

  for (int col = 0; col < p; ++col) {
    // Pivot
    int max_row = col;
    double max_val = std::abs(A[col * (p + 1) + col]);
    for (int row = col + 1; row < p; ++row) {
      double val = std::abs(A[row * (p + 1) + col]);
      if (val > max_val) { max_val = val; max_row = row; }
    }
    if (max_row != col) {
      for (int j = 0; j <= p; ++j)
        std::swap(A[col * (p + 1) + j], A[max_row * (p + 1) + j]);
    }

    double pivot = A[col * (p + 1) + col];
    if (std::abs(pivot) < 1e-15) {
      beta.assign(p, 0.0);
      return;
    }

    for (int row = col + 1; row < p; ++row) {
      double factor = A[row * (p + 1) + col] / pivot;
      for (int j = col; j <= p; ++j)
        A[row * (p + 1) + j] -= factor * A[col * (p + 1) + j];
    }
  }

  // Back substitution
  beta.resize(p);
  for (int i = p - 1; i >= 0; --i) {
    double sum = A[i * (p + 1) + p];
    for (int j = i + 1; j < p; ++j)
      sum -= A[i * (p + 1) + j] * beta[j];
    beta[i] = sum / A[i * (p + 1) + i];
  }
}

} // namespace

//' American Option 4D Monte Carlo (Double-Heston Quanto, Longstaff-Schwartz)
//'
//' Prices an American option on a foreign equity quoted in domestic currency
//' under double Heston stochastic volatility using Longstaff-Schwartz
//' least-squares Monte Carlo with \code{std::thread} parallelism for
//' the forward simulation pass.
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
//' @param n_paths Number of Monte Carlo paths (e.g. 1e5).
//' @param n_steps Number of exercise dates / time steps (e.g. 50).
//' @param seed RNG seed for reproducibility.
//' @param n_basis Number of polynomial basis functions in \code{P/K} for LSM
//'   regression (default 4). The full basis is
//'   \code{1, P/K, ..., (P/K)^(n_basis-1), v_S/theta_S, v_X/theta_X}.
//' @param n_threads Number of threads for forward simulation (0 = auto-detect;
//'   default 0). Set to 1 for fully reproducible results across machines.
//'
//' @return A list with components \code{price} (option price) and
//'   \code{std_error} (Monte Carlo standard error estimate). For tidy-data
//'   workflows, convert this to a one-row data frame/tibble with explicit
//'   columns such as model, exercise, price, and std_error.
//'
//' @details
//' The forward simulation uses Euler-Maruyama with full truncation, identical
//' to \code{\link{mc_european_heston_4d}}. The product \eqn{P_t = S_t X_t} is
//' stored at each exercise date.
//'
//' The backward Longstaff-Schwartz pass regresses continuation values on
//' polynomial basis functions of P and the current variance levels v_S, v_X
//' at each exercise date. The basis is \code{1, P/K, ..., (P/K)^(n_basis-1),
//' v_S/theta_S, v_X/theta_X} (n_basis+2 regressors total). The regression
//' determines the optimal exercise boundary. Only in-the-money paths are
//' used in the regression.
//'
//' @examples
//' mc_american_heston_4d(100, 20, 2000, 1, 0.05, 0.02, 0.01,
//'   2, 0.04, 0.5, 0.04, 1.5, 0.02, 0.3, 0.02,
//'   0.3, -0.7, -0.5, 0, 0, 0, "put", 50000, 50, 42, 4)
//'
//' @export
// [[Rcpp::export]]
Rcpp::List mc_american_heston_4d(
    double s_0, double x_0, double k, double tau,
    double r_d, double r_f, double q,
    double kappa_s, double theta_s, double xi_s, double v_s0,
    double kappa_x, double theta_x, double xi_x, double v_x0,
    double rho_sx, double rho_sv, double rho_xv,
    double rho_svx, double rho_xvs, double rho_vsvx,
    String type,
    int n_paths, int n_steps, int seed, int n_basis = 4, int n_threads = 0) {

  // Input validation
  if (s_0 <= 0 || x_0 <= 0 || k <= 0) stop("s_0, x_0, and k must be > 0");
  if (tau <= 0)    stop("tau must be > 0");
  if (kappa_s <= 0 || theta_s <= 0 || xi_s <= 0 || v_s0 <= 0)
    stop("kappa_s, theta_s, xi_s, and v_s0 must be > 0");
  if (kappa_x <= 0 || theta_x <= 0 || xi_x <= 0 || v_x0 <= 0)
    stop("kappa_x, theta_x, xi_x, and v_x0 must be > 0");
  if (n_paths < 1) stop("n_paths must be >= 1");
  if (n_steps < 1) stop("n_steps must be >= 1");
  if (n_basis < 2 || n_basis > 8) stop("n_basis must be between 2 and 8");

  // Guard against memory exhaustion: 3 path arrays (P, vs, vx) each n_paths * (n_steps+1)
  const double alloc_bytes = 3.0 * static_cast<double>(n_paths) * (n_steps + 1) * sizeof(double);
  if (alloc_bytes > 2e9)
    stop("Requested allocation (3 * n_paths * n_steps) exceeds 2 GB safety limit. "
         "Reduce n_paths or n_steps.");

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

  if (2.0 * kappa_s * theta_s <= xi_s * xi_s)
    Rcpp::warning("Feller condition violated for equity variance process.");
  if (2.0 * kappa_x * theta_x <= xi_x * xi_x)
    Rcpp::warning("Feller condition violated for FX variance process.");

  // Correlation matrix and Cholesky
  double R[4][4] = {
    {1.0,    rho_sx,  rho_sv,  rho_svx},
    {rho_sx, 1.0,     rho_xvs, rho_xv},
    {rho_sv, rho_xvs, 1.0,     rho_vsvx},
    {rho_svx,rho_xv,  rho_vsvx,1.0}
  };

  double L[4][4];
  if (!cholesky4_am(R, L))
    stop("Correlation matrix is not positive semi-definite.");

  const double dt = tau / n_steps;
  const double sqrt_dt = std::sqrt(dt);
  const double inv_k = 1.0 / k;
  const double inv_theta_s = 1.0 / theta_s;
  const double inv_theta_x = 1.0 / theta_x;

  std::vector<double> disc_steps(n_steps + 1);
  for (int h = 0; h <= n_steps; ++h)
    disc_steps[h] = std::exp(-r_d * h * dt);


  // ======================================================================= //
  // Forward pass: simulate paths and store product P = S*X at each step
  // Storage: n_paths x (n_steps+1), column-major
  // ======================================================================= //
  std::vector<double> P_paths(static_cast<size_t>(n_paths) * (n_steps + 1));
  std::vector<double> vs_paths(static_cast<size_t>(n_paths) * (n_steps + 1));
  std::vector<double> vx_paths(static_cast<size_t>(n_paths) * (n_steps + 1));

  // Set initial values
  const double P_0 = s_0 * x_0;
  for (int p = 0; p < n_paths; ++p) {
    P_paths[p]  = P_0;
    vs_paths[p] = v_s0;
    vx_paths[p] = v_x0;
  }

  // Parallel forward simulation
  unsigned int n_thr = (n_threads > 0)
    ? static_cast<unsigned int>(n_threads)
    : std::thread::hardware_concurrency();
  if (n_thr == 0) n_thr = 1;
  if (n_thr > static_cast<unsigned int>(n_paths)) n_thr = n_paths;

  auto simulate_forward = [&](unsigned int tid, int path_begin, int path_end) {
    std::mt19937_64 rng(static_cast<uint64_t>(seed) + tid * 1000003ULL);
    std::normal_distribution<double> randn(0.0, 1.0);

    for (int p = path_begin; p < path_end; ++p) {
      double log_s = std::log(s_0);
      double log_x = std::log(x_0);
      double vs = v_s0;
      double vx = v_x0;

      for (int t = 0; t < n_steps; ++t) {
        double z1 = randn(rng), z2 = randn(rng);
        double z3 = randn(rng), z4 = randn(rng);

        double w0 = L[0][0]*z1;
        double w1 = L[1][0]*z1 + L[1][1]*z2;
        double w2 = L[2][0]*z1 + L[2][1]*z2 + L[2][2]*z3;
        double w3 = L[3][0]*z1 + L[3][1]*z2 + L[3][2]*z3 + L[3][3]*z4;

        double vs_pos = std::max(vs, 0.0);
        double vx_pos = std::max(vx, 0.0);
        double sqrt_vs = std::sqrt(vs_pos);
        double sqrt_vx = std::sqrt(vx_pos);

        double drift_s = r_f - q - rho_sx * sqrt_vs * sqrt_vx;

        log_s += (drift_s - 0.5 * vs_pos) * dt + sqrt_vs * sqrt_dt * w0;
        log_x += (r_d - r_f - 0.5 * vx_pos) * dt + sqrt_vx * sqrt_dt * w1;
        vs    += kappa_s * (theta_s - vs_pos) * dt + xi_s * sqrt_vs * sqrt_dt * w2;
        vx    += kappa_x * (theta_x - vx_pos) * dt + xi_x * sqrt_vx * sqrt_dt * w3;

        // Store product and variances at step t+1
        double S_t = std::exp(log_s);
        double X_t = std::exp(log_x);
        size_t col = static_cast<size_t>(t + 1) * n_paths;
        P_paths [p + col] = S_t * X_t;
        vs_paths[p + col] = vs;
        vx_paths[p + col] = vx;
      }
    }
  };

  {
    std::vector<std::thread> threads;
    threads.reserve(n_thr);
    int chunk = n_paths / n_thr;
    int rem = n_paths % n_thr;
    int start = 0;
    for (unsigned int t = 0; t < n_thr; ++t) {
      int end = start + chunk + (static_cast<int>(t) < rem ? 1 : 0);
      threads.emplace_back(simulate_forward, t, start, end);
      start = end;
    }
    for (auto& th : threads) th.join();
  }

  // ======================================================================= //
  // Backward pass: Longstaff-Schwartz
  // ======================================================================= //

  // Cash flow matrix: for each path, the (discounted to time 0) payoff
  // from optimal exercise. Initialize with terminal payoff.
  std::vector<double> cashflow(n_paths);
  std::vector<int> exercise_time(n_paths, n_steps);

  for (int p = 0; p < n_paths; ++p) {
    double Pt = P_paths[p + static_cast<size_t>(n_steps) * n_paths];
    cashflow[p] = is_call ? std::max(Pt - k, 0.0) : std::max(k - Pt, 0.0);
  }

  const int n_reg = n_basis + 2;
  std::vector<int> itm_indices;
  itm_indices.reserve(n_paths / 2);
  std::vector<double> Xmat;
  std::vector<double> Y;
  std::vector<double> beta;

  // Backward from n_steps-1 to 1 (step 0 is t=0, no exercise there)
  for (int step = n_steps - 1; step >= 1; --step) {
    const size_t step_offset = static_cast<size_t>(step) * n_paths;

    // Find ITM paths at this step
    itm_indices.clear();

    for (int p = 0; p < n_paths; ++p) {
      double Pt = P_paths[p + step_offset];
      double exercise_val = is_call ? std::max(Pt - k, 0.0) : std::max(k - Pt, 0.0);
      if (exercise_val > 0.0) itm_indices.push_back(p);
    }

    if (itm_indices.empty()) continue;

    int n_itm = static_cast<int>(itm_indices.size());

    // Build design matrix X (n_itm x n_reg) and response Y
    Xmat.resize(static_cast<size_t>(n_itm) * n_reg, 0.0); // column-major
    Y.resize(n_itm);

    for (int i = 0; i < n_itm; ++i) {
      int p = itm_indices[i];
      double Pt  = P_paths [p + step_offset];
      double vst = std::max(vs_paths[p + step_offset], 0.0);
      double vxt = std::max(vx_paths[p + step_offset], 0.0);
      double p_over_k = Pt * inv_k;

      // Polynomial basis: 1, P/K, (P/K)^2, ...
      double basis = 1.0;
      for (int b = 0; b < n_basis; ++b) {
        Xmat[i + b * n_itm] = basis;
        basis *= p_over_k;
      }
      // Variance regressors: normalized by long-run mean for scale
      Xmat[i + n_basis       * n_itm] = vst * inv_theta_s;
      Xmat[i + (n_basis + 1) * n_itm] = vxt * inv_theta_x;

      // Continuation value: discounted future cashflow from this path
      int steps_fwd = exercise_time[p] - step;
      double disc_fwd = disc_steps[steps_fwd];
      Y[i] = cashflow[p] * disc_fwd;
    }

    // Solve regression
    beta.clear();
    if (n_itm >= n_reg) {
      solve_normal_equations(Xmat, Y, n_itm, n_reg, beta);
    } else {
      continue; // not enough ITM paths for regression
    }

    // Compare exercise vs continuation for ITM paths
    for (int i = 0; i < n_itm; ++i) {
      int p = itm_indices[i];
      double Pt  = P_paths [p + step_offset];
      double vst = std::max(vs_paths[p + step_offset], 0.0);
      double vxt = std::max(vx_paths[p + step_offset], 0.0);
      double p_over_k = Pt * inv_k;
      double exercise_val = is_call ? std::max(Pt - k, 0.0) : std::max(k - Pt, 0.0);

      // Estimated continuation
      double continuation = 0.0;
      double basis = 1.0;
      for (int b = 0; b < n_basis; ++b) {
        continuation += beta[b] * basis;
        basis *= p_over_k;
      }
      continuation += beta[n_basis]       * (vst * inv_theta_s);
      continuation += beta[n_basis + 1]   * (vxt * inv_theta_x);

      if (exercise_val > continuation) {
        cashflow[p] = exercise_val;
        exercise_time[p] = step;
      }
    }
  }

  // ======================================================================= //
  // Compute price: average of discounted cashflows
  // ======================================================================= //
  double sum_disc_cf = 0.0;
  double sum_disc_cf_sq = 0.0;

  for (int p = 0; p < n_paths; ++p) {
    double disc_cf = cashflow[p] * disc_steps[exercise_time[p]];
    sum_disc_cf += disc_cf;
    sum_disc_cf_sq += disc_cf * disc_cf;
  }

  double price = sum_disc_cf / n_paths;
  double variance = sum_disc_cf_sq / n_paths - price * price;
  double std_error = std::sqrt(std::max(variance, 0.0) / n_paths);

  return Rcpp::List::create(
    Rcpp::Named("price") = price,
    Rcpp::Named("std_error") = std_error
  );
}
