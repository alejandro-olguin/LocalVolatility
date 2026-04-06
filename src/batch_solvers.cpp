// [[Rcpp::interfaces(r, cpp)]]

#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <cmath>
#include <string>
#include <thread>
#include <functional>

#include "thomas_factored.h"

using namespace Rcpp;

// ============================================================================
// Batch European solver: prices multiple European options sharing one sigma
// surface in a single backward PDE sweep with pre-factored Thomas.
// ============================================================================

namespace {

// Raw-pointer linear interpolation (no R allocations)
inline double interp_raw(double s0, double s_min, double ds,
                         const double* prices, int n_nodes) {
  if (s0 <= s_min) return prices[0];
  double s_max = s_min + ds * (n_nodes - 1);
  if (s0 >= s_max) return prices[n_nodes - 1];
  double pos = (s0 - s_min) / ds;
  int i = static_cast<int>(pos);
  double w = pos - i;
  return (1.0 - w) * prices[i] + w * prices[i + 1];
}

} // anonymous namespace


//' Batch-price European options under local volatility (1D)
//'
//' Prices multiple European options that share the same local volatility surface
//' and PDE grid in a single backward Crank-Nicolson sweep. Much faster than
//' calling \code{\link{european_option_lv}} in a loop because the tridiagonal
//' LHS is factored once per time step and reused across all options.
//'
//' @param spots Numeric vector of spot prices (one per option).
//' @param strikes Numeric vector of strike prices.
//' @param taus Numeric vector of times to expiry (years).
//' @param r_ds Numeric vector of domestic risk-free rates.
//' @param r_fs Numeric vector of foreign risk-free rates or cost-of-carry rates.
//' @param sigma Local volatility matrix of size \code{(n_s + 1) x (n_t + 1)}.
//' @param types Character vector, each element \code{"call"} or \code{"put"}.
//' @param s_min,s_max Domain bounds for the asset grid.
//' @param n_s Number of spatial intervals (\code{n_s + 1} nodes).
//' @param n_t Number of time steps.
//'
//' @return Numeric vector of option prices (same length as \code{strikes}).
//'
//' @details
//' All options must lie on the same spatial domain \code{[s_min, s_max]} and
//' use the same sigma surface. The solver steps backward from \code{max(taus)}
//' and extracts each option's price at the time step closest to its maturity.
//' The tridiagonal Crank-Nicolson system is factored once per time step (the
//' LHS depends only on sigma and the grid, not on the strike), and each
//' option's RHS is solved via forward/back substitution — O(n_s) per option
//' per time step instead of O(2 * n_s) for a full Thomas solve.
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector batch_price_european_lv(
    const Rcpp::NumericVector& spots,
    const Rcpp::NumericVector& strikes,
    const Rcpp::NumericVector& taus,
    const Rcpp::NumericVector& r_ds,
    const Rcpp::NumericVector& r_fs,
    const Rcpp::NumericMatrix& sigma,
    const Rcpp::StringVector& types,
    double s_min, double s_max,
    int n_s, int n_t) {

  const int n_opt = strikes.size();

  // --- Validation ---
  if (spots.size() != n_opt || taus.size() != n_opt || r_ds.size() != n_opt ||
      r_fs.size() != n_opt || types.size() != n_opt)
    Rcpp::stop("All input vectors must have the same length.");
  if (n_s < 2 || n_t < 1) Rcpp::stop("n_s >= 2 and n_t >= 1 required.");
  if (s_max <= s_min) Rcpp::stop("s_max must be > s_min.");
  if (sigma.nrow() != n_s + 1 || sigma.ncol() != n_t + 1)
    Rcpp::stop("sigma must have dimensions (n_s+1) x (n_t+1).");

  const double ds = (s_max - s_min) / n_s;
  const int N = n_s - 1; // interior nodes

  // Find max tau to determine the time grid
  double tau_max = 0.0;
  for (int j = 0; j < n_opt; ++j)
    if (taus[j] > tau_max) tau_max = taus[j];
  if (tau_max <= 0.0) Rcpp::stop("All taus must be > 0.");

  const double dt = tau_max / n_t;

  // Parse types and validate
  std::vector<bool> is_call(n_opt);
  for (int j = 0; j < n_opt; ++j) {
    const char* ts = CHAR(STRING_ELT(types, j));
    if (ts[0] == 'c' || ts[0] == 'C') is_call[j] = true;
    else if (ts[0] == 'p' || ts[0] == 'P') is_call[j] = false;
    else Rcpp::stop("type must be 'call' or 'put' at index %d.", j + 1);
  }

  // For each option, find which time step index corresponds to its tau.
  // tau_j maps to step index: n_t - round(tau_j / dt).
  // We extract the price at time step n where (n_t - n) * dt ≈ tau_j.
  // Equivalently, we extract at backward step (n_t - extraction_step[j]).
  std::vector<int> extract_step(n_opt); // backward step index at which to extract
  for (int j = 0; j < n_opt; ++j) {
    int steps_j = static_cast<int>(std::round(taus[j] / dt));
    steps_j = std::max(1, std::min(steps_j, n_t));
    extract_step[j] = steps_j;
  }

  // We need to use per-option r_d and r_f. Since CN coefficients depend on r_d and q(=r_f),
  // different rates mean different tridiagonal systems. Group by (r_d, r_f) for factorization.
  // For simplicity (and because in calibration all options share the same rates), we check
  // if rates are uniform. If so, use the fast shared-factorization path.
  // Otherwise, fall back to per-option solves (still zero-alloc).

  bool uniform_rates = true;
  double r_d0 = r_ds[0], r_f0 = r_fs[0];
  for (int j = 1; j < n_opt; ++j) {
    if (r_ds[j] != r_d0 || r_fs[j] != r_f0) { uniform_rates = false; break; }
  }

  // --- Pre-compute S_i values ---
  std::vector<double> S(n_s + 1);
  for (int i = 0; i <= n_s; ++i) S[i] = s_min + i * ds;

  // --- Allocate per-option price grids ---
  // old_prices[j] and new_prices[j] are (n_s+1) arrays for option j
  std::vector<std::vector<double>> old_p(n_opt, std::vector<double>(n_s + 1));
  std::vector<std::vector<double>> new_p(n_opt, std::vector<double>(n_s + 1));
  std::vector<bool> extracted(n_opt, false);

  // Result vector
  Rcpp::NumericVector result(n_opt);

  // --- Terminal condition (at tau_max) for each option ---
  for (int j = 0; j < n_opt; ++j) {
    double k = strikes[j];
    for (int i = 0; i <= n_s; ++i) {
      old_p[j][i] = is_call[j] ? std::max(S[i] - k, 0.0)
                                : std::max(k - S[i], 0.0);
    }
  }

  // --- Shared workspace for Thomas ---
  // Tridiagonal bands (interior nodes only, size N = n_s - 1)
  std::vector<double> a_vec(N), b_vec(N), c_vec(N);
  std::vector<double> cp(N), inv_denom(N);
  // Per-option RHS and solution (allocated once, reused)
  std::vector<double> d_vec(N), x_vec(N);

  // --- Backward time stepping ---
  for (int n = n_t - 1; n >= 0; --n) {
    double t_rem = (n_t - n) * dt; // time remaining at step n

    // Build tridiagonal LHS for this time step (shared across all options if uniform rates)
    if (uniform_rates) {
      double rd = r_d0, q = r_f0;
      for (int i = 1; i < n_s; ++i) {
        double sig_n = sigma(i, n);
        double alpha_n = 0.5 * sig_n * sig_n * S[i] * S[i];
        double diff_coeff = alpha_n / (ds * ds);
        double conv_coeff = (rd - q) * S[i] / (2.0 * ds);

        a_vec[i - 1] = -0.5 * dt * (diff_coeff - conv_coeff);
        b_vec[i - 1] =  1.0 + 0.5 * dt * (2.0 * diff_coeff + rd);
        c_vec[i - 1] = -0.5 * dt * (diff_coeff + conv_coeff);
      }
      thomas::factor(a_vec.data(), b_vec.data(), c_vec.data(),
                     cp.data(), inv_denom.data(), N);
    }

    // Process each active option
    for (int j = 0; j < n_opt; ++j) {
      if (extracted[j]) continue;

      double rd = r_ds[j], q = r_fs[j], k = strikes[j];
      double spot_j = spots[j];

      // Boundaries at this time step
      double u_bound, l_bound;
      if (is_call[j]) {
        u_bound = 0.0;
        l_bound = std::max(s_max * std::exp(-q * t_rem)
                           - k * std::exp(-rd * t_rem), 0.0);
      } else {
        u_bound = k * std::exp(-rd * t_rem);
        l_bound = 0.0;
      }

      new_p[j][0]    = u_bound;
      new_p[j][n_s]  = l_bound;

      // Build RHS (explicit part)
      for (int i = 1; i < n_s; ++i) {
        double sig_np1 = sigma(i, n + 1);
        double alpha_np = 0.5 * sig_np1 * sig_np1 * S[i] * S[i];
        double diff_coeff = alpha_np / (ds * ds);
        double conv_coeff = (rd - q) * S[i] / (2.0 * ds);

        d_vec[i - 1] =
            0.5 * dt * (diff_coeff - conv_coeff) * old_p[j][i - 1]
          + (1.0 - 0.5 * dt * (2.0 * diff_coeff + rd)) * old_p[j][i]
          + 0.5 * dt * (diff_coeff + conv_coeff) * old_p[j][i + 1];
      }

      // Boundary contributions
      if (uniform_rates) {
        d_vec[0]     -= a_vec[0]     * u_bound;
        d_vec[N - 1] -= c_vec[N - 1] * l_bound;
      } else {
        // Need to build LHS for this option's rates
        for (int i = 1; i < n_s; ++i) {
          double sig_n = sigma(i, n);
          double alpha_n = 0.5 * sig_n * sig_n * S[i] * S[i];
          double diff_coeff = alpha_n / (ds * ds);
          double conv_coeff = (rd - q) * S[i] / (2.0 * ds);

          a_vec[i - 1] = -0.5 * dt * (diff_coeff - conv_coeff);
          b_vec[i - 1] =  1.0 + 0.5 * dt * (2.0 * diff_coeff + rd);
          c_vec[i - 1] = -0.5 * dt * (diff_coeff + conv_coeff);
        }
        d_vec[0]     -= a_vec[0]     * u_bound;
        d_vec[N - 1] -= c_vec[N - 1] * l_bound;
      }

      // Zero the touching off-diagonals before solving
      double a0_save = a_vec[0], cN_save = c_vec[N - 1];

      // Solve
      if (uniform_rates) {
        // Use pre-factored solve (shared LHS)
        // Adjust d for zeroed off-diagonals: a[0] contribution already subtracted,
        // but factor used original a[0]. We need to adjust:
        // The factored solve expects the system as-is. Since a[0]=0 by convention
        // in the tridiagonal (boundary row removed), and we already subtracted
        // a_vec[0]*u_bound from d, we solve with factored a,b,c directly.
        thomas::solve_factored(a_vec.data(), cp.data(), inv_denom.data(),
                               d_vec.data(), x_vec.data(), N);
      } else {
        a_vec[0] = 0.0;
        c_vec[N - 1] = 0.0;
        thomas::solve(a_vec.data(), b_vec.data(), c_vec.data(),
                      d_vec.data(), x_vec.data(), cp.data(), N);
        a_vec[0] = a0_save;
        c_vec[N - 1] = cN_save;
      }

      for (int i = 1; i < n_s; ++i)
        new_p[j][i] = x_vec[i - 1];

      std::swap(old_p[j], new_p[j]);

      // Check if we should extract this option's price at this time step
      // extract_step[j] = number of backward steps for this option's tau
      // We are at backward step (n_t - n)
      int bstep = n_t - n;
      if (bstep == extract_step[j]) {
        result[j] = interp_raw(spot_j, s_min, ds, old_p[j].data(), n_s + 1);
        extracted[j] = true;
      }
    }
  }

  // Handle any options not yet extracted (shouldn't happen, but safety)
  for (int j = 0; j < n_opt; ++j) {
    if (!extracted[j]) {
      result[j] = interp_raw(spots[j], s_min, ds, old_p[j].data(), n_s + 1);
    }
  }

  return result;
}


// ============================================================================
// Batch American solver: prices multiple American options sharing one sigma
// surface, parallelized across options via OpenMP.
// ============================================================================

//' Batch-price American options under local volatility (1D)
//'
//' Prices multiple American options that share the same local volatility surface
//' and PDE grid. Each option is an independent backward Crank-Nicolson solve
//' with penalty-projection for early exercise, parallelized across options
//' via OpenMP.
//'
//' @param spots Numeric vector of spot prices (one per option).
//' @param strikes Numeric vector of strike prices.
//' @param taus Numeric vector of times to expiry (years).
//' @param r_ds Numeric vector of domestic risk-free rates.
//' @param r_fs Numeric vector of foreign risk-free rates or cost-of-carry rates.
//' @param sigma Local volatility matrix of size \code{(n_s + 1) x (n_t + 1)}.
//' @param types Character vector, each element \code{"call"} or \code{"put"}.
//' @param s_min,s_max Domain bounds for the asset grid.
//' @param n_s Number of spatial intervals (\code{n_s + 1} nodes).
//' @param n_t Number of time steps.
//' @param lambda Penalty parameter (> 0, typically 1e4 to 1e6).
//' @param tolerance Convergence tolerance for penalty iterations.
//'
//' @return Numeric vector of option prices (same length as \code{strikes}).
//'
//' @details
//' Each option is solved independently using the same sigma surface and grid.
//' The solver uses Crank-Nicolson with a penalty-projection fixed-point
//' iteration for the American constraint, capped at 200 iterations per time
//' step. Options are solved in parallel across CPU cores using std::thread.
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector batch_price_american_lv(
    const Rcpp::NumericVector& spots,
    const Rcpp::NumericVector& strikes,
    const Rcpp::NumericVector& taus,
    const Rcpp::NumericVector& r_ds,
    const Rcpp::NumericVector& r_fs,
    const Rcpp::NumericMatrix& sigma,
    const Rcpp::StringVector& types,
    double s_min, double s_max,
    int n_s, int n_t,
    double lambda, double tolerance) {

  const int n_opt = strikes.size();

  // --- Validation ---
  if (spots.size() != n_opt || taus.size() != n_opt || r_ds.size() != n_opt ||
      r_fs.size() != n_opt || types.size() != n_opt)
    Rcpp::stop("All input vectors must have the same length.");
  if (n_s < 2 || n_t < 1) Rcpp::stop("n_s >= 2 and n_t >= 1 required.");
  if (s_max <= s_min) Rcpp::stop("s_max must be > s_min.");
  if (lambda <= 0.0) Rcpp::stop("lambda must be > 0.");
  if (tolerance <= 0.0) Rcpp::stop("tolerance must be > 0.");
  if (sigma.nrow() != n_s + 1 || sigma.ncol() != n_t + 1)
    Rcpp::stop("sigma must have dimensions (n_s+1) x (n_t+1).");

  const double ds = (s_max - s_min) / n_s;
  const int N = n_s - 1;
  const int max_iter = 200;

  // Parse types (before parallel region — StringVector not thread-safe)
  std::vector<bool> is_call(n_opt);
  for (int j = 0; j < n_opt; ++j) {
    const char* ts = CHAR(STRING_ELT(types, j));
    if (ts[0] == 'c' || ts[0] == 'C') is_call[j] = true;
    else if (ts[0] == 'p' || ts[0] == 'P') is_call[j] = false;
    else Rcpp::stop("type must be 'call' or 'put' at index %d.", j + 1);
  }

  // Pre-compute S values
  std::vector<double> S(n_s + 1);
  for (int i = 0; i <= n_s; ++i) S[i] = s_min + i * ds;

  // Copy sigma to contiguous column-major raw array for thread-safe access
  // Rcpp::NumericMatrix is column-major: sigma(i,j) = data[i + j*(n_s+1)]
  const int sig_rows = n_s + 1;
  const int sig_cols = n_t + 1;
  std::vector<double> sig_raw(sig_rows * sig_cols);
  for (int j = 0; j < sig_cols; ++j)
    for (int i = 0; i < sig_rows; ++i)
      sig_raw[i + j * sig_rows] = sigma(i, j);

  // Copy scalar parameters for thread safety
  std::vector<double> spots_v(n_opt), strikes_v(n_opt), taus_v(n_opt),
                      rds_v(n_opt), rfs_v(n_opt);
  for (int j = 0; j < n_opt; ++j) {
    spots_v[j]   = spots[j];
    strikes_v[j] = strikes[j];
    taus_v[j]    = taus[j];
    rds_v[j]     = r_ds[j];
    rfs_v[j]     = r_fs[j];
  }

  std::vector<double> results(n_opt, 0.0);
  std::vector<int> warn_flags(n_opt, 0); // non-convergence warnings

  // Worker function: solves a range of options [j_begin, j_end)
  auto solve_range = [&](int j_begin, int j_end) {
    // Per-thread workspace (allocated once, reused across options)
    std::vector<double> old_prices(n_s + 1), new_prices(n_s + 1), payoff(n_s + 1);
    std::vector<double> a(N), b(N), c(N), d(N);
    std::vector<double> a0(N), b0(N), c0(N), d0(N);
    std::vector<double> bb(N), dd(N);
    std::vector<double> cp(N), x(N), temp(N), pv(N);

    for (int j = j_begin; j < j_end; ++j) {
      double k     = strikes_v[j];
      double tau_j = taus_v[j];
      double rd    = rds_v[j];
      double q     = rfs_v[j];
      double spot  = spots_v[j];
      bool   call  = is_call[j];

      if (tau_j <= 0.0) {
        results[j] = call ? std::max(spot - k, 0.0) : std::max(k - spot, 0.0);
        continue;
      }

      double dt_j = tau_j / n_t;

      // Terminal condition
      for (int i = 0; i <= n_s; ++i) {
        payoff[i] = call ? std::max(S[i] - k, 0.0) : std::max(k - S[i], 0.0);
        old_prices[i] = payoff[i];
      }

      // Backward time stepping
      for (int n = n_t - 1; n >= 0; --n) {
        double t_rem = (n_t - n) * dt_j;

        // Boundaries
        double u_bound, l_bound;
        if (call) {
          u_bound = 0.0;
          l_bound = std::max(s_max * std::exp(-q * t_rem)
                             - k * std::exp(-rd * t_rem), 0.0);
        } else {
          u_bound = k;       // American put at S=0: immediate exercise
          l_bound = 0.0;
        }

        new_prices[0]    = u_bound;
        new_prices[n_s]  = l_bound;

        // Build tridiagonal system
        // sigma column n for implicit level, n+1 for explicit level
        for (int i = 1; i < n_s; ++i) {
          double sig_n   = sig_raw[i + n * sig_rows];
          double sig_np1 = sig_raw[i + (n + 1) * sig_rows];

          double alpha_n  = 0.5 * sig_n * sig_n * S[i] * S[i];
          double alpha_np = 0.5 * sig_np1 * sig_np1 * S[i] * S[i];

          double diff_n  = alpha_n / (ds * ds);
          double conv    = (rd - q) * S[i] / (2.0 * ds);
          double diff_np = alpha_np / (ds * ds);

          a[i - 1] = -0.5 * dt_j * (diff_n - conv);
          b[i - 1] =  1.0 + 0.5 * dt_j * (2.0 * diff_n + rd);
          c[i - 1] = -0.5 * dt_j * (diff_n + conv);

          d[i - 1] =
              0.5 * dt_j * (diff_np - conv) * old_prices[i - 1]
            + (1.0 - 0.5 * dt_j * (2.0 * diff_np + rd)) * old_prices[i]
            + 0.5 * dt_j * (diff_np + conv) * old_prices[i + 1];
        }

        // Boundary contributions
        d[0]     -= a[0]     * u_bound;
        d[N - 1] -= c[N - 1] * l_bound;
        a[0]      = 0.0;
        c[N - 1]  = 0.0;

        // Save base system for penalty iteration
        std::copy(a.begin(), a.end(), a0.begin());
        std::copy(b.begin(), b.end(), b0.begin());
        std::copy(c.begin(), c.end(), c0.begin());
        std::copy(d.begin(), d.end(), d0.begin());

        // Initialize penalty iteration from previous time step
        for (int i = 0; i < N; ++i)
          temp[i] = old_prices[i + 1];

        double error = 1e30;
        int iter = 0;

        while (error > tolerance && iter < max_iter) {
          // Obstacle violation indicators
          for (int i = 0; i < N; ++i)
            pv[i] = (temp[i] < payoff[i + 1]) ? lambda : 0.0;

          // Rebuild system with penalty (from base, not accumulated)
          std::copy(b0.begin(), b0.end(), bb.begin());
          std::copy(d0.begin(), d0.end(), dd.begin());
          for (int i = 0; i < N; ++i) {
            bb[i] += pv[i];
            dd[i] += pv[i] * payoff[i + 1];
          }

          // Solve
          thomas::solve(a0.data(), bb.data(), c0.data(),
                        dd.data(), x.data(), cp.data(), N);

          // Convergence check
          double maxdiff = 0.0;
          for (int i = 0; i < N; ++i) {
            double diff = std::fabs(x[i] - temp[i]);
            if (diff > maxdiff) maxdiff = diff;
          }
          error = maxdiff;
          std::swap(x, temp); // temp now holds the latest solution
          ++iter;
        }

        if (iter >= max_iter) warn_flags[j] = n; // record which step failed

        // Write interior solution with obstacle projection
        for (int i = 1; i < n_s; ++i)
          new_prices[i] = std::max(temp[i - 1], payoff[i]);

        new_prices[0]    = u_bound;
        new_prices[n_s]  = l_bound;

        std::swap(old_prices, new_prices);
      }

      // Interpolate at spot
      results[j] = interp_raw(spot, s_min, ds, old_prices.data(), n_s + 1);
    }
  };

  // Distribute options across threads
  unsigned int n_threads = std::thread::hardware_concurrency();
  if (n_threads == 0) n_threads = 1;
  if (static_cast<int>(n_threads) > n_opt) n_threads = n_opt;

  if (n_threads <= 1) {
    // Single-threaded path
    solve_range(0, n_opt);
  } else {
    std::vector<std::thread> threads;
    threads.reserve(n_threads);
    int chunk = n_opt / n_threads;
    int remainder = n_opt % n_threads;
    int start = 0;
    for (unsigned int t = 0; t < n_threads; ++t) {
      int end = start + chunk + (static_cast<int>(t) < remainder ? 1 : 0);
      threads.emplace_back(solve_range, start, end);
      start = end;
    }
    for (auto& th : threads) th.join();
  }

  // Issue warnings outside parallel region
  for (int j = 0; j < n_opt; ++j) {
    if (warn_flags[j] > 0) {
      Rcpp::warning("Option %d: penalty iteration did not converge at time step %d.",
                    j + 1, warn_flags[j]);
    }
  }

  Rcpp::NumericVector out(n_opt);
  for (int j = 0; j < n_opt; ++j) out[j] = results[j];
  return out;
}
