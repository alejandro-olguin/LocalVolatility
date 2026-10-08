// [[Rcpp::interfaces(r, cpp)]]

#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <cmath>
#include <string>
#include <thread>
#include <functional>

#include "thomas_factored.h"
#include "sigma_time.h"
#include <limits>

using namespace Rcpp;

// ============================================================================
// Batch European solver: prices multiple European options sharing one sigma
// surface; one exact-step backward Crank-Nicolson sweep per option, threaded.
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
//' @param sigma Local volatility matrix with \code{n_s + 1} rows. Its columns
//'   are uniformly spaced in calendar time on \code{[0, t_max]}: column 1 is the
//'   valuation date and the last column is \code{t_max}. Any number of columns
//'   is accepted; values are linearly interpolated in time.
//' @param types Character vector, each element \code{"call"} or \code{"put"}.
//' @param s_min,s_max Domain bounds for the asset grid.
//' @param n_s Number of spatial intervals (\code{n_s + 1} nodes).
//' @param n_t Number of time steps per option (each option uses its own exact
//'   step \code{tau_j / n_t}, so maturities are never rounded).
//' @param t_max Calendar horizon spanned by the sigma columns (years). Defaults
//'   to \code{max(taus)} when not positive.
//' @param n_threads Number of threads for parallel option solves
//'   (0 = auto-detect; default 0).
//'
//' @return Numeric vector of option prices (same length as \code{strikes}).
//'
//' @details
//' All options must lie on the same spatial domain \code{[s_min, s_max]} and
//' use the same sigma surface. Each option is an independent backward
//' Crank-Nicolson sweep from its own maturity with step \code{tau_j / n_t};
//' at every time level the local volatility is read at the corresponding
//' calendar time by linear interpolation between the two neighbouring sigma
//' columns. Options are solved in parallel across CPU cores with
//' \code{std::thread}; a failed linear solve yields \code{NaN} with a warning
//' instead of aborting.
//'
//' @examples
//' Sigma <- matrix(0.2, nrow = 51, ncol = 51)
//' batch_price_european_lv(
//'   spots   = rep(100, 3),
//'   strikes = c(90, 100, 110),
//'   taus    = rep(1, 3),
//'   r_ds    = rep(0.05, 3),
//'   r_fs    = rep(0.02, 3),
//'   sigma   = Sigma,
//'   types   = rep("call", 3),
//'   s_min = 1, s_max = 400, n_s = 50, n_t = 50)
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
    int n_s, int n_t,
    double t_max = -1.0, int n_threads = 0) {

  const int n_opt = strikes.size();

  // --- Validation ---
  if (spots.size() != n_opt || taus.size() != n_opt || r_ds.size() != n_opt ||
      r_fs.size() != n_opt || types.size() != n_opt)
    Rcpp::stop("All input vectors must have the same length.");
  if (n_s < 2 || n_t < 1) Rcpp::stop("n_s >= 2 and n_t >= 1 required.");
  if (s_max <= s_min) Rcpp::stop("s_max must be > s_min.");
  if (n_threads < 0) Rcpp::stop("n_threads must be >= 0.");
  if (sigma.nrow() != n_s + 1 || sigma.ncol() < 1)
    Rcpp::stop("sigma must have n_s+1 rows (columns = calendar-time nodes on [0, t_max]).");

  const double ds = (s_max - s_min) / n_s;
  const int N = n_s - 1; // interior nodes

  // Calendar horizon spanned by the sigma columns (defaults to max(taus)).
  double tau_max = 0.0;
  for (int j = 0; j < n_opt; ++j)
    if (taus[j] > tau_max) tau_max = taus[j];
  if (tau_max <= 0.0) Rcpp::stop("All taus must be > 0.");
  if (t_max <= 0.0) t_max = tau_max;

  // Parse types and validate (before the parallel region)
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

  // Copy sigma to a contiguous column-major raw array (thread-safe access)
  const int sig_rows = n_s + 1;
  const int sig_cols = sigma.ncol();
  std::vector<double> sig_raw(static_cast<std::size_t>(sig_rows) * sig_cols);
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
  std::vector<int>    err_flags(n_opt, 0);

  // Worker: each option is an independent backward Crank-Nicolson sweep with
  // its own exact step dt_j = tau_j / n_t; sigma is read at calendar time.
  auto solve_range = [&](int j_begin, int j_end) {
    std::vector<double> old_prices(n_s + 1), new_prices(n_s + 1);
    std::vector<double> a(N), b(N), c(N), d(N), cp(N), x(N);

    for (int j = j_begin; j < j_end; ++j) {
      try {
        const double k     = strikes_v[j];
        const double tau_j = taus_v[j];
        const double rd    = rds_v[j];
        const double q     = rfs_v[j];
        const double spot  = spots_v[j];
        const bool   call  = is_call[j];

        if (tau_j <= 0.0) {
          results[j] = call ? std::max(spot - k, 0.0) : std::max(k - spot, 0.0);
          continue;
        }

        const double dt_j = tau_j / n_t;

        for (int i = 0; i <= n_s; ++i)
          old_prices[i] = call ? std::max(S[i] - k, 0.0) : std::max(k - S[i], 0.0);

        for (int n = n_t - 1; n >= 0; --n) {
          const double t_rem = (n_t - n) * dt_j;
          const lvtime::ColLerp c_imp = lvtime::col_lerp(n * dt_j,       t_max, sig_cols);
          const lvtime::ColLerp c_exp = lvtime::col_lerp((n + 1) * dt_j, t_max, sig_cols);

          double u_bound, l_bound;
          if (call) {
            u_bound = 0.0;
            l_bound = std::max(s_max * std::exp(-q * t_rem)
                               - k * std::exp(-rd * t_rem), 0.0);
          } else {
            u_bound = k * std::exp(-rd * t_rem);
            l_bound = 0.0;
          }

          for (int i = 1; i < n_s; ++i) {
            const double sig_n   = lvtime::at_raw(sig_raw.data(), sig_rows, i, c_imp);
            const double sig_np1 = lvtime::at_raw(sig_raw.data(), sig_rows, i, c_exp);
            const double diff_n  = 0.5 * sig_n   * sig_n   * S[i] * S[i] / (ds * ds);
            const double diff_np = 0.5 * sig_np1 * sig_np1 * S[i] * S[i] / (ds * ds);
            const double conv    = (rd - q) * S[i] / (2.0 * ds);

            a[i - 1] = -0.5 * dt_j * (diff_n - conv);
            b[i - 1] =  1.0 + 0.5 * dt_j * (2.0 * diff_n + rd);
            c[i - 1] = -0.5 * dt_j * (diff_n + conv);

            d[i - 1] =
                0.5 * dt_j * (diff_np - conv) * old_prices[i - 1]
              + (1.0 - 0.5 * dt_j * (2.0 * diff_np + rd)) * old_prices[i]
              + 0.5 * dt_j * (diff_np + conv) * old_prices[i + 1];
          }

          d[0]     -= a[0]     * u_bound;
          d[N - 1] -= c[N - 1] * l_bound;
          a[0]      = 0.0;
          c[N - 1]  = 0.0;

          thomas::solve(a.data(), b.data(), c.data(), d.data(), x.data(), cp.data(), N);

          new_prices[0]   = u_bound;
          new_prices[n_s] = l_bound;
          for (int i = 1; i < n_s; ++i) new_prices[i] = x[i - 1];
          std::swap(old_prices, new_prices);
        }

        results[j] = interp_raw(spot, s_min, ds, old_prices.data(), n_s + 1);
      } catch (...) {
        // Never let an exception escape a worker thread (std::terminate).
        results[j]   = std::numeric_limits<double>::quiet_NaN();
        err_flags[j] = 1;
      }
    }
  };

  unsigned int n_thr = (n_threads > 0)
    ? static_cast<unsigned int>(n_threads)
    : std::thread::hardware_concurrency();
  if (n_thr == 0) n_thr = 1;
  if (static_cast<int>(n_thr) > n_opt) n_thr = n_opt;

  if (n_thr <= 1) {
    solve_range(0, n_opt);
  } else {
    std::vector<std::thread> threads;
    threads.reserve(n_thr);
    const int chunk = n_opt / n_thr, remainder = n_opt % n_thr;
    int start = 0;
    for (unsigned int t = 0; t < n_thr; ++t) {
      const int end = start + chunk + (static_cast<int>(t) < remainder ? 1 : 0);
      threads.emplace_back(solve_range, start, end);
      start = end;
    }
    for (auto& th : threads) th.join();
  }

  for (int j = 0; j < n_opt; ++j)
    if (err_flags[j])
      Rcpp::warning("Option %d: linear solve failed (zero pivot); price set to NaN.", j + 1);

  Rcpp::NumericVector result(n_opt);
  for (int j = 0; j < n_opt; ++j) result[j] = results[j];
  return result;
}


// ============================================================================
// Batch American solver: prices multiple American options sharing one sigma
// surface, parallelized across options via std::thread.
// ============================================================================

//' Batch-price American options under local volatility (1D)
//'
//' Prices multiple American options that share the same local volatility surface
//' and PDE grid. Each option is an independent backward Crank-Nicolson solve
//' with penalty-projection for early exercise, parallelized across options
//' via \code{std::thread}.
//'
//' @param spots Numeric vector of spot prices (one per option).
//' @param strikes Numeric vector of strike prices.
//' @param taus Numeric vector of times to expiry (years).
//' @param r_ds Numeric vector of domestic risk-free rates.
//' @param r_fs Numeric vector of foreign risk-free rates or cost-of-carry rates.
//' @param sigma Local volatility matrix with \code{n_s + 1} rows. Its columns
//'   are uniformly spaced in calendar time on \code{[0, t_max]}: column 1 is the
//'   valuation date and the last column is \code{t_max}. Any number of columns
//'   is accepted; values are linearly interpolated in time.
//' @param types Character vector, each element \code{"call"} or \code{"put"}.
//' @param s_min,s_max Domain bounds for the asset grid.
//' @param n_s Number of spatial intervals (\code{n_s + 1} nodes).
//' @param n_t Number of time steps per option (exact step \code{tau_j / n_t}).
//' @param lambda Penalty parameter (> 0, typically 1e4 to 1e6).
//' @param tolerance Convergence tolerance for penalty iterations.
//' @param n_threads Number of threads for parallel option solves
//'   (0 = auto-detect; default 0).
//' @param t_max Calendar horizon spanned by the sigma columns (years). Defaults
//'   to \code{max(taus)} when not positive.
//'
//' @return Numeric vector of option prices (same length as \code{strikes}).
//'
//' @details
//' Each option is solved independently using the same sigma surface and grid;
//' the local volatility is read at calendar time (see \code{sigma}).
//' The solver uses Crank-Nicolson with a penalty-projection fixed-point
//' iteration for the American constraint, capped at 200 iterations per time
//' step. Options are solved in parallel across CPU cores using \code{std::thread}.
//' For tidy-data workflows, bind the returned price vector to contract metadata
//' so each row corresponds to one option and each column to one variable.
//'
//' @examples
//' Sigma <- matrix(0.2, nrow = 51, ncol = 51)
//' batch_price_american_lv(
//'   spots   = rep(100, 3),
//'   strikes = c(90, 100, 110),
//'   taus    = rep(1, 3),
//'   r_ds    = rep(0.05, 3),
//'   r_fs    = rep(0.02, 3),
//'   sigma   = Sigma,
//'   types   = rep("put", 3),
//'   s_min = 1, s_max = 400, n_s = 50, n_t = 50,
//'   lambda = 1e4, tolerance = 1e-8,
//'   n_threads = 0)
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
    double lambda, double tolerance, int n_threads = 0, double t_max = -1.0) {

  const int n_opt = strikes.size();

  // --- Validation ---
  if (spots.size() != n_opt || taus.size() != n_opt || r_ds.size() != n_opt ||
      r_fs.size() != n_opt || types.size() != n_opt)
    Rcpp::stop("All input vectors must have the same length.");
  if (n_s < 2 || n_t < 1) Rcpp::stop("n_s >= 2 and n_t >= 1 required.");
  if (s_max <= s_min) Rcpp::stop("s_max must be > s_min.");
  if (lambda <= 0.0) Rcpp::stop("lambda must be > 0.");
  if (tolerance <= 0.0) Rcpp::stop("tolerance must be > 0.");
  if (n_threads < 0) Rcpp::stop("n_threads must be >= 0.");
  if (sigma.nrow() != n_s + 1 || sigma.ncol() < 1)
    Rcpp::stop("sigma must have n_s+1 rows (columns = calendar-time nodes on [0, t_max]).");

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
  const int sig_cols = sigma.ncol();
  std::vector<double> sig_raw(static_cast<std::size_t>(sig_rows) * sig_cols);
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

  double tau_max = 0.0;
  for (int j = 0; j < n_opt; ++j) if (taus_v[j] > tau_max) tau_max = taus_v[j];
  if (t_max <= 0.0) t_max = tau_max;

  std::vector<double> results(n_opt, 0.0);
  std::vector<int> warn_flags(n_opt, -1); // time step of first non-convergence
  std::vector<int> err_flags(n_opt, 0);

  // Worker function: solves a range of options [j_begin, j_end)
  auto solve_range = [&](int j_begin, int j_end) {
    // Per-thread workspace (allocated once, reused across options)
    std::vector<double> old_prices(n_s + 1), new_prices(n_s + 1), payoff(n_s + 1);
    std::vector<double> a(N), b(N), c(N), d(N);
    std::vector<double> a0(N), b0(N), c0(N), d0(N);
    std::vector<double> bb(N), dd(N);
    std::vector<double> cp(N), x(N), temp(N), pv(N);

    for (int j = j_begin; j < j_end; ++j) {
     try {
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
        const lvtime::ColLerp c_imp = lvtime::col_lerp(n * dt_j,       t_max, sig_cols);
        const lvtime::ColLerp c_exp = lvtime::col_lerp((n + 1) * dt_j, t_max, sig_cols);

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

        // Build tridiagonal system; sigma read at calendar time t_n (implicit)
        // and t_{n+1} (explicit) by interpolation between columns.
        for (int i = 1; i < n_s; ++i) {
          double sig_n   = lvtime::at_raw(sig_raw.data(), sig_rows, i, c_imp);
          double sig_np1 = lvtime::at_raw(sig_raw.data(), sig_rows, i, c_exp);

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
     } catch (...) {
      results[j]   = std::numeric_limits<double>::quiet_NaN();
      err_flags[j] = 1;
     }
    }
  };

  // Distribute options across threads
  unsigned int n_thr = (n_threads > 0)
    ? static_cast<unsigned int>(n_threads)
    : std::thread::hardware_concurrency();
  if (n_thr == 0) n_thr = 1;
  if (static_cast<int>(n_thr) > n_opt) n_thr = n_opt;

  if (n_thr <= 1) {
    // Single-threaded path
    solve_range(0, n_opt);
  } else {
    std::vector<std::thread> threads;
    threads.reserve(n_thr);
    int chunk = n_opt / n_thr;
    int remainder = n_opt % n_thr;
    int start = 0;
    for (unsigned int t = 0; t < n_thr; ++t) {
      int end = start + chunk + (static_cast<int>(t) < remainder ? 1 : 0);
      threads.emplace_back(solve_range, start, end);
      start = end;
    }
    for (auto& th : threads) th.join();
  }

  // Issue warnings outside parallel region
  for (int j = 0; j < n_opt; ++j) {
    if (err_flags[j])
      Rcpp::warning("Option %d: linear solve failed (zero pivot); price set to NaN.", j + 1);
    else if (warn_flags[j] >= 0)
      Rcpp::warning("Option %d: penalty iteration did not converge at time step %d.",
                    j + 1, warn_flags[j]);
  }

  Rcpp::NumericVector out(n_opt);
  for (int j = 0; j < n_opt; ++j) out[j] = results[j];
  return out;
}
