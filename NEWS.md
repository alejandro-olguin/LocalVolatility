# LocalVolatility NEWS

## 3.2.0 (2026-04-14)

### Documentation

- **Vignette** (`vignettes/local-volatility.Rmd`): Added a comprehensive
  JSS-standard vignette with the mathematical framework, software design
  overview, and fully reproducible examples for all solver families.
- **Package-level help page** (`R/LocalVolatility-package.R`): Added a
  structured function index with sections for 1D, 2D, Heston, and Monte Carlo
  solvers, plus full reference list.
- **`inst/CITATION`**: Added formal bibentry so users can cite the package.
- **Runnable examples**: Added `@examples` blocks to all 15 exported functions.

### Robustness

- **2 GB memory guard** in `mc_american_heston_4d()`: Added a pre-allocation
  check that stops with a clear error if `n_paths * n_steps` would exceed 2 GB,
  preventing silent out-of-memory crashes.

### Code quality

- Removed an unused variable warning from the MC American solver.
- Removed `.DS_Store` files and R CMD check build artifacts from tracking
  (added patterns to `.Rbuildignore` and `.gitignore`).

## 3.1.0 (2026-04-06)

### New features

- **Heston stochastic volatility solvers**: Three new exported functions for
  pricing European and American options under the Heston (1993) model:
  - `heston_cf()`: Closed-form European price via the Heston characteristic
    function with Simpson's rule numerical integration. Uses the numerically
    stable "little Heston trap" formulation (Albrecher et al. 2007).
  - `european_option_heston()`: 2D finite-difference PDE solver on (S, v)
    using Yanenko operator splitting with Tavella-Randall grids.
  - `american_option_heston()`: Same PDE solver with penalty-projection
    for early exercise.
- **4D Monte Carlo (double-Heston quanto)**: Two new exported functions for
  pricing options on a foreign equity with stochastic volatility on both
  equity and FX:
  - `mc_european_heston_4d()`: European MC with `std::thread` parallelism,
    antithetic variates, Euler-Maruyama with full truncation, and Cholesky
    decomposition for correlated 4D Brownian motion.
  - `mc_american_heston_4d()`: American MC via Longstaff-Schwartz
    least-squares regression with parallel forward simulation.
- All Heston/MC functions validate inputs and warn when the Feller condition
  (2κθ > ξ²) is violated.

### Testing

- Added `test_heston.R` with 26 tests: CF put-call parity, positivity,
  deep ITM/OTM, input validation, Feller warning, PDE vs CF convergence
  (call/put), PDE put-call parity, American ≥ European.
- Added `test_mc_heston_4d.R` with 18 tests: positivity, reproducibility,
  put-call parity (within MC error), MC vs 2D PDE convergence (xi→0),
  SE scaling with n_paths, American ≥ European, input validation.
- Total test count: 115 tests across 8 files.

## 3.0.0 (2026-04-05)

### BREAKING CHANGES — Parameter and file renames

This release homogenises naming across all exported functions. **Users must
update their code** as follows:

| Function | Old parameter | New parameter |
|----------|--------------|---------------|
| `european_option_lv` | `q` | `r_f` |
| `american_option_lv` | `q` | `r_f` |
| `batch_price_american_lv` | `qs` | `r_fs` |
| `european_option_cf_2d` | `n` (ADR ratio) | `adr_ratio` |
| `european_option_cf_2d` | order: `..., rho, n, type` | order: `..., rho, type, adr_ratio` |
| `tavella_randall` | `x0` | `x_0` |
| `tavella_randall` | `n` | `n_grid` |

**Migration guide:**
- In 1D calls, replace `q = ...` with `r_f = ...` (same PDE slot: cost-of-carry).
- In `batch_price_american_lv`, replace `qs = ...` with `r_fs = ...`.
- In `european_option_cf_2d`, move `type` before `adr_ratio` and rename `n` to `adr_ratio`.
- In `tavella_randall`, rename `x0` to `x_0` and `n` to `n_grid`.

### Naming improvements
- Source files renamed from `1D_*` / `2D_*` prefixes to consistent
  `{exercise}_{variant}_{dim}.cpp` pattern (e.g. `european_lv_1d.cpp`).
- Standardised `@param tolerance` description across all American solvers.
- Standardised `@param r_f` / `@param r_fs` descriptions across 1D and batch solvers.

## 2.2.0 (2026-04-01)

### New features

- **Batch European solver** (`batch_price_european_lv`): Prices multiple European
  options in a single C++ call sharing one sigma surface. Uses pre-factored Thomas
  algorithm — the tridiagonal LHS is factored once per time step, then each option
  only needs an O(N) forward/back substitution. Supports mixed strikes, maturities,
  rates, and option types.
- **Batch American solver** (`batch_price_american_lv`): Prices multiple American
  options in a single C++ call with `std::thread` parallelism. Each thread gets a
  pre-allocated workspace (tridiagonal bands, penalty buffers) with zero R heap
  allocations in the hot path.
- **Zero-allocation Thomas algorithm** (`src/thomas_factored.h`): Header-only
  routines operating on raw `double*`. Three variants: `thomas::solve()` (full
  solve), `thomas::factor()` (LU factorization), `thomas::solve_factored()`
  (substitution with pre-factored LHS).

### Performance

- Batch European: ~20-50x faster than serial loop (shared factorization + extraction
  at per-option tau).
- Batch American: ~5-8x faster than serial loop (thread parallelism + zero-alloc
  Thomas). Measured 5.2x speedup on gradient evaluation with 25 options.

## 2.1.0 (2026-04-01)

### Bug fixes

- **Cross-derivative denominator in Pass 2 (all 2D solvers):** Fixed a copy-paste
  error where the mixed-derivative stencil denominator in the X-implicit pass
  duplicated two terms instead of using cross-products. On non-uniform
  (Tavella-Randall) grids this produced pricing error proportional to grid
  non-uniformity times rho. Uniform grids were unaffected.
- **1D interpolation with non-zero s_min:** Both 1D solvers called the wrong
  `linear_interpolation` overload, defaulting `s_min = 0`. Fixed to pass the
  user-supplied `s_min`, which matters whenever the grid does not start at zero.
- **1D American call far-field boundary:** Was a constant undiscounted
  `max(s_max - k, 0)` for all time steps. Changed to the time-dependent
  discounted European asymptote; the penalty method handles early exercise
  separately.
- **1D American put boundary at S = 0:** Was `K * exp(-r_d * tau_rem)` (European
  put value). Changed to `K` (immediate exercise value), which is the correct
  Dirichlet condition for an American put at the lower boundary.
- **Missing `exp(-r_f * tau)` in 2D far-field boundaries:** All six 2D solvers
  used `S * X * exp(-q * tau)` on the product leg. Under the domestic measure the
  correct asymptote is `S * X * exp(-(r_f + q) * tau) - K * exp(-r_d * tau)`.
  The missing foreign-rate discount could produce boundary error proportional to
  `1 - exp(-r_f * tau)`.

### Robustness

- Penalty iterations now have a hard cap of 200 with an `Rcpp::warning()` on
  non-convergence (all three American solvers).
- Added sigma-matrix dimension validation: 1D solvers check `(n_s+1) x (n_t+1)`;
  2D local-vol solvers check `(n_s+1) x n_t` and `(n_x+1) x n_t`.
- Added `rho in [-1, 1]` validation to all 2D PDE solvers (previously only the
  closed-form function checked).
- Added `sigma >= 0`, `lambda > 0`, and `tolerance > 0` guards where applicable.

### Code quality

- Homogenised all comments and error messages to English (previously mixed
  Spanish/English in the 1D American solver).
- Replaced `U = clone(Unew)` with `std::swap` in European 2D solvers.
- Removed unused `iters` counter from 2D American solvers.
- Rewrote cross-derivative denominator in factored form
  `(h_im1 + h_i) * (hx[j-1] + hx[j])` for clarity.
- Consistent code style, indentation, and variable naming across all eight
  source files.

### Documentation

- Sigma-matrix size requirements are now precise in all roxygen blocks
  (no longer "commonly ... or ...").
- Added `@return` and `@details` with LaTeX formulae to every exported function.
- Added `@export` to `american_option_2d` (was missing).
- Updated README with corrected boundary formulae and expanded testing section.

### Testing

- Added `test_batch_solvers.R` with 11 tests for batch European and American solvers:
  correctness vs single solver (calls, puts, mixed taus, mixed types), positivity,
  input validation (mismatched lengths, wrong sigma dimensions, invalid lambda).
- Expanded from 1 test file / 2 assertions to 6 test files / 57 tests:
  - `test_1d_european.R` -- Black-Scholes convergence (call/put), put-call
    parity, non-zero s_min, sigma-dimension rejection.
  - `test_1d_american.R` -- American >= European (put), non-dividend call
    equals European, input validation.
  - `test_closed_form.R` -- Put-call parity, tau = 0 intrinsic, sigma = 0 edge
    case, correlation extremes, input rejection.
  - `test_tavella_randall.R` -- Endpoint exactness, monotonicity, node
    concentration near x0, input rejection.
  - `test_2d_local_equals_constant.R` -- Extended with put options, high-r_f
    PDE vs closed-form comparison, American finite/positive checks.

## 2.0.0 (2025-10-24)

- 2D PDE drift and boundary alignment
  - Adopted domestic-measure S-drift consistently across 2D solvers:
    \(\mu_S = r_f - q - \rho\,\sigma_S\,\sigma_X\).
  - Updated far-field Dirichlet boundaries for product options to use
    \(\exp(-q\,\tau)\) on the \(S\cdot X\) term and \(\exp(-r_d\,\tau)\) on \(K\).
- American option solvers (2D) use penalty projection with split penalty per sub-step for stability.
- Added equivalence tests ensuring local-vol matrices filled with constants reproduce constant-vol prices (2D, European and American).
- Documentation overhaul: detailed roxygen blocks colocated with exported C++ functions; clarified grid sizes, drifts, boundaries, and examples.
- Minor typo fixes in parameter docs; standardized parameter names.

## 1.0 (2024-11-04)

- Initial release with pricing engines for 1D/2D options under local volatility.
- Utilities: Tavella-Randall nonuniform grid, Thomas algorithm, bilinear interpolation.
