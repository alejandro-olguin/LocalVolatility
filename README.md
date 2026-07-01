# LocalVolatility

Pricing engines for options under local and stochastic volatility, implemented in C++ via Rcpp and exposed to R. Includes European and American exercise, operator splitting, and robust boundary conditions.

- Fast finite-difference solvers (Crank-Nicolson in 1D; Yanenko splitting in 2D)
- **Heston stochastic volatility**: closed-form (characteristic function) and PDE solvers for European and American options
- **4D Monte Carlo** (double-Heston quanto): European and American (Longstaff-Schwartz) with `std::thread` parallelism, antithetic variates, and full 4×4 correlation via Cholesky
- Nonuniform Tavella-Randall grids
- Mixed derivative with correlation in 2D
- Penalty projection for American options
- Local volatility via space-time matrices; constant-vol 2D variant for parity checks
- Closed-form Black-Scholes for ADR products (2D)
- **Batch solvers** for calibration workloads: price many options in a single C++ call with shared grid setup, pre-factored Thomas algorithm, and `std::thread` parallelism

## Installation

This is an R package. From an R session:

```r
install.packages("devtools")     # if needed
devtools::install()               # from this repository folder
# or, for development:
devtools::load_all()
```

macOS users need a working C++ toolchain for Rcpp (e.g., Xcode + Command Line Tools and the proper R toolchain).

## Quick start

### 1D European (local volatility)

```r
library(LocalVolatility)

s_0 <- 100; k <- 100; tau <- 1
r_d <- 0.05; r_f <- 0.02
n_s <- 200; n_t <- 200
s_min <- 1; s_max <- 400

# Sigma matrix: (n_s + 1) x (n_t + 1)
Sigma <- matrix(0.2, nrow = n_s + 1, ncol = n_t + 1)
price_eur_1d <- european_option_lv(s_0, k, tau, r_d, r_f, Sigma, "call",
                                   s_min, s_max, n_s, n_t)
price_eur_1d
```

### 1D American (local volatility)

```r
lambda <- 1e4; tol <- 1e-8
price_am_1d <- american_option_lv(s_0, k, tau, r_d, r_f, Sigma, "put",
                                  s_min, s_max, n_s, n_t, lambda, tol)
price_am_1d
```

### 2D European (constant volatility)

```r
s_0 <- 100; x_0 <- 20; k <- 2000; tau <- 1
r_d <- 0.05; r_f <- 0.02; q <- 0.01; rho <- 0.3
sigma_s <- 0.2; sigma_x <- 0.15
s_min <- 10; s_max <- 300
x_min <-  1; x_max <-  60
n_s <- 40; n_x <- 40; n_t <- 50; alpha <- 3

price_eur_2d_cv <- european_option_2d(s_0, x_0, k, tau, r_d, r_f, q,
                                      sigma_s, sigma_x, rho, "call",
                                      s_min, s_max, x_min, x_max,
                                      n_s, n_x, n_t, alpha)
price_eur_2d_cv
```

### 2D European (local volatility)

```r
# Sigma matrices: (n_s + 1) x n_t and (n_x + 1) x n_t
SigmaS <- matrix(sigma_s, nrow = n_s + 1, ncol = n_t)
SigmaX <- matrix(sigma_x, nrow = n_x + 1, ncol = n_t)

price_eur_2d_lv <- european_option_lv_2d(s_0, x_0, k, tau, r_d, r_f, q,
                                         SigmaS, SigmaX, rho, "call",
                                         s_min, s_max, x_min, x_max,
                                         n_s, n_x, n_t, alpha)
price_eur_2d_lv
```

### 2D American (constant vs local volatility)

```r
lambda <- 1e4; tol <- 1e-8
price_am_2d_cv <- american_option_2d(s_0, x_0, k, tau, r_d, r_f, q,
                                     sigma_s, sigma_x, rho, "call",
                                     s_min, s_max, x_min, x_max,
                                     n_s, n_x, n_t, alpha, lambda, tol)

price_am_2d_lv <- american_option_lv_2d(s_0, x_0, k, tau, r_d, r_f, q,
                                        SigmaS, SigmaX, rho, "call",
                                        s_min, s_max, x_min, x_max,
                                        n_s, n_x, n_t, alpha, lambda, tol)

c(constant = price_am_2d_cv, local = price_am_2d_lv)
```

### Batch pricing (calibration)

```r
# Price multiple European options in one call
spots  <- rep(100, 5)
strikes <- c(90, 95, 100, 105, 110)
taus   <- rep(1.0, 5)
r_ds   <- rep(0.05, 5); r_fs <- rep(0.02, 5)
types  <- rep("call", 5)
Sigma  <- matrix(0.2, nrow = n_s + 1, ncol = n_t + 1)

prices <- batch_price_european_lv(spots, strikes, taus, r_ds, r_fs,
                                   Sigma, types, s_min, s_max, n_s, n_t)

# Price multiple American options in one call (parallel)
prices_am <- batch_price_american_lv(spots, strikes, taus, r_ds, r_fs,
                                      Sigma, types, s_min, s_max, n_s, n_t,
                                      lambda = 1e4, tolerance = 1e-8,
                                      n_threads = 0)
```

### Heston stochastic volatility

```r
# Heston parameters
kappa <- 2; theta <- 0.04; xi <- 0.5; rho_h <- -0.7; v_0 <- 0.04

# Closed-form (characteristic function)
price_heston_cf <- heston_cf(100, 100, 1, 0.05, 0, kappa, theta, xi, rho_h, v_0, "call")

# European PDE
price_heston_pde <- european_option_heston(
  100, 100, 1, 0.05, 0, kappa, theta, xi, rho_h, v_0, "call",
  s_min = 20, s_max = 300, v_min = 0.001, v_max = 1.0,
  n_s = 80, n_v = 40, n_t = 100, alpha = 3)

# American PDE
price_heston_am <- american_option_heston(
  100, 100, 1, 0.05, 0, kappa, theta, xi, rho_h, v_0, "put",
  s_min = 20, s_max = 300, v_min = 0.001, v_max = 1.0,
  n_s = 80, n_v = 40, n_t = 100, alpha = 3, lambda = 1e4, tolerance = 1e-8)

c(cf = price_heston_cf, pde = price_heston_pde, american = price_heston_am)
```

### 4D Monte Carlo (double-Heston quanto)

```r
# European MC: equity and FX each with Heston stochastic vol
mc_eu <- mc_european_heston_4d(
  s_0 = 100, x_0 = 20, k = 2000, tau = 1,
  r_d = 0.05, r_f = 0.02, q = 0.01,
  kappa_s = 2, theta_s = 0.04, xi_s = 0.5, v_s0 = 0.04,    # equity Heston
  kappa_x = 1.5, theta_x = 0.02, xi_x = 0.3, v_x0 = 0.02,  # FX Heston
  rho_sx = 0.3, rho_sv = -0.7, rho_xv = -0.5,               # key correlations
  rho_svx = 0, rho_xvs = 0, rho_vsvx = 0,                   # cross (usually ~0)
  type = "call", n_paths = 500000, n_steps = 100, seed = 42)
c(price = mc_eu$price, se = mc_eu$std_error)

# American MC (Longstaff-Schwartz)
mc_am <- mc_american_heston_4d(
  s_0 = 100, x_0 = 20, k = 2000, tau = 1,
  r_d = 0.05, r_f = 0.02, q = 0.01,
  kappa_s = 2, theta_s = 0.04, xi_s = 0.5, v_s0 = 0.04,
  kappa_x = 1.5, theta_x = 0.02, xi_x = 0.3, v_x0 = 0.02,
  rho_sx = 0.3, rho_sv = -0.7, rho_xv = -0.5,
  rho_svx = 0, rho_xvs = 0, rho_vsvx = 0,
  type = "put", n_paths = 100000, n_steps = 50, seed = 42, n_basis = 4)
c(price = mc_am$price, se = mc_am$std_error)
```

## Tidy-data conventions

The package is designed to fit tidy workflows:

- Vectorised batch APIs (`batch_price_european_lv()`, `batch_price_american_lv()`)
  return one value per input contract, preserving input order.
- Build analysis tables with one row per contract and one column per variable
  (for example: `spot`, `strike`, `tau`, `type`, `model_price`).
- Scalar/list outputs can be wrapped as one-row records with explicit fields.

```r
contracts <- data.frame(
  spot = spots,
  strike = strikes,
  tau = taus,
  r_d = r_ds,
  r_f = r_fs,
  type = types
)
contracts$model_price <- prices
head(contracts)
```

### Closed-form (ADR)

```r
price_cf <- european_option_cf_2d(s_0, x_0, k, tau, r_d, r_f, q,
                                   sigma_s, sigma_x, rho, "call", 1)
price_cf
```

## Mathematics

Let \(\tau = T - t\) denote time-to-maturity; prices are computed backward in time.

### 1D European (local volatility)

Let \(v = v(S,t)\). Under the domestic money-market measure:

\[ v_t + \tfrac{1}{2}\,\sigma(S,t)^2\,S^2\,v_{SS} + (r_d - r_f) S v_S - r_d v = 0, \quad v(S,T) = \max(\pm(S-K),0). \]

Boundary asymptotes for calls: \(S\to 0: v\to 0\); \(S\to \infty: v\to S e^{-r_f\tau} - K e^{-r_d\tau}\)_+.

### 1D American (local volatility)

Linear complementarity problem (LCP): with payoff \(\phi\),

\[ \min\{ -\mathcal{L} v,\ v - \phi \} = 0,\quad \mathcal{L}v := v_t + \tfrac{1}{2}\sigma^2 S^2 v_{SS} + (r_d-r_f) S v_S - r_d v. \]

We solve via a penalty-projection scheme with iteration cap and convergence monitoring.

### 2D European (domestic measure)

Let \(V = V(S,X,t)\). With equity local vol \(\sigma_S\) and FX local vol \(\sigma_X\), correlation \(\rho\),

\[ V_t + \tfrac{1}{2}(\sigma_S S)^2 V_{SS} + \tfrac{1}{2}(\sigma_X X)^2 V_{XX} + \rho\,\sigma_S\,\sigma_X\,S X\, V_{SX}
  + \underbrace{\mu_S}_{r_f - q - \rho\sigma_S\sigma_X} S V_S + (r_d - r_f) X V_X - r_d V = 0. \]

Terminal condition: \(V(S,X,T) = \max(\pm(SX-K),0)\). Far-field Dirichlet boundaries use
\( (SX)^+ e^{-(r_f + q)\tau} - K e^{-r_d\tau} \) for calls and the corresponding put asymptotes.

### 2D American

Same differential operator, with the LCP \(\min\{-\mathcal{L}V, V-\phi\}=0\). We use Yanenko operator
splitting (implicit in one dimension per sub-step) and a penalty projection; the penalty is split across
sub-steps for stability.

### Heston stochastic volatility

Let \(V = V(S,v,t)\) where \(v\) is the instantaneous variance. The Heston (1993) PDE is:

\[ V_t + \tfrac{1}{2}v S^2 V_{SS} + \rho\xi v S V_{Sv} + \tfrac{1}{2}\xi^2 v V_{vv} + (r_d - q) S V_S + \kappa(\theta - v) V_v - r_d V = 0. \]

Parameters: \(\kappa\) (mean reversion speed), \(\theta\) (long-run variance), \(\xi\) (vol-of-vol),
\(\rho\) (stock-variance correlation), \(v_0\) (initial variance). Terminal condition:
\(V(S,v,T) = \max(\pm(S-K),0)\).

The closed-form pricer `heston_cf` uses the characteristic function with the "little Heston trap"
formulation (Albrecher et al. 2007) for numerical stability. The PDE solvers use Yanenko splitting
on \((S,v)\) with the same infrastructure as the 2D product-option solvers.

## Numerics

- **Grids:** Tavella-Randall nonuniform grids with parameter `alpha` (returns `n_grid+1` nodes for consistency with the PDE code).
- **Splitting:** Yanenko (ADI-like), implicit in S then implicit in X.
- **Linear solves:** Thomas algorithm for tridiagonal systems.
- **American:** Penalty method with `lambda` (typically 1e4-1e6) and relative tolerance, capped at 200 iterations per time step with a warning on non-convergence.
- **Interpolation:** Linear (1D) or bilinear (2D) at the spot price on the computed grid.

## Input validation

All solvers validate inputs at entry:

- Grid parameters: `n_s >= 2` (1D) or `>= 3` (2D), `n_t >= 1`, bounds increasing
- Sigma matrix dimensions must match grid: `(n_s+1) x (n_t+1)` for 1D, `(n_s+1) x n_t` for 2D local-vol
- Correlation: `rho` in `[-1, 1]` (2D solvers)
- Volatilities: `sigma >= 0` (constant-vol solvers)
- American parameters: `lambda > 0`, `tolerance > 0`

## Testing

115 tests across 8 test files verify correctness:

| Test file | Coverage |
|-----------|----------|
| `test_1d_european.R` | Black-Scholes convergence (call/put), put-call parity, non-zero s_min, input validation |
| `test_1d_american.R` | American >= European (put), non-dividend call = European, input validation |
| `test_closed_form.R` | Put-call parity, tau=0 / sigma=0 edge cases, correlation extremes, input rejection |
| `test_tavella_randall.R` | Endpoint exactness, monotonicity, node concentration, input rejection |
| `test_2d_local_equals_constant.R` | Local-vol = constant-vol (European/American, call/put), PDE vs closed-form with high r_f |
| `test_batch_solvers.R` | Batch vs single solver (European/American), mixed taus/types, input validation |
| `test_heston.R` | CF put-call parity, PDE vs CF convergence (call/put), American >= European, input validation |
| `test_mc_heston_4d.R` | MC positivity, reproducibility, put-call parity, MC vs PDE (xi→0), SE scaling, American >= European |

Run tests with:

```r
devtools::test()
```

## License

GPL (>= 2). See `DESCRIPTION`.
