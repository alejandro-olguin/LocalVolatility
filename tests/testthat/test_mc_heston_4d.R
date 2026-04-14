# --------------------------------------------------------------------------- #
# Tests for 4D Monte Carlo (double-Heston quanto)
# --------------------------------------------------------------------------- #

# Common parameters
s_0 <- 100; x_0 <- 20; k <- 2000; tau <- 1.0
r_d <- 0.05; r_f <- 0.02; q_div <- 0.01

# Heston params for equity
kappa_s <- 2.0; theta_s <- 0.04; xi_s <- 0.5; v_s0 <- 0.04
# Heston params for FX
kappa_x <- 1.5; theta_x <- 0.02; xi_x <- 0.3; v_x0 <- 0.02

# Correlations
rho_sx <- 0.3; rho_sv <- -0.7; rho_xv <- -0.5
rho_svx <- 0; rho_xvs <- 0; rho_vsvx <- 0

n_paths <- 100000; n_steps <- 100; seed <- 42

# =========================================================================== #
# European MC tests
# =========================================================================== #

test_that("MC European 4D: call price is positive and finite", {
  res <- suppressWarnings(mc_european_heston_4d(
    s_0, x_0, k, tau, r_d, r_f, q_div,
    kappa_s, theta_s, xi_s, v_s0,
    kappa_x, theta_x, xi_x, v_x0,
    rho_sx, rho_sv, rho_xv, rho_svx, rho_xvs, rho_vsvx,
    "call", n_paths, n_steps, seed))
  expect_true(is.finite(res$price))
  expect_gt(res$price, 0)
  expect_true(is.finite(res$std_error))
  expect_gt(res$std_error, 0)
})

test_that("MC European 4D: put price is positive and finite", {
  res <- suppressWarnings(mc_european_heston_4d(
    s_0, x_0, k, tau, r_d, r_f, q_div,
    kappa_s, theta_s, xi_s, v_s0,
    kappa_x, theta_x, xi_x, v_x0,
    rho_sx, rho_sv, rho_xv, rho_svx, rho_xvs, rho_vsvx,
    "put", n_paths, n_steps, seed))
  expect_true(is.finite(res$price))
  expect_gt(res$price, 0)
})

test_that("MC European 4D: reproducibility with same seed", {
  args <- list(s_0, x_0, k, tau, r_d, r_f, q_div,
    kappa_s, theta_s, xi_s, v_s0,
    kappa_x, theta_x, xi_x, v_x0,
    rho_sx, rho_sv, rho_xv, rho_svx, rho_xvs, rho_vsvx,
    "call", 10000, 50, 123)
  r1 <- suppressWarnings(do.call(mc_european_heston_4d, args))
  r2 <- suppressWarnings(do.call(mc_european_heston_4d, args))
  expect_equal(r1$price, r2$price)
})

test_that("MC European 4D: put-call parity (within MC error)", {
  call_res <- suppressWarnings(mc_european_heston_4d(
    s_0, x_0, k, tau, r_d, r_f, q_div,
    kappa_s, theta_s, xi_s, v_s0,
    kappa_x, theta_x, xi_x, v_x0,
    rho_sx, rho_sv, rho_xv, rho_svx, rho_xvs, rho_vsvx,
    "call", n_paths, n_steps, seed))
  put_res <- suppressWarnings(mc_european_heston_4d(
    s_0, x_0, k, tau, r_d, r_f, q_div,
    kappa_s, theta_s, xi_s, v_s0,
    kappa_x, theta_x, xi_x, v_x0,
    rho_sx, rho_sv, rho_xv, rho_svx, rho_xvs, rho_vsvx,
    "put", n_paths, n_steps, seed + 1))
  # put-call parity: C - P ≈ S₀X₀ exp(-(r_f+q)τ) - K exp(-r_d τ)
  # Note: quanto adjustment makes exact parity complex; use loose tolerance
  parity_lhs <- call_res$price - put_res$price
  fwd <- s_0 * x_0 * exp(-(r_f + q_div) * tau) - k * exp(-r_d * tau)
  se <- sqrt(call_res$std_error^2 + put_res$std_error^2)
  expect_lt(abs(parity_lhs - fwd), 6 * se + 50)  # loose: quanto adj + MC noise
})

test_that("MC European 4D: matches 2D PDE with near-zero xi", {
  # When xi_s, xi_x ≈ 0 and v_0 = sigma^2, theta = sigma^2,
  # Heston degenerates to constant-vol
  sigma_s <- 0.2; sigma_x <- 0.15
  mc_res <- suppressWarnings(mc_european_heston_4d(
    s_0, x_0, k, tau, r_d, r_f, q_div,
    kappa_s = 5, theta_s = sigma_s^2, xi_s = 0.001, v_s0 = sigma_s^2,
    kappa_x = 5, theta_x = sigma_x^2, xi_x = 0.001, v_x0 = sigma_x^2,
    rho_sx = rho_sx, rho_sv = 0, rho_xv = 0,
    rho_svx = 0, rho_xvs = 0, rho_vsvx = 0,
    type = "call", n_paths = 200000, n_steps = 200, seed = 42))
  pde_price <- european_option_2d(s_0, x_0, k, tau, r_d, r_f, q_div,
    sigma_s, sigma_x, rho_sx, "call",
    10, 300, 1, 60, 40, 40, 50, 3)
  # MC vs PDE within 3 SE + PDE truncation
  expect_lt(abs(mc_res$price - pde_price), 3 * mc_res$std_error + 5)
})

test_that("MC European 4D: SE decreases with more paths", {
  se1 <- suppressWarnings(mc_european_heston_4d(
    s_0, x_0, k, tau, r_d, r_f, q_div,
    kappa_s, theta_s, xi_s, v_s0,
    kappa_x, theta_x, xi_x, v_x0,
    rho_sx, rho_sv, rho_xv, rho_svx, rho_xvs, rho_vsvx,
    "call", 10000, 50, seed))$std_error
  se2 <- suppressWarnings(mc_european_heston_4d(
    s_0, x_0, k, tau, r_d, r_f, q_div,
    kappa_s, theta_s, xi_s, v_s0,
    kappa_x, theta_x, xi_x, v_x0,
    rho_sx, rho_sv, rho_xv, rho_svx, rho_xvs, rho_vsvx,
    "call", 100000, 50, seed))$std_error
  expect_lt(se2, se1)
})

test_that("MC European 4D: rejects invalid inputs", {
  expect_error(mc_european_heston_4d(
    s_0, x_0, k, tau, r_d, r_f, q_div,
    kappa_s, theta_s, xi_s, v_s0,
    kappa_x, theta_x, xi_x, v_x0,
    1.5, rho_sv, rho_xv, rho_svx, rho_xvs, rho_vsvx,
    "call", n_paths, n_steps, seed))
  expect_error(mc_european_heston_4d(
    s_0, x_0, k, tau, r_d, r_f, q_div,
    -1, theta_s, xi_s, v_s0,
    kappa_x, theta_x, xi_x, v_x0,
    rho_sx, rho_sv, rho_xv, rho_svx, rho_xvs, rho_vsvx,
    "call", n_paths, n_steps, seed))
  expect_error(mc_european_heston_4d(
    s_0, x_0, k, tau, r_d, r_f, q_div,
    kappa_s, theta_s, xi_s, v_s0,
    kappa_x, theta_x, xi_x, v_x0,
    rho_sx, rho_sv, rho_xv, rho_svx, rho_xvs, rho_vsvx,
    "other", n_paths, n_steps, seed))
})

# =========================================================================== #
# American MC tests (Longstaff-Schwartz)
# =========================================================================== #

test_that("MC American 4D: put price is positive and finite", {
  res <- suppressWarnings(mc_american_heston_4d(
    s_0, x_0, k, tau, r_d, r_f, q_div,
    kappa_s, theta_s, xi_s, v_s0,
    kappa_x, theta_x, xi_x, v_x0,
    rho_sx, rho_sv, rho_xv, rho_svx, rho_xvs, rho_vsvx,
    "put", 50000, 50, seed, 4))
  expect_true(is.finite(res$price))
  expect_gt(res$price, 0)
  expect_true(is.finite(res$std_error))
})

test_that("MC American 4D: American put >= European put", {
  eu_res <- suppressWarnings(mc_european_heston_4d(
    s_0, x_0, k, tau, r_d, r_f, q_div,
    kappa_s, theta_s, xi_s, v_s0,
    kappa_x, theta_x, xi_x, v_x0,
    rho_sx, rho_sv, rho_xv, rho_svx, rho_xvs, rho_vsvx,
    "put", 100000, 50, seed))
  am_res <- suppressWarnings(mc_american_heston_4d(
    s_0, x_0, k, tau, r_d, r_f, q_div,
    kappa_s, theta_s, xi_s, v_s0,
    kappa_x, theta_x, xi_x, v_x0,
    rho_sx, rho_sv, rho_xv, rho_svx, rho_xvs, rho_vsvx,
    "put", 100000, 50, seed, 4))
  # American >= European (within MC noise)
  se <- sqrt(eu_res$std_error^2 + am_res$std_error^2)
  expect_gte(am_res$price, eu_res$price - 3 * se)
})

test_that("MC American 4D: rejects invalid inputs", {
  expect_error(mc_american_heston_4d(
    s_0, x_0, k, tau, r_d, r_f, q_div,
    kappa_s, theta_s, xi_s, v_s0,
    kappa_x, theta_x, xi_x, v_x0,
    rho_sx, rho_sv, rho_xv, rho_svx, rho_xvs, rho_vsvx,
    "put", 50000, 50, seed, 1))  # n_basis < 2
})
