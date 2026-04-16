# --------------------------------------------------------------------------- #
# Tests for Heston stochastic volatility solvers
# --------------------------------------------------------------------------- #

# Common parameters (Heston 1993 benchmark)
# Note: these violate Feller condition (2*kappa*theta = 0.16 < xi^2 = 0.25)
s_0   <- 100
k     <- 100
tau   <- 1.0
r_d   <- 0.05
q_div <- 0.0
kappa <- 2.0
theta <- 0.04
xi    <- 0.5
rho   <- -0.7
v_0   <- 0.04

# PDE grid parameters
s_min <- 20;  s_max <- 300
v_min <- 0.001; v_max <- 1.0
n_s <- 80; n_v <- 40; n_t <- 100
alpha <- 3.0

# =========================================================================== #
# Closed-form (characteristic function) tests
# =========================================================================== #

test_that("Heston CF: put-call parity holds", {
  call_cf <- suppressWarnings(
    heston_cf(s_0, k, tau, r_d, q_div, kappa, theta, xi, rho, v_0, "call"))
  put_cf <- suppressWarnings(
    heston_cf(s_0, k, tau, r_d, q_div, kappa, theta, xi, rho, v_0, "put"))
  parity <- call_cf - put_cf - (s_0 * exp(-q_div * tau) - k * exp(-r_d * tau))
  expect_lt(abs(parity), 1e-6)
})

test_that("Heston CF: call price is positive and finite", {
  price <- suppressWarnings(
    heston_cf(s_0, k, tau, r_d, q_div, kappa, theta, xi, rho, v_0, "call"))
  expect_true(is.finite(price))
  expect_gt(price, 0)
})

test_that("Heston CF: put price is positive and finite", {
  price <- suppressWarnings(
    heston_cf(s_0, k, tau, r_d, q_div, kappa, theta, xi, rho, v_0, "put"))
  expect_true(is.finite(price))
  expect_gt(price, 0)
})

test_that("Heston CF: deep ITM call is close to intrinsic", {
  price <- suppressWarnings(
    heston_cf(200, 100, tau, r_d, q_div, kappa, theta, xi, rho, v_0, "call"))
  intrinsic <- 200 * exp(-q_div * tau) - 100 * exp(-r_d * tau)
  expect_gt(price, intrinsic * 0.95)
})

test_that("Heston CF: deep OTM call is near zero", {
  price <- suppressWarnings(
    heston_cf(50, 150, tau, r_d, q_div, kappa, theta, xi, rho, v_0, "call"))
  expect_lt(price, 1.0)
})

test_that("Heston CF: rejects invalid inputs", {
  expect_error(heston_cf(s_0, k, tau, r_d, q_div, kappa, theta, xi, 1.5, v_0, "call"))
  expect_error(heston_cf(s_0, k, tau, r_d, q_div, -1, theta, xi, rho, v_0, "call"))
  expect_error(heston_cf(s_0, k, tau, r_d, q_div, kappa, theta, -0.1, rho, v_0, "call"))
  expect_error(heston_cf(s_0, k, tau, r_d, q_div, kappa, theta, xi, rho, -0.01, "call"))
  expect_error(heston_cf(s_0, k, tau, r_d, q_div, kappa, theta, xi, rho, v_0, "other"))
})

test_that("Heston CF: warns when Feller condition violated", {
  expect_warning(
    heston_cf(s_0, k, tau, r_d, q_div, kappa, theta, xi, rho, v_0, "call"),
    "Feller")
})

# =========================================================================== #
# European PDE tests
# =========================================================================== #

test_that("Heston PDE: European call matches closed-form within 2%", {
  price_pde <- suppressWarnings(european_option_heston(
    s_0, k, tau, r_d, q_div, kappa, theta, xi, rho, v_0, "call",
    s_min, s_max, v_min, v_max, n_s, n_v, n_t, alpha))
  price_cf <- suppressWarnings(
    heston_cf(s_0, k, tau, r_d, q_div, kappa, theta, xi, rho, v_0, "call"))
  expect_lt(abs(price_pde - price_cf) / price_cf, 0.02)
})

test_that("Heston PDE: European put matches closed-form within 2%", {
  price_pde <- suppressWarnings(european_option_heston(
    s_0, k, tau, r_d, q_div, kappa, theta, xi, rho, v_0, "put",
    s_min, s_max, v_min, v_max, n_s, n_v, n_t, alpha))
  price_cf <- suppressWarnings(
    heston_cf(s_0, k, tau, r_d, q_div, kappa, theta, xi, rho, v_0, "put"))
  expect_lt(abs(price_pde - price_cf) / price_cf, 0.02)
})

test_that("Heston PDE: European put-call parity (loose)", {
  call_pde <- suppressWarnings(european_option_heston(
    s_0, k, tau, r_d, q_div, kappa, theta, xi, rho, v_0, "call",
    s_min, s_max, v_min, v_max, n_s, n_v, n_t, alpha))
  put_pde <- suppressWarnings(european_option_heston(
    s_0, k, tau, r_d, q_div, kappa, theta, xi, rho, v_0, "put",
    s_min, s_max, v_min, v_max, n_s, n_v, n_t, alpha))
  parity <- call_pde - put_pde - (s_0 * exp(-q_div * tau) - k * exp(-r_d * tau))
  expect_lt(abs(parity), 1.0)
})

test_that("Heston PDE: prices are positive and finite", {
  price <- suppressWarnings(european_option_heston(
    s_0, k, tau, r_d, q_div, kappa, theta, xi, rho, v_0, "call",
    s_min, s_max, v_min, v_max, n_s, n_v, n_t, alpha))
  expect_true(is.finite(price))
  expect_gt(price, 0)
})

test_that("Heston PDE: rejects invalid inputs", {
  expect_error(european_option_heston(
    s_0, k, tau, r_d, q_div, kappa, theta, xi, 1.5, v_0, "call",
    s_min, s_max, v_min, v_max, n_s, n_v, n_t, alpha))
  expect_error(european_option_heston(
    s_0, k, tau, r_d, q_div, kappa, theta, xi, rho, v_0, "call",
    s_min, s_max, -0.01, v_max, n_s, n_v, n_t, alpha))
  expect_error(european_option_heston(
    s_0, k, tau, r_d, q_div, kappa, theta, xi, rho, v_0, "call",
    s_min, s_max, v_min, v_max, 2, n_v, n_t, alpha))
})

# =========================================================================== #
# American PDE tests
# =========================================================================== #

test_that("Heston American put >= European put", {
  eu <- suppressWarnings(european_option_heston(
    s_0, k, tau, r_d, q_div, kappa, theta, xi, rho, v_0, "put",
    s_min, s_max, v_min, v_max, n_s, n_v, n_t, alpha))
  am <- suppressWarnings(american_option_heston(
    s_0, k, tau, r_d, q_div, kappa, theta, xi, rho, v_0, "put",
    s_min, s_max, v_min, v_max, n_s, n_v, n_t, alpha, 1e4, 1e-8))
  expect_gte(am, eu - 0.01)
})

test_that("Heston American: price is positive and finite", {
  price <- suppressWarnings(american_option_heston(
    s_0, k, tau, r_d, q_div, kappa, theta, xi, rho, v_0, "put",
    s_min, s_max, v_min, v_max, n_s, n_v, n_t, alpha, 1e4, 1e-8))
  expect_true(is.finite(price))
  expect_gt(price, 0)
})

test_that("Heston American: rejects invalid inputs", {
  expect_error(american_option_heston(
    s_0, k, tau, r_d, q_div, kappa, theta, xi, rho, v_0, "put",
    s_min, s_max, v_min, v_max, n_s, n_v, n_t, alpha, -1, 1e-8))
  expect_error(american_option_heston(
    s_0, k, tau, r_d, q_div, kappa, theta, xi, rho, v_0, "put",
    s_min, s_max, v_min, v_max, n_s, n_v, n_t, alpha, 1e4, -1))
})
