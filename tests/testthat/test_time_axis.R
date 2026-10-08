# Calendar-time convention for sigma matrices, exact maturities in the batch
# solvers, and far-field boundary consistency. These tests guard against the
# three behaviours fixed in 3.3.0: relative-time / time-to-maturity indexing of
# sigma columns, maturity rounding to the shared time grid, and the extra r_f
# discount in the 2D call asymptote.

bs_call <- function(S, K, tau, r, q, sig) {
  d1 <- (log(S / K) + (r - q + 0.5 * sig^2) * tau) / (sig * sqrt(tau))
  d2 <- d1 - sig * sqrt(tau)
  S * exp(-q * tau) * pnorm(d1) - K * exp(-r * tau) * pnorm(d2)
}
iv_call <- function(p, S, K, tau, r, q)
  uniroot(function(s) bs_call(S, K, tau, r, q, s) - p, c(1e-4, 5), tol = 1e-10)$root

# Step local vol: 0.40 for calendar t < 0.25, 0.10 afterwards, on [0, 1].
S0 <- 100; K <- 100; r <- 0.03; t_max <- 1
n_s <- 300; n_t <- 300
tgrid   <- seq(0, t_max, length.out = n_t + 1)
sig_cal <- matrix(ifelse(tgrid < 0.25, 0.40, 0.10), n_s + 1, n_t + 1, byrow = TRUE)
iv_1y   <- sqrt(0.25 * 0.16 + 0.75 * 0.01)

test_that("batch European reads sigma in calendar time (step-vol test)", {
  p <- batch_price_european_lv(rep(S0, 2), c(K, K), c(0.25, 1), c(r, r), c(0, 0),
                               sig_cal, c("call", "call"), 1, 400, n_s, n_t)
  expect_lt(abs(iv_call(p[1], S0, K, 0.25, r, 0) - 0.40), 0.004)
  expect_lt(abs(iv_call(p[2], S0, K, 1.00, r, 0) - iv_1y), 0.004)
})

test_that("batch American reads sigma in calendar time (step-vol test)", {
  p <- batch_price_american_lv(rep(S0, 2), c(K, K), c(0.25, 1), c(r, r), c(0, 0),
                               sig_cal, c("call", "call"), 1, 400, n_s, n_t,
                               1e6, 1e-8)
  expect_lt(abs(iv_call(p[1], S0, K, 0.25, r, 0) - 0.40), 0.004)
  expect_lt(abs(iv_call(p[2], S0, K, 1.00, r, 0) - iv_1y), 0.004)
})

test_that("2D LV solvers read a shared sigma matrix in calendar time", {
  n <- 50; nt <- 150
  tg <- seq(0, t_max, length.out = nt)
  SS <- matrix(ifelse(tg < 0.25, 0.40, 0.10), n + 1, nt, byrow = TRUE)
  XX <- matrix(1e-4, n + 1, nt)
  pe <- european_option_lv_2d(S0, 1, K, 0.25, r, r, 0, SS, XX, 0, "call",
                              1, 400, 0.5, 2, n, n, nt, alpha = 0.5 * S0,
                              sigma_t_max = t_max)
  pa <- american_option_lv_2d(S0, 1, K, 0.25, r, r, 0, SS, XX, 0, "call",
                              1, 400, 0.5, 2, n, n, nt, alpha = 0.5 * S0,
                              lambda = 500, tolerance = 1e-6, sigma_t_max = t_max)
  expect_lt(abs(iv_call(pe, S0, K, 0.25, r, 0) - 0.40), 0.01)
  expect_lt(abs(iv_call(pa, S0, K, 0.25, r, 0) - 0.40), 0.01)
  # Default sigma_t_max = tau: a matrix describing exactly [0, tau].
  SSo <- matrix(ifelse(seq(0, 0.25, length.out = nt) < 0.25, 0.40, 0.10),
                n + 1, nt, byrow = TRUE)
  po <- european_option_lv_2d(S0, 1, K, 0.25, r, r, 0, SSo, XX, 0, "call",
                              1, 400, 0.5, 2, n, n, nt, alpha = 0.5 * S0)
  expect_lt(abs(iv_call(po, S0, K, 0.25, r, 0) - 0.40), 0.01)
})

test_that("batch European matches single solver with mixed taus and time-varying sigma", {
  taus <- c(0.25, 0.6, 1.0)
  batch <- batch_price_european_lv(rep(S0, 3), rep(K, 3), taus, rep(r, 3), rep(0, 3),
                                   sig_cal, rep("call", 3), 1, 400, n_s, n_t)
  for (i in seq_along(taus)) {
    tt     <- seq(0, taus[i], length.out = n_t + 1)
    sig_i  <- matrix(ifelse(tt < 0.25, 0.40, 0.10), n_s + 1, n_t + 1, byrow = TRUE)
    single <- european_option_lv(S0, K, taus[i], r, 0, sig_i, "call", 1, 400, n_s, n_t)
    # The shared grid interpolates the step between columns; the per-option
    # matrix places it on a node, so agreement is at discretisation level.
    expect_lt(abs(batch[i] - single) / single, 5e-3)
  }
})

test_that("batch European does not round short maturities to the shared grid", {
  x0 <- 1.26; sig <- matrix(0.09, 401, 351)
  taus <- c(7 / 365, 14 / 365, 1.5)
  p <- batch_price_european_lv(rep(x0, 3), rep(1.26, 3), taus, rep(0.05, 3), rep(0.045, 3),
                               sig, rep("call", 3), 0.3 * x0, 3 * x0, 400, 350)
  for (i in 1:3)
    expect_lt(abs(iv_call(p[i], x0, 1.26, taus[i], 0.05, 0.045) - 0.09), 0.004)
})

test_that("2D European call far-field boundary is consistent with the closed form at 5y", {
  s0 <- 30; x0 <- 1.25; K <- 38; tau <- 5
  rd <- 0.05; rf <- 0.045; q <- 0.03; sS <- 0.25; sX <- 0.08; rho <- -0.1
  cf  <- european_option_cf_2d(s0, x0, K, tau, rd, rf, q, sS, sX, rho, "call", 1)
  pde <- european_option_2d(s0, x0, K, tau, rd, rf, q, sS, sX, rho, "call",
                            0.01 * s0, 8 * s0, 0.01 * x0, 8 * x0, 120, 120, 150,
                            alpha = 0.5 * s0)
  expect_lt(abs(pde - cf) / cf, 0.01)
})

test_that("2D American price is insensitive to the penalty parameter once scaled", {
  s0 <- 30; x0 <- 1.25; K <- 38; tau <- 1
  SS <- matrix(0.25, 61, 60); XX <- matrix(0.08, 61, 60)
  p <- sapply(c(100, 1e4, 1e6), function(l)
    american_option_lv_2d(s0, x0, K, tau, 0.05, 0.045, 0.06, SS, XX, -0.1, "call",
                          0.01 * s0, 4 * s0, 0.01 * x0, 4 * x0, 60, 60, 60,
                          alpha = 0.5 * s0, lambda = l, tolerance = 1e-7))
  expect_lt(max(p) - min(p), 1e-3)
})
