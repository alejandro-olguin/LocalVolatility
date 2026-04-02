test_that("1D European call vs Black-Scholes closed-form", {
  s_0 <- 100; k <- 100; tau <- 1.0; r_d <- 0.05; q <- 0.02
  sigma <- 0.2
  s_min <- 0; s_max <- 300
  n_s <- 200; n_t <- 200

  # Constant vol matrix: (n_s+1) x (n_t+1)
  SigmaMat <- matrix(sigma, nrow = n_s + 1, ncol = n_t + 1)

  pde <- european_option_lv(s_0, k, tau, r_d, q, SigmaMat, "call",
                             s_min, s_max, n_s, n_t)

  # Black-Scholes analytical
  d1 <- (log(s_0 / k) + (r_d - q + 0.5 * sigma^2) * tau) / (sigma * sqrt(tau))
  d2 <- d1 - sigma * sqrt(tau)
  bs <- s_0 * exp(-q * tau) * pnorm(d1) - k * exp(-r_d * tau) * pnorm(d2)

  expect_true(is.finite(pde))
  expect_lt(abs(pde - bs) / bs, 0.005)  # within 0.5%
})

test_that("1D European put vs Black-Scholes closed-form", {
  s_0 <- 100; k <- 110; tau <- 0.5; r_d <- 0.03; q <- 0.01
  sigma <- 0.25
  s_min <- 0; s_max <- 300
  n_s <- 200; n_t <- 200

  SigmaMat <- matrix(sigma, nrow = n_s + 1, ncol = n_t + 1)

  pde <- european_option_lv(s_0, k, tau, r_d, q, SigmaMat, "put",
                             s_min, s_max, n_s, n_t)

  d1 <- (log(s_0 / k) + (r_d - q + 0.5 * sigma^2) * tau) / (sigma * sqrt(tau))
  d2 <- d1 - sigma * sqrt(tau)
  bs <- k * exp(-r_d * tau) * pnorm(-d2) - s_0 * exp(-q * tau) * pnorm(-d1)

  expect_true(is.finite(pde))
  expect_lt(abs(pde - bs) / bs, 0.005)
})

test_that("1D European put-call parity holds", {
  s_0 <- 100; k <- 100; tau <- 1.0; r_d <- 0.05; q <- 0.02
  sigma <- 0.2
  s_min <- 0; s_max <- 300
  n_s <- 200; n_t <- 200

  SigmaMat <- matrix(sigma, nrow = n_s + 1, ncol = n_t + 1)

  call_p <- european_option_lv(s_0, k, tau, r_d, q, SigmaMat, "call",
                                s_min, s_max, n_s, n_t)
  put_p  <- european_option_lv(s_0, k, tau, r_d, q, SigmaMat, "put",
                                s_min, s_max, n_s, n_t)

  parity <- call_p - put_p - (s_0 * exp(-q * tau) - k * exp(-r_d * tau))
  expect_lt(abs(parity), 0.05)  # PDE tolerance
})

test_that("1D European with non-zero s_min works correctly", {
  s_0 <- 100; k <- 100; tau <- 1.0; r_d <- 0.05; q <- 0.0
  sigma <- 0.2
  s_min <- 50; s_max <- 200
  n_s <- 150; n_t <- 150

  SigmaMat <- matrix(sigma, nrow = n_s + 1, ncol = n_t + 1)

  price <- european_option_lv(s_0, k, tau, r_d, q, SigmaMat, "call",
                               s_min, s_max, n_s, n_t)

  expect_true(is.finite(price))
  expect_gt(price, 0)
})

test_that("1D European rejects invalid sigma dimensions", {
  s_0 <- 100; k <- 100; tau <- 1.0; r_d <- 0.05; q <- 0.0
  n_s <- 50; n_t <- 50

  # Wrong dimensions: should be (n_s+1) x (n_t+1) = 51 x 51
  bad_sigma <- matrix(0.2, nrow = 50, ncol = 50)

  expect_error(
    european_option_lv(s_0, k, tau, r_d, q, bad_sigma, "call", 0, 200, n_s, n_t),
    "sigma"
  )
})
