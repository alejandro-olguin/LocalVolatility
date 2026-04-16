test_that("1D American put >= European put", {
  s_0 <- 100; k <- 100; tau <- 1.0; r_d <- 0.05; r_f <- 0.02
  sigma <- 0.25
  s_min <- 0; s_max <- 300
  n_s <- 150; n_t <- 150

  SigmaMat <- matrix(sigma, nrow = n_s + 1, ncol = n_t + 1)

  eur <- european_option_lv(s_0, k, tau, r_d, r_f, SigmaMat, "put",
                             s_min, s_max, n_s, n_t)
  am  <- american_option_lv(s_0, k, tau, r_d, r_f, SigmaMat, "put",
                             s_min, s_max, n_s, n_t,
                             lambda = 1e4, tolerance = 1e-8)

  expect_true(is.finite(eur) && is.finite(am))
  expect_gte(am, eur - 1e-8)
})

test_that("1D American call on non-dividend stock equals European", {
  s_0 <- 100; k <- 100; tau <- 1.0; r_d <- 0.05; r_f <- 0.0
  sigma <- 0.2
  s_min <- 0; s_max <- 300
  n_s <- 150; n_t <- 150

  SigmaMat <- matrix(sigma, nrow = n_s + 1, ncol = n_t + 1)

  eur <- european_option_lv(s_0, k, tau, r_d, r_f, SigmaMat, "call",
                             s_min, s_max, n_s, n_t)
  am  <- american_option_lv(s_0, k, tau, r_d, r_f, SigmaMat, "call",
                             s_min, s_max, n_s, n_t,
                             lambda = 1e4, tolerance = 1e-8)

  expect_true(is.finite(eur) && is.finite(am))
  # With r_f=0, American call = European call
  expect_lt(abs(am - eur) / eur, 0.01)
})

test_that("1D American put price is positive and finite", {
  s_0 <- 100; k <- 110; tau <- 0.5; r_d <- 0.03; r_f <- 0.01
  sigma <- 0.3
  s_min <- 0; s_max <- 250
  n_s <- 100; n_t <- 100

  SigmaMat <- matrix(sigma, nrow = n_s + 1, ncol = n_t + 1)

  price <- american_option_lv(s_0, k, tau, r_d, r_f, SigmaMat, "put",
                               s_min, s_max, n_s, n_t,
                               lambda = 1e4, tolerance = 1e-8)

  expect_true(is.finite(price))
  expect_gt(price, 0)
})

test_that("1D American rejects invalid inputs", {
  n_s <- 50; n_t <- 50
  bad_sigma <- matrix(0.2, nrow = 50, ncol = 50)  # wrong dims

  expect_error(
    american_option_lv(100, 100, 1.0, 0.05, 0.0, bad_sigma, "call",
                        0, 200, n_s, n_t, 1e4, 1e-8),
    "sigma"
  )
})
