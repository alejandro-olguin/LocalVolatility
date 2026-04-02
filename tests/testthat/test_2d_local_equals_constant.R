test_that("2D local-vol equals constant-vol (European call)", {
  s_0 <- 100; x_0 <- 20; k <- 2000; tau <- 1.0
  r_d <- 0.05; r_f <- 0.02; q <- 0.01; rho <- 0.3
  sigma_s <- 0.2; sigma_x <- 0.15
  s_min <- 10; s_max <- 300; x_min <- 1; x_max <- 60
  n_s <- 40; n_x <- 40; n_t <- 50; alpha <- 3.0

  SigmaS <- matrix(sigma_s, nrow = n_s + 1, ncol = n_t)
  SigmaX <- matrix(sigma_x, nrow = n_x + 1, ncol = n_t)

  price_cv <- european_option_2d(
    s_0, x_0, k, tau, r_d, r_f, q, sigma_s, sigma_x, rho, "call",
    s_min, s_max, x_min, x_max, n_s, n_x, n_t, alpha
  )
  price_lv <- european_option_lv_2d(
    s_0, x_0, k, tau, r_d, r_f, q, SigmaS, SigmaX, rho, "call",
    s_min, s_max, x_min, x_max, n_s, n_x, n_t, alpha
  )

  expect_true(is.finite(price_cv) && is.finite(price_lv))
  expect_lt(abs(price_cv - price_lv), 1e-6)
})

test_that("2D local-vol equals constant-vol (European put)", {
  s_0 <- 100; x_0 <- 20; k <- 2000; tau <- 1.0
  r_d <- 0.05; r_f <- 0.02; q <- 0.01; rho <- 0.3
  sigma_s <- 0.2; sigma_x <- 0.15
  s_min <- 10; s_max <- 300; x_min <- 1; x_max <- 60
  n_s <- 40; n_x <- 40; n_t <- 50; alpha <- 3.0

  SigmaS <- matrix(sigma_s, nrow = n_s + 1, ncol = n_t)
  SigmaX <- matrix(sigma_x, nrow = n_x + 1, ncol = n_t)

  price_cv <- european_option_2d(
    s_0, x_0, k, tau, r_d, r_f, q, sigma_s, sigma_x, rho, "put",
    s_min, s_max, x_min, x_max, n_s, n_x, n_t, alpha
  )
  price_lv <- european_option_lv_2d(
    s_0, x_0, k, tau, r_d, r_f, q, SigmaS, SigmaX, rho, "put",
    s_min, s_max, x_min, x_max, n_s, n_x, n_t, alpha
  )

  expect_true(is.finite(price_cv) && is.finite(price_lv))
  expect_lt(abs(price_cv - price_lv), 1e-6)
})

test_that("2D local-vol equals constant-vol (American call)", {
  s_0 <- 100; x_0 <- 20; k <- 2000; tau <- 1.0
  r_d <- 0.05; r_f <- 0.02; q <- 0.01; rho <- 0.3
  sigma_s <- 0.2; sigma_x <- 0.15
  s_min <- 10; s_max <- 300; x_min <- 1; x_max <- 60
  n_s <- 40; n_x <- 40; n_t <- 50; alpha <- 3.0
  lambda <- 5; tol <- 1e-8

  SigmaS <- matrix(sigma_s, nrow = n_s + 1, ncol = n_t)
  SigmaX <- matrix(sigma_x, nrow = n_x + 1, ncol = n_t)

  price_cv <- american_option_2d(
    s_0, x_0, k, tau, r_d, r_f, q, sigma_s, sigma_x, rho, "call",
    s_min, s_max, x_min, x_max, n_s, n_x, n_t, alpha, lambda, tol
  )
  price_lv <- american_option_lv_2d(
    s_0, x_0, k, tau, r_d, r_f, q, SigmaS, SigmaX, rho, "call",
    s_min, s_max, x_min, x_max, n_s, n_x, n_t, alpha, lambda, tol
  )

  expect_true(is.finite(price_cv) && is.finite(price_lv))
  expect_lt(abs(price_cv - price_lv), 1e-6)
})

test_that("2D American prices are finite and positive", {
  s_0 <- 100; x_0 <- 20; k <- 2000; tau <- 1.0
  r_d <- 0.05; r_f <- 0.02; q <- 0.01; rho <- 0.3
  sigma_s <- 0.2; sigma_x <- 0.15
  s_min <- 10; s_max <- 300; x_min <- 1; x_max <- 60
  n_s <- 40; n_x <- 40; n_t <- 50; alpha <- 3.0
  lambda <- 5; tol <- 1e-8

  am_call <- american_option_2d(
    s_0, x_0, k, tau, r_d, r_f, q, sigma_s, sigma_x, rho, "call",
    s_min, s_max, x_min, x_max, n_s, n_x, n_t, alpha, lambda, tol
  )
  am_put <- american_option_2d(
    s_0, x_0, k, tau, r_d, r_f, q, sigma_s, sigma_x, rho, "put",
    s_min, s_max, x_min, x_max, n_s, n_x, n_t, alpha, lambda, tol
  )

  expect_true(is.finite(am_call) && am_call > 0)
  expect_true(is.finite(am_put) && am_put > 0)
})

test_that("2D European PDE matches closed-form with high r_f", {
  s_0 <- 100; x_0 <- 20; k <- 2000; tau <- 1.0
  r_d <- 0.05; r_f <- 0.10; q <- 0.01; rho <- 0.3
  sigma_s <- 0.2; sigma_x <- 0.15
  s_min <- 10; s_max <- 300; x_min <- 1; x_max <- 60
  n_s <- 60; n_x <- 60; n_t <- 80; alpha <- 3.0

  pde <- european_option_2d(
    s_0, x_0, k, tau, r_d, r_f, q, sigma_s, sigma_x, rho, "call",
    s_min, s_max, x_min, x_max, n_s, n_x, n_t, alpha
  )
  cf <- european_option_cf_2d(s_0, x_0, k, tau, r_d, r_f, q,
                               sigma_s, sigma_x, rho, 1, "call")

  expect_true(is.finite(pde) && is.finite(cf))
  # PDE should be reasonably close to closed-form (within 1% of price)
  expect_lt(abs(pde - cf) / cf, 0.01)
})
