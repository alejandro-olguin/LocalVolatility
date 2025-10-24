test_that("2D local-vol equals constant-vol when local matrices are constant", {
  # Parameters
  s_0   <- 100
  x_0   <- 20
  k     <- 2000
  tau   <- 1.0
  r_d   <- 0.05
  r_f   <- 0.02
  q     <- 0.01
  rho   <- 0.3
  type  <- "call"

  # Vols
  sigma_s <- 0.2
  sigma_x <- 0.15

  # Grids
  s_min <- 10;  s_max <- 300
  x_min <-  1;  x_max <-  60
  n_s <- 40; n_x <- 40
  n_t <- 50
  alpha <- 3.0

  # Build local-vol matrices filled with the constant vol values
  # Dimensions follow the C++ access pattern sigma_s(i, n) and sigma_x(j, n),
  # where i in 0..n_s and j in 0..n_x are spatial indices, and n in 0..n_t-1 is time index.
  SigmaS <- matrix(sigma_s, nrow = n_s + 1, ncol = n_t)
  SigmaX <- matrix(sigma_x, nrow = n_x + 1, ncol = n_t)

  # European option: constant-vol vs local-vol
  price_cv_eur <- european_option_2d(
    s_0, x_0, k, tau, r_d, r_f, q,
    sigma_s, sigma_x, rho, type,
    s_min, s_max, x_min, x_max,
    n_s, n_x, n_t, alpha
  )

  price_lv_eur <- european_option_lv_2d(
    s_0, x_0, k, tau, r_d, r_f, q,
    SigmaS, SigmaX, rho, type,
    s_min, s_max, x_min, x_max,
    n_s, n_x, n_t, alpha
  )

  expect_true(is.finite(price_cv_eur) && is.finite(price_lv_eur))
  expect_lt(abs(price_cv_eur - price_lv_eur), 1e-6)

  # American option: constant-vol vs local-vol (penalty method)
  lambda <- 5
  tol    <- 1e-8

  price_cv_am <- american_option_2d(
    s_0, x_0, k, tau, r_d, r_f, q,
    sigma_s, sigma_x, rho, type,
    s_min, s_max, x_min, x_max,
    n_s, n_x, n_t, alpha,
    lambda, tol
  )

  price_lv_am <- american_option_lv_2d(
    s_0, x_0, k, tau, r_d, r_f, q,
    SigmaS, SigmaX, rho, type,
    s_min, s_max, x_min, x_max,
    n_s, n_x, n_t, alpha,
    lambda, tol
  )

  expect_true(is.finite(price_cv_am) && is.finite(price_lv_am))
  expect_lt(abs(price_cv_am - price_lv_am), 1e-6)
})
