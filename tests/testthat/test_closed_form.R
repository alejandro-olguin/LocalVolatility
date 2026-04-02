test_that("Closed-form put-call parity holds", {
  s_0 <- 100; x_0 <- 20; k <- 2000; tau <- 1.0
  r_d <- 0.05; r_f <- 0.02; q <- 0.01; rho <- 0.3
  sigma_s <- 0.2; sigma_x <- 0.15; n <- 1

  call_price <- european_option_cf_2d(s_0, x_0, k, tau, r_d, r_f, q,
                                       sigma_s, sigma_x, rho, n, "call")
  put_price  <- european_option_cf_2d(s_0, x_0, k, tau, r_d, r_f, q,
                                       sigma_s, sigma_x, rho, n, "put")

  adr_0 <- (s_0 * x_0) / n
  parity <- call_price - put_price - (adr_0 * exp(-q * tau) - k * exp(-r_d * tau))
  expect_lt(abs(parity), 1e-10)
})

test_that("Closed-form tau=0 returns intrinsic", {
  s_0 <- 100; x_0 <- 25; k <- 2000; tau <- 0.0
  r_d <- 0.05; r_f <- 0.02; q <- 0.01; rho <- 0.3
  sigma_s <- 0.2; sigma_x <- 0.15; n <- 1

  call_price <- european_option_cf_2d(s_0, x_0, k, tau, r_d, r_f, q,
                                       sigma_s, sigma_x, rho, n, "call")
  put_price  <- european_option_cf_2d(s_0, x_0, k, tau, r_d, r_f, q,
                                       sigma_s, sigma_x, rho, n, "put")

  adr_0 <- (s_0 * x_0) / n
  expect_equal(call_price, max(adr_0 - k, 0))
  expect_equal(put_price,  max(k - adr_0, 0))
})

test_that("Closed-form sigma=0 returns discounted intrinsic", {
  s_0 <- 100; x_0 <- 25; k <- 2000; tau <- 1.0
  r_d <- 0.05; r_f <- 0.02; q <- 0.01; rho <- 0.0
  n <- 1

  call_price <- european_option_cf_2d(s_0, x_0, k, tau, r_d, r_f, q,
                                       0.0, 0.0, rho, n, "call")

  adr_0 <- (s_0 * x_0) / n
  fwd <- adr_0 * exp((r_d - q) * tau)
  expected <- exp(-r_d * tau) * max(fwd - k, 0)
  expect_equal(call_price, expected)
})

test_that("Closed-form rejects invalid inputs", {
  expect_error(european_option_cf_2d(100, 20, 2000, 1, 0.05, 0.02, 0.01,
                                      0.2, 0.15, 1.5, 1, "call"))  # rho > 1
  expect_error(european_option_cf_2d(0, 20, 2000, 1, 0.05, 0.02, 0.01,
                                      0.2, 0.15, 0.3, 1, "call"))  # s_0 = 0
  expect_error(european_option_cf_2d(100, 20, 2000, 1, 0.05, 0.02, 0.01,
                                      0.2, 0.15, 0.3, 0, "call"))  # n = 0
})

test_that("Closed-form correlation extremes produce valid prices", {
  s_0 <- 100; x_0 <- 20; k <- 2000; tau <- 1.0
  r_d <- 0.05; r_f <- 0.02; q <- 0.01
  sigma_s <- 0.2; sigma_x <- 0.15; n <- 1

  for (rho in c(-1.0, 0.0, 1.0)) {
    call_price <- european_option_cf_2d(s_0, x_0, k, tau, r_d, r_f, q,
                                         sigma_s, sigma_x, rho, n, "call")
    expect_true(is.finite(call_price))
    expect_gte(call_price, 0)
  }
})
