# --- European batch solver tests ---

test_that("Batch European matches single solver (uniform rates, calls)", {
  s_0 <- 100; n_s <- 200; n_t <- 200; s_min <- 1; s_max <- 400
  Sigma <- matrix(0.2, nrow = n_s + 1, ncol = n_t + 1)

  k <- c(90, 100, 110)
  tau <- rep(1.0, 3)
  r_d <- rep(0.05, 3); r_f <- rep(0.02, 3)

  single <- vapply(seq_along(k), function(i) {
    european_option_lv(s_0, k[i], tau[i], r_d[i], r_f[i], Sigma, "call",
                       s_min, s_max, n_s, n_t)
  }, numeric(1))

  batch <- batch_price_european_lv(rep(s_0, 3), k, tau, r_d, r_f, Sigma,
                                    rep("call", 3), s_min, s_max, n_s, n_t)

  expect_equal(length(batch), 3L)
  for (i in seq_along(k)) {
    expect_lt(abs(batch[i] - single[i]) / single[i], 1e-4)
  }
})

test_that("Batch European matches single solver (puts)", {
  s_0 <- 100; n_s <- 200; n_t <- 200; s_min <- 1; s_max <- 400
  Sigma <- matrix(0.25, nrow = n_s + 1, ncol = n_t + 1)

  k <- c(90, 100, 110)
  tau <- rep(1.0, 3)
  r_d <- rep(0.05, 3); r_f <- rep(0.02, 3)

  single <- vapply(seq_along(k), function(i) {
    european_option_lv(s_0, k[i], tau[i], r_d[i], r_f[i], Sigma, "put",
                       s_min, s_max, n_s, n_t)
  }, numeric(1))

  batch <- batch_price_european_lv(rep(s_0, 3), k, tau, r_d, r_f, Sigma,
                                    rep("put", 3), s_min, s_max, n_s, n_t)

  for (i in seq_along(k)) {
    expect_lt(abs(batch[i] - single[i]) / max(single[i], 1e-10), 1e-4)
  }
})

test_that("Batch European handles mixed taus", {
  s_0 <- 100; n_s <- 200; n_t <- 200; s_min <- 1; s_max <- 400
  Sigma <- matrix(0.2, nrow = n_s + 1, ncol = n_t + 1)

  k <- c(100, 100, 100)
  tau <- c(0.25, 0.5, 1.0)
  r_d <- rep(0.05, 3); r_f <- rep(0.02, 3)

  single <- vapply(seq_along(k), function(i) {
    european_option_lv(s_0, k[i], tau[i], r_d[i], r_f[i], Sigma, "call",
                       s_min, s_max, n_s, n_t)
  }, numeric(1))

  batch <- batch_price_european_lv(rep(s_0, 3), k, tau, r_d, r_f, Sigma,
                                    rep("call", 3), s_min, s_max, n_s, n_t)

  # Shorter taus have slightly different dt alignment, allow 0.1% relative error
  for (i in seq_along(k)) {
    expect_lt(abs(batch[i] - single[i]) / single[i], 1e-3)
  }
})

test_that("Batch European rejects mismatched input lengths", {
  n_s <- 50; n_t <- 50; s_min <- 1; s_max <- 400
  Sigma <- matrix(0.2, nrow = n_s + 1, ncol = n_t + 1)

  expect_error(
    batch_price_european_lv(c(100, 100), c(100), c(1), c(0.05), c(0.02),
                             Sigma, c("call"), s_min, s_max, n_s, n_t),
    "same length"
  )
})

test_that("Batch European rejects wrong sigma dimensions", {
  n_s <- 50; n_t <- 50; s_min <- 1; s_max <- 400
  Sigma_bad <- matrix(0.2, nrow = n_s, ncol = n_t + 1) # wrong rows

  expect_error(
    batch_price_european_lv(100, 100, 1, 0.05, 0.02, Sigma_bad, "call",
                             s_min, s_max, n_s, n_t),
    "dimensions"
  )
})


# --- American batch solver tests ---

test_that("Batch American matches single solver (puts)", {
  s_0 <- 100; n_s <- 200; n_t <- 200; s_min <- 1; s_max <- 400
  Sigma <- matrix(0.2, nrow = n_s + 1, ncol = n_t + 1)
  lambda <- 1e4; tol <- 1e-8

  k <- c(90, 100, 110)
  tau <- rep(1.0, 3)
  r_d <- rep(0.05, 3); r_f <- rep(0.02, 3)

  single <- vapply(seq_along(k), function(i) {
    american_option_lv(s_0, k[i], tau[i], r_d[i], r_f[i], Sigma, "put",
                       s_min, s_max, n_s, n_t, lambda, tol)
  }, numeric(1))

  batch <- batch_price_american_lv(rep(s_0, 3), k, tau, r_d, r_f, Sigma,
                                    rep("put", 3), s_min, s_max, n_s, n_t,
                                    lambda, tol)

  expect_equal(length(batch), 3L)
  for (i in seq_along(k)) {
    expect_equal(batch[i], single[i], tolerance = 1e-10)
  }
})

test_that("Batch American matches single solver (calls)", {
  s_0 <- 100; n_s <- 200; n_t <- 200; s_min <- 1; s_max <- 400
  Sigma <- matrix(0.2, nrow = n_s + 1, ncol = n_t + 1)
  lambda <- 1e4; tol <- 1e-8

  k <- c(90, 100, 110)
  tau <- rep(1.0, 3)
  r_d <- rep(0.05, 3); r_f <- rep(0.02, 3)

  single <- vapply(seq_along(k), function(i) {
    american_option_lv(s_0, k[i], tau[i], r_d[i], r_f[i], Sigma, "call",
                       s_min, s_max, n_s, n_t, lambda, tol)
  }, numeric(1))

  batch <- batch_price_american_lv(rep(s_0, 3), k, tau, r_d, r_f, Sigma,
                                    rep("call", 3), s_min, s_max, n_s, n_t,
                                    lambda, tol)

  for (i in seq_along(k)) {
    expect_equal(batch[i], single[i], tolerance = 1e-10)
  }
})

test_that("Batch American handles mixed types", {
  s_0 <- 100; n_s <- 200; n_t <- 200; s_min <- 1; s_max <- 400
  Sigma <- matrix(0.2, nrow = n_s + 1, ncol = n_t + 1)
  lambda <- 1e4; tol <- 1e-8

  k <- c(100, 100)
  tau <- c(1.0, 1.0)
  r_d <- rep(0.05, 2); r_f <- rep(0.02, 2)
  types <- c("call", "put")

  single <- vapply(seq_along(k), function(i) {
    american_option_lv(s_0, k[i], tau[i], r_d[i], r_f[i], Sigma, types[i],
                       s_min, s_max, n_s, n_t, lambda, tol)
  }, numeric(1))

  batch <- batch_price_american_lv(rep(s_0, 2), k, tau, r_d, r_f, Sigma, types,
                                    s_min, s_max, n_s, n_t, lambda, tol)

  for (i in seq_along(k)) {
    expect_equal(batch[i], single[i], tolerance = 1e-10)
  }
})

test_that("Batch American prices are positive and finite", {
  s_0 <- 100; n_s <- 100; n_t <- 100; s_min <- 1; s_max <- 300
  Sigma <- matrix(0.3, nrow = n_s + 1, ncol = n_t + 1)
  lambda <- 1e4; tol <- 1e-8

  k <- seq(80, 120, by = 5)
  n_opt <- length(k)
  tau <- rep(0.5, n_opt)
  r_d <- rep(0.05, n_opt); r_f <- rep(0.02, n_opt)

  batch <- batch_price_american_lv(rep(s_0, n_opt), k, tau, r_d, r_f, Sigma,
                                    rep("put", n_opt), s_min, s_max, n_s, n_t,
                                    lambda, tol)

  expect_true(all(is.finite(batch)))
  expect_true(all(batch > 0))
})

test_that("Batch American rejects invalid lambda", {
  n_s <- 50; n_t <- 50; s_min <- 1; s_max <- 400
  Sigma <- matrix(0.2, nrow = n_s + 1, ncol = n_t + 1)

  expect_error(
    batch_price_american_lv(100, 100, 1, 0.05, 0.02, Sigma, "put",
                             s_min, s_max, n_s, n_t, -1, 1e-8),
    "lambda"
  )
})

test_that("Batch American rejects mismatched lengths", {
  n_s <- 50; n_t <- 50; s_min <- 1; s_max <- 400
  Sigma <- matrix(0.2, nrow = n_s + 1, ncol = n_t + 1)

  expect_error(
    batch_price_american_lv(c(100, 100), c(100), c(1), c(0.05), c(0.02),
                             Sigma, c("put"), s_min, s_max, n_s, n_t, 1e4, 1e-8),
    "same length"
  )
})
