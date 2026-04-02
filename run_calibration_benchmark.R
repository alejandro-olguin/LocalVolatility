#!/usr/bin/env Rscript
# Calibration benchmark: batch solvers vs serial solvers
# Uses saved calibration data from the fitting project

library(LocalVolatility)
library(akima)
library(lbfgsb3c)

# --- Load previous calibration data ---
results_file <- "/Users/alejandro/Documents/Projects/Fitting Local Volatility/results/calibration_results_20251110_154314.RData"
load(results_file)
cat("Loaded calibration data from:", results_file, "\n")
cat("Stock options:", nrow(stock_data_sparse), " | FX options:", nrow(fx_data_sparse), "\n")
cat("Grid: n_s =", n_s, "n_t =", n_t, "\n\n")

# --- Helper functions (from calibration_functions.R) ---
get_knot_points_xy <- function(n_obs, tau) {
  t_points <- length(tau)
  u_points <- floor(n_obs / t_points)
  list(u_points = u_points, t_points = t_points)
}

get_knots_xy <- function(n_obs, tau, u_min, u_max) {
  t_knots <- sort(unique(tau))
  knot_points_xy <- get_knot_points_xy(n_obs, t_knots)
  u_lo <- min(u_min, u_max); u_hi <- max(u_min, u_max)
  u_knots <- seq(u_lo, u_hi, length.out = max(2L, knot_points_xy$u_points))
  list(u = u_knots, t = t_knots)
}

get_knots_z_as_matrix <- function(knots_z, knot_points_xy) {
  matrix(knots_z, nrow = knot_points_xy$u_points, ncol = knot_points_xy$t_points)
}

get_local_volatility <- function(knots_z, knots_xy, grid_xy) {
  sigma_mat <- akima::bicubic(x = knots_xy$u, y = knots_xy$t, z = knots_z,
                               x0 = grid_xy$u, y0 = grid_xy$t) |>
    dplyr::bind_cols() |> reshape2::acast(x ~ y, value.var = "z")
  sigma_mat <- pmin(pmax(sigma_mat, 1e-4), 100)
  sigma_mat
}

get_sigma_interpolated <- function(sigma_knots, knot_points_xy, knots_xy, grid_xy) {
  knots_matrix <- get_knots_z_as_matrix(sigma_knots, knot_points_xy)
  get_local_volatility(knots_matrix, knots_xy, grid_xy)
}

get_grid_xy <- function(u_min, u_max, t_min, t_max, n_u, n_t) {
  u_grid <- seq(u_min, u_max, length = n_u + 1)
  t_grid <- seq(t_min, t_max, length = n_t + 1)
  tidyr::expand_grid(u = u_grid, t = t_grid)
}

# Tikhonov penalty (roughness)
tikhonov_penalty <- function(sigma, weight) {
  if (weight <= 0) return(0)
  R <- nrow(sigma); C <- ncol(sigma)
  pen <- 0; edges <- 0
  for (i in 1:R) for (j in 2:C) { pen <- pen + (sigma[i,j] - sigma[i,j-1])^2; edges <- edges + 1 }
  for (j in 1:C) for (i in 2:R) { pen <- pen + (sigma[i,j] - sigma[i-1,j])^2; edges <- edges + 1 }
  weight * pen / edges
}

# ============================================================================
# STOCK (American) calibration
# ============================================================================
cat("=" |> rep(70) |> paste(collapse=""), "\n")
cat("STOCK LOCAL VOLATILITY CALIBRATION (American options)\n")
cat("=" |> rep(70) |> paste(collapse=""), "\n\n")

data_s <- stock_data_sparse
n_obs_s <- nrow(data_s)
tau_s <- sort(unique(data_s$tau))
knot_points_xy_s <- list(u_points = length(knots_xy_s$u), t_points = length(knots_xy_s$t))

grid_xy_s <- get_grid_xy(s_min, s_max, 0, max(data_s$tau), n_s, n_t)
sigma_init_s <- rep(config$lv$initial_guess, knot_points_xy_s$u_points * knot_points_xy_s$t_points)
tw_s <- tikhonov_weight

cat("Knot grid:", knot_points_xy_s$u_points, "x", knot_points_xy_s$t_points,
    "=", length(sigma_init_s), "parameters\n")
cat("Options:", n_obs_s, "| Tikhonov weight:", tw_s, "\n\n")

# --- Loss: SERIAL (old approach) ---
loss_serial_american <- function(sigma_knots, data, knot_points_xy, knots_xy, grid_xy,
                                  u_min, u_max, n_u, n_t, lambda, tolerance, tikhonov_weight) {
  sigma <- get_sigma_interpolated(sigma_knots, knot_points_xy, knots_xy, grid_xy)
  SCALE <- 1e6
  loss <- 0
  for (i in seq_len(nrow(data))) {
    model <- american_option_lv(data$spot[i], data$k[i], data$tau[i], data$r_d[i], data$q[i],
                                 sigma, data$type[i], u_min, u_max, n_u, n_t, lambda, tolerance)
    loss <- loss + SCALE * data$w[i] * (model - data$price[i])^2
  }
  loss + tikhonov_penalty(sigma, tikhonov_weight)
}

# --- Loss: BATCH (new approach) ---
loss_batch_american <- function(sigma_knots, data, knot_points_xy, knots_xy, grid_xy,
                                 u_min, u_max, n_u, n_t, lambda, tolerance, tikhonov_weight) {
  sigma <- get_sigma_interpolated(sigma_knots, knot_points_xy, knots_xy, grid_xy)
  SCALE <- 1e6
  models <- batch_price_american_lv(data$spot, data$k, data$tau, data$r_d, data$q,
                                     sigma, data$type, u_min, u_max, n_u, n_t, lambda, tolerance)
  diffs <- models - data$price
  loss <- sum(SCALE * data$w * diffs^2)
  loss + tikhonov_penalty(sigma, tikhonov_weight)
}

# --- Gradient (central differences, shared by both) ---
make_gradient <- function(loss_fn, ...) {
  function(sigma_knots, ...) {
    args <- list(...)
    h <- 1e-5
    vapply(seq_along(sigma_knots), function(i) {
      hi <- h * (1 + abs(sigma_knots[i]))
      v_plus <- sigma_knots; v_plus[i] <- v_plus[i] + hi
      v_minus <- sigma_knots; v_minus[i] <- v_minus[i] - hi
      f_plus <- do.call(loss_fn, c(list(v_plus), args))
      f_minus <- do.call(loss_fn, c(list(v_minus), args))
      (f_plus - f_minus) / (2 * hi)
    }, numeric(1))
  }
}

# --- Benchmark single loss + gradient evaluation ---
common_args_s <- list(
  data = data_s, knot_points_xy = knot_points_xy_s, knots_xy = knots_xy_s,
  grid_xy = grid_xy_s, u_min = s_min, u_max = s_max, n_u = n_s, n_t = n_t,
  lambda = lambda, tolerance = tolerance, tikhonov_weight = tw_s
)

cat("--- Benchmarking single loss evaluation ---\n")
t_serial <- system.time(do.call(loss_serial_american, c(list(sigma_init_s), common_args_s)))
t_batch  <- system.time(do.call(loss_batch_american,  c(list(sigma_init_s), common_args_s)))
cat("Serial:", t_serial["elapsed"], "s | Batch:", t_batch["elapsed"], "s | Speedup:",
    round(t_serial["elapsed"] / t_batch["elapsed"], 1), "x\n\n")

cat("--- Benchmarking single gradient evaluation ---\n")
grad_serial <- make_gradient(loss_serial_american)
grad_batch  <- make_gradient(loss_batch_american)

t_gs <- system.time(g_s <- do.call(grad_serial, c(list(sigma_init_s), common_args_s)))
t_gb <- system.time(g_b <- do.call(grad_batch,  c(list(sigma_init_s), common_args_s)))
cat("Serial:", t_gs["elapsed"], "s | Batch:", t_gb["elapsed"], "s | Speedup:",
    round(t_gs["elapsed"] / t_gb["elapsed"], 1), "x\n")
cat("Gradient max diff:", max(abs(g_s - g_b)), "\n\n")

# --- Run optimization with BATCH solver ---
cat("--- Running BATCH American calibration (lbfgsb3) ---\n")
t_opt <- system.time({
  stock_knots_batch <- lbfgsb3c::lbfgsb3(
    sigma_init_s,
    fn = loss_batch_american,
    gr = grad_batch,
    data = data_s, knot_points_xy = knot_points_xy_s, knots_xy = knots_xy_s,
    grid_xy = grid_xy_s, u_min = s_min, u_max = s_max, n_u = n_s, n_t = n_t,
    lambda = lambda, tolerance = tolerance, tikhonov_weight = tw_s,
    lower = 0.001, upper = 5,
    control = list(trace = 1, maxit = 1000)
  )
})
cat("\nOptimization time:", t_opt["elapsed"], "s\n")
cat("Convergence:", stock_knots_batch$convergence, "| Final loss:", stock_knots_batch$f, "\n")

# --- Validate: compute model prices with fitted LV ---
sigma_fitted_s <- get_sigma_interpolated(stock_knots_batch$par, knot_points_xy_s, knots_xy_s, grid_xy_s)

stock_model_prices <- batch_price_american_lv(
  data_s$spot, data_s$k, data_s$tau, data_s$r_d, data_s$q,
  sigma_fitted_s, data_s$type, s_min, s_max, n_s, n_t, lambda, tolerance
)

cat("\n--- Stock fit quality ---\n")
price_errors <- stock_model_prices - data_s$price
rel_errors <- price_errors / data_s$price
cat("Price RMSE:", sqrt(mean(price_errors^2)), "\n")
cat("Mean abs rel error:", mean(abs(rel_errors)) * 100, "%\n")
cat("Max abs rel error:", max(abs(rel_errors)) * 100, "%\n")

# ============================================================================
# FX (European) calibration
# ============================================================================
cat("\n\n")
cat("=" |> rep(70) |> paste(collapse=""), "\n")
cat("FX LOCAL VOLATILITY CALIBRATION (European options)\n")
cat("=" |> rep(70) |> paste(collapse=""), "\n\n")

data_fx <- fx_data_sparse
n_obs_fx <- nrow(data_fx)
tau_fx <- sort(unique(data_fx$tau))
knot_points_xy_fx_local <- list(u_points = length(knots_xy_fx$u), t_points = length(knots_xy_fx$t))
grid_xy_fx_local <- get_grid_xy(x_min, x_max, 0, max(data_fx$tau), n_x, n_t)
sigma_init_fx <- rep(0.15, knot_points_xy_fx_local$u_points * knot_points_xy_fx_local$t_points)
tw_fx <- 10

cat("Knot grid:", knot_points_xy_fx_local$u_points, "x", knot_points_xy_fx_local$t_points,
    "=", length(sigma_init_fx), "parameters\n")
cat("Options:", n_obs_fx, "| Tikhonov weight:", tw_fx, "\n\n")

# --- Loss: BATCH European ---
loss_batch_european <- function(sigma_knots, data, knot_points_xy, knots_xy, grid_xy,
                                 u_min, u_max, n_u, n_t, tikhonov_weight) {
  sigma <- get_sigma_interpolated(sigma_knots, knot_points_xy, knots_xy, grid_xy)
  SCALE <- 1e6
  models <- batch_price_european_lv(data$spot, data$k, data$tau, data$r_d, data$r_f,
                                     sigma, data$type, u_min, u_max, n_u, n_t)
  diffs <- models - data$price
  loss <- sum(SCALE * data$w * diffs^2)
  loss + tikhonov_penalty(sigma, tikhonov_weight)
}

grad_batch_eur <- make_gradient(loss_batch_european)

common_args_fx <- list(
  data = data_fx, knot_points_xy = knot_points_xy_fx_local, knots_xy = knots_xy_fx,
  grid_xy = grid_xy_fx_local, u_min = x_min, u_max = x_max, n_u = n_x, n_t = n_t,
  tikhonov_weight = tw_fx
)

# --- Benchmark ---
cat("--- Benchmarking single loss evaluation ---\n")
t_eur <- system.time(do.call(loss_batch_european, c(list(sigma_init_fx), common_args_fx)))
cat("Batch:", t_eur["elapsed"], "s\n\n")

# --- Run optimization ---
cat("--- Running BATCH European calibration (lbfgsb3) ---\n")
t_opt_fx <- system.time({
  fx_knots_batch <- lbfgsb3c::lbfgsb3(
    sigma_init_fx,
    fn = loss_batch_european,
    gr = grad_batch_eur,
    data = data_fx, knot_points_xy = knot_points_xy_fx_local, knots_xy = knots_xy_fx,
    grid_xy = grid_xy_fx_local, u_min = x_min, u_max = x_max, n_u = n_x, n_t = n_t,
    tikhonov_weight = tw_fx,
    lower = 0.0001, upper = 5,
    control = list(trace = 1, factr = 1e-6, pgtol = 1e-5, maxit = 1000)
  )
})
cat("\nOptimization time:", t_opt_fx["elapsed"], "s\n")
cat("Convergence:", fx_knots_batch$convergence, "| Final loss:", fx_knots_batch$f, "\n")

# --- Validate ---
sigma_fitted_fx <- get_sigma_interpolated(fx_knots_batch$par, knot_points_xy_fx_local, knots_xy_fx, grid_xy_fx_local)

fx_model_prices <- batch_price_european_lv(
  data_fx$spot, data_fx$k, data_fx$tau, data_fx$r_d, data_fx$r_f,
  sigma_fitted_fx, data_fx$type, x_min, x_max, n_x, n_t
)

cat("\n--- FX fit quality ---\n")
price_errors_fx <- fx_model_prices - data_fx$price
nonzero <- data_fx$price > 1e-6
rel_errors_fx <- price_errors_fx[nonzero] / data_fx$price[nonzero]
cat("Price RMSE:", sqrt(mean(price_errors_fx^2)), "\n")
cat("Mean abs rel error:", mean(abs(rel_errors_fx)) * 100, "%\n")
cat("Max abs rel error:", max(abs(rel_errors_fx)) * 100, "%\n")

cat("\n\n--- SUMMARY ---\n")
cat("Stock (American) optimization:", round(t_opt["elapsed"], 1), "s\n")
cat("FX (European) optimization:", round(t_opt_fx["elapsed"], 1), "s\n")
cat("Total:", round(t_opt["elapsed"] + t_opt_fx["elapsed"], 1), "s\n")
