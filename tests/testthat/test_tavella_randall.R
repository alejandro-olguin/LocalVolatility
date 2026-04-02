test_that("Tavella-Randall endpoints are exact", {
  grid <- tavella_randall(100, 3.0, 10, 300, 50)
  expect_equal(grid[1], 10)
  expect_equal(grid[51], 300)
  expect_length(grid, 51)  # n+1 points
})

test_that("Tavella-Randall grid is monotonically increasing", {
  grid <- tavella_randall(100, 3.0, 10, 300, 50)
  diffs <- diff(grid)
  expect_true(all(diffs > 0))
})

test_that("Tavella-Randall concentrates nodes near x0", {
  grid <- tavella_randall(100, 3.0, 10, 300, 100)
  diffs <- diff(grid)

  # Find the step containing x0=100
  idx_near <- which(grid >= 100)[1]
  step_near <- diffs[max(1, idx_near - 1)]

  # Step near boundaries should be larger
  step_boundary <- diffs[1]
  expect_gt(step_boundary, step_near)
})

test_that("Tavella-Randall rejects invalid inputs", {
  expect_error(tavella_randall(100, 3.0, 300, 10, 50))   # x_max <= x_min
  expect_error(tavella_randall(100, -1.0, 10, 300, 50))  # alpha <= 0
  expect_error(tavella_randall(100, 3.0, 10, 300, 0))    # n < 1
})
