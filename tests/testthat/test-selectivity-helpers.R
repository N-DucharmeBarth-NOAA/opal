# Tests for selectivity helper functions
library(RTMB)

# Define project paths
proj_dir <- this.path::this.proj()
r_dir <- file.path(proj_dir, "code", "rtmb", "R")

# Source helper functions
suppressMessages({
  source(file.path(r_dir, "selectivity.R"))
  source(file.path(r_dir, "functions.R"))
  source(file.path(r_dir, "length-weight.R"))
})

# Test sel_logistic ----

test_that("sel_logistic returns ~0.5 at inflection point when a = 0", {
  lens <- seq(50, 200, by = 1)  # Use finer bins to get closer to exact inflection
  # a = 0 places inflection at mean(lens)
  sel <- sel_logistic(lens, par = c(0, 0, 0, 0, 0, 0))
  mean_len <- mean(lens)
  idx <- which.min(abs(lens - mean_len))
  expect_equal(sel[idx], 0.5, tolerance = 1e-3)  # Allow small tolerance for bin discretization
})

test_that("sel_logistic values are in (0, 1]", {
  lens <- seq(10, 250, by = 1)
  # Use real-line parameterization: a = -2 shifts peak left, b = 1 increases width
  sel <- sel_logistic(lens, par = c(-2, 1, 0, 0, 0, 0))
  expect_true(all(sel > 0 & sel <= 1))
})

test_that("sel_logistic approaches 1 for large lengths", {
  lens <- seq(50, 300, by = 2)
  sel <- sel_logistic(lens, par = c(-2, 1, 0, 0, 0, 0))
  expect_true(sel[length(sel)] > 0.95)
})

test_that("sel_logistic approaches 0 for small lengths", {
  lens <- seq(10, 250, by = 2)
  sel <- sel_logistic(lens, par = c(2, 1, 0, 0, 0, 0))
  expect_true(sel[1] < 0.05)
})

test_that("sel_logistic output length matches input", {
  lens <- seq(50, 200, by = 5)
  sel <- sel_logistic(lens, par = c(0, 0.5, 0, 0, 0, 0))
  expect_equal(length(sel), length(lens))
})

test_that("sel_logistic is monotonically increasing", {
  lens <- seq(50, 200, by = 1)
  sel <- sel_logistic(lens, par = c(0, 0.5, 0, 0, 0, 0))
  expect_true(all(diff(sel) > 0))
})

test_that("sel_logistic inflection point shifts with parameter a", {
  lens <- seq(50, 200, by = 1)
  mu <- mean(lens)
  sd <- sd(lens)
  
  # a = -1 shifts inflection left, a = 1 shifts inflection right
  sel_left <- sel_logistic(lens, par = c(-1, 0.5, 0, 0, 0, 0))
  sel_right <- sel_logistic(lens, par = c(1, 0.5, 0, 0, 0, 0))
  
  # Find the index closest to 0.5 for each curve (inflection point)
  idx_infl_left <- which.min(abs(sel_left - 0.5))
  idx_infl_right <- which.min(abs(sel_right - 0.5))
  
  # The inflection point should shift right
  expect_true(lens[idx_infl_left] < lens[idx_infl_right])
})

test_that("sel_logistic width increases with parameter b", {
  lens <- seq(50, 200, by = 1)
  
  # b = 0 gives narrow curve, b = 1 gives wider curve
  sel_narrow <- sel_logistic(lens, par = c(0, 0, 0, 0, 0, 0))
  sel_wide <- sel_logistic(lens, par = c(0, 1, 0, 0, 0, 0))
  
  # At inflection point, width is determined by b
  # The wide curve should have a shallower slope
  mean_len <- mean(lens)
  idx <- which.min(abs(lens - mean_len))
  
  # The wide curve's slope should be smaller
  slope_narrow <- diff(sel_narrow[c(idx, idx + 1)])
  slope_wide <- diff(sel_wide[c(idx, idx + 1)])
  expect_true(slope_narrow > slope_wide)
})

# Test sel_double_normal ----

test_that("sel_double_normal produces valid curve with real-line parameters", {
  lens <- seq(20, 200, by = 2)
  # All parameters on real line: a=0 (peak at mean), b=0 (neutral plateau), 
  # c=0, d=0 (neutral widths), e=-9, f=-9 (very low initial/final selectivity)
  sel <- sel_double_normal(lens, par = c(0, 0, 0, 0, -9, -9))
  expect_true(all(sel >= 0 & sel <= 1, na.rm = TRUE))
  expect_equal(length(sel), length(lens))
})

test_that("sel_double_normal peak is near maximum", {
  lens <- seq(20, 200, by = 2)
  # Peak location a = 0 places peak near mean(lens)
  sel <- sel_double_normal(lens, par = c(0, 0, 1, 1, -9, -9))
  mean_len <- mean(lens)
  peak_idx <- which.min(abs(lens - mean_len))
  expect_true(sel[peak_idx] > 0.9)
})

test_that("sel_double_normal ascending limb increases up to peak", {
  lens <- seq(20, 200, by = 2)
  sel <- sel_double_normal(lens, par = c(0, 0, 1, 1, -9, -9))
  mean_len <- mean(lens)
  peak_idx <- which.min(abs(lens - mean_len))
  asc <- sel[1:peak_idx]
  # Allow small numerical tolerance in monotonicity check
  expect_true(all(diff(asc) >= -1e-10))
})

test_that("sel_double_normal works with moderate initial/final selectivity", {
  lens <- seq(20, 200, by = 2)
  # e = -2 (logistic ~ 0.12), f = -1 (logistic ~ 0.27)
  sel <- sel_double_normal(lens, par = c(0, 0, 1, 1, -2, -1))
  expect_true(all(sel >= 0 & sel <= 1, na.rm = TRUE))
  # With e = -2, initial selectivity should be non-trivial
  expect_true(sel[1] > 0.05)
  expect_true(sel[1] < 0.25)
})

test_that("sel_double_normal output length matches input", {
  lens <- seq(30, 180, by = 5)
  sel <- sel_double_normal(lens, par = c(0, 0, 1, 1, -9, -9))
  expect_equal(length(sel), length(lens))
})

test_that("sel_double_normal peak location shifts with parameter a", {
  lens <- seq(20, 200, by = 2)
  mu <- mean(lens)
  sd <- sd(lens)
  
  # a = -1 shifts peak left, a = 1 shifts peak right
  sel_left <- sel_double_normal(lens, par = c(-1, 0, 1, 1, -9, -9))
  sel_right <- sel_double_normal(lens, par = c(1, 0, 1, 1, -9, -9))
  
  idx_max_left <- which.max(sel_left)
  idx_max_right <- which.max(sel_right)
  
  expect_true(lens[idx_max_left] < lens[idx_max_right])
})

test_that("sel_double_normal ascending width changes with parameter c", {
  lens <- seq(20, 200, by = 2)
  
  # c = 0 (narrow), c = 1 (wide)
  sel_narrow <- sel_double_normal(lens, par = c(0, 0, 0, 1, -9, -9))
  sel_wide <- sel_double_normal(lens, par = c(0, 0, 1, 1, -9, -9))
  
  # The curve with c = 1 should have a broader ascending limb
  mean_len <- mean(lens)
  peak_idx <- which.min(abs(lens - mean_len))
  
  # At half the peak range, the narrow curve should be lower
  half_peak_idx <- peak_idx / 2
  expect_true(sel_narrow[half_peak_idx] < sel_wide[half_peak_idx])
})

# Test get_pla ----

test_that("get_pla columns sum to 1", {
  len_lower <- seq(20, 198, by = 2)
  len_upper <- seq(22, 200, by = 2)
  mu_a <- 30 + (180 - 30) * (1 - exp(-0.2 * (0:30)))
  sd_a <- 0.1 * mu_a
  pla <- get_pla(len_lower, len_upper, mu_a, sd_a)
  col_sums <- colSums(pla)
  expect_equal(col_sums, rep(1, length(mu_a)), tolerance = 1e-4)
})

test_that("get_pla values are non-negative", {
  len_lower <- seq(20, 198, by = 2)
  len_upper <- seq(22, 200, by = 2)
  mu_a <- 30 + (180 - 30) * (1 - exp(-0.2 * (0:30)))
  sd_a <- 0.1 * mu_a
  pla <- get_pla(len_lower, len_upper, mu_a, sd_a)
  expect_true(all(pla >= 0))
})

test_that("get_pla has correct dimensions", {
  len_lower <- seq(20, 198, by = 2)
  len_upper <- seq(22, 200, by = 2)
  mu_a <- 30 + (180 - 30) * (1 - exp(-0.2 * (0:30)))
  sd_a <- 0.1 * mu_a
  pla <- get_pla(len_lower, len_upper, mu_a, sd_a)
  expect_equal(dim(pla), c(length(len_lower), length(mu_a)))
})

test_that("get_pla concentrates probability around mean length at age", {
  len_lower <- seq(20, 198, by = 2)
  len_upper <- seq(22, 200, by = 2)
  len_mid <- (len_lower + len_upper) / 2
  mu_a <- 30 + (180 - 30) * (1 - exp(-0.2 * (0:30)))
  sd_a <- 0.1 * mu_a
  pla <- get_pla(len_lower, len_upper, mu_a, sd_a)
  # For each age, the bin containing the mean length should have high probability
  for (a in seq_along(mu_a)) {
    closest_bin <- which.min(abs(len_mid - mu_a[a]))
    expect_true(pla[closest_bin, a] > 0.01,
                info = sprintf("Age %d: expected high probability near mean length %.1f", a, mu_a[a]))
  }
})
