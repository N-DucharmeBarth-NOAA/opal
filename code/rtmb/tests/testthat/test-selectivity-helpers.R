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

test_that("sel_logistic returns 0.5 at inflection point", {
  lens <- seq(50, 200, by = 2)
  sel <- sel_logistic(lens, a = 120, b = 15)
  expect_equal(sel[which(lens == 120)], 0.5, tolerance = 1e-6)
})

test_that("sel_logistic values are in (0, 1]", {
  lens <- seq(10, 250, by = 1)
  sel <- sel_logistic(lens, a = 120, b = 15)
  expect_true(all(sel > 0 & sel <= 1))
})

test_that("sel_logistic approaches 1 for large lengths", {
  lens <- seq(50, 300, by = 2)
  sel <- sel_logistic(lens, a = 120, b = 15)
  expect_true(sel[length(sel)] > 0.999)
})

test_that("sel_logistic approaches 0 for small lengths", {
  lens <- seq(10, 250, by = 2)
  sel <- sel_logistic(lens, a = 120, b = 15)
  expect_true(sel[1] < 0.001)
})

test_that("sel_logistic output length matches input", {
  lens <- seq(50, 200, by = 5)
  sel <- sel_logistic(lens, a = 100, b = 20)
  expect_equal(length(sel), length(lens))
})

test_that("sel_logistic is monotonically increasing", {
  lens <- seq(50, 200, by = 1)
  sel <- sel_logistic(lens, a = 120, b = 15)
  expect_true(all(diff(sel) > 0))
})

# Test sel_double_normal ----

test_that("sel_double_normal produces valid curve with logit-space e and f", {
  lens <- seq(20, 200, by = 2)
  sel <- sel_double_normal(lens, a = 120, b = 0, c = 3, d = 4, e = -9, f = -9)
  expect_true(all(sel >= 0 & sel <= 1, na.rm = TRUE))
  expect_equal(length(sel), length(lens))
})

test_that("sel_double_normal peak is near maximum", {
  lens <- seq(20, 200, by = 2)
  sel <- sel_double_normal(lens, a = 120, b = 0, c = 3, d = 4, e = -9, f = -9)
  peak_idx <- which(lens == 120)
  expect_true(sel[peak_idx] > 0.9)
})

test_that("sel_double_normal ascending limb increases up to peak", {
  lens <- seq(20, 200, by = 2)
  sel <- sel_double_normal(lens, a = 120, b = 0, c = 3, d = 4, e = -9, f = -9)
  peak_idx <- which(lens == 120)
  asc <- sel[1:peak_idx]
  # Allow small numerical tolerance in monotonicity check
  expect_true(all(diff(asc) >= -1e-10))
})

test_that("sel_double_normal works with moderate e and f", {
  lens <- seq(20, 200, by = 2)
  sel <- sel_double_normal(lens, a = 130, b = 0, c = 3, d = 3, e = -2, f = -1)
  expect_true(all(sel >= 0 & sel <= 1, na.rm = TRUE))
  # With e = -2 (logistic ~ 0.12), initial selectivity should be non-trivial
  expect_true(sel[1] > 0.05)
})

test_that("sel_double_normal output length matches input", {
  lens <- seq(30, 180, by = 5)
  sel <- sel_double_normal(lens, a = 100, b = 0, c = 3, d = 4, e = -9, f = -9)
  expect_equal(length(sel), length(lens))
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
