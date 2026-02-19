# Tests for get_selectivity and SS3 parameter conversion
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

# Test get_selectivity ----

test_that("get_selectivity returns correct dimensions", {
  mu_a <- 30 + (180 - 30) * (1 - exp(-0.2 * (0:39)))
  sd_a <- 0.1 * mu_a
  len_lower <- seq(10, 198, by = 2)
  len_upper <- seq(12, 200, by = 2)
  len_mid   <- seq(11, 199, by = 2)
  data <- list(
    n_fishery = 3, n_year = 10, n_age = 40,
    sel_type_f = c(1L, 2L, 1L)
  )
  par_sel <- matrix(0, nrow = 3, ncol = 6)
  par_sel[1, ] <- c(0, 0, 0, 0, 0, 0)           # logistic
  par_sel[2, ] <- c(0, 0, 0, 0, -9, -9)          # double-normal
  par_sel[3, ] <- c(-0.5, 0.5, 0, 0, 0, 0)       # logistic
  sel <- get_selectivity(data, par_sel, mu_a, sd_a, len_lower, len_upper, len_mid)
  expect_equal(dim(sel), c(3, 10, 40))
})

test_that("get_selectivity values are in [0, 1]", {
  mu_a <- 30 + (180 - 30) * (1 - exp(-0.2 * (0:39)))
  sd_a <- 0.1 * mu_a
  len_lower <- seq(10, 198, by = 2)
  len_upper <- seq(12, 200, by = 2)
  len_mid   <- seq(11, 199, by = 2)
  data <- list(
    n_fishery = 2, n_year = 5, n_age = 40,
    sel_type_f = c(1L, 2L)
  )
  par_sel <- matrix(0, nrow = 2, ncol = 6)
  par_sel[1, ] <- c(0, 0, 0, 0, 0, 0)
  par_sel[2, ] <- c(0, 0, 1, 1, -9, -9)
  sel <- get_selectivity(data, par_sel, mu_a, sd_a, len_lower, len_upper, len_mid)
  expect_true(all(sel >= 0 & sel <= 1))
})

test_that("get_selectivity produces valid selectivity values", {
  mu_a <- 30 + (180 - 30) * (1 - exp(-0.2 * (0:39)))
  sd_a <- 0.1 * mu_a
  len_lower <- seq(10, 198, by = 2)
  len_upper <- seq(12, 200, by = 2)
  len_mid   <- seq(11, 199, by = 2)
  data <- list(
    n_fishery = 3, n_year = 10, n_age = 40,
    sel_type_f = c(1L, 2L, 2L)
  )
  par_sel <- matrix(0, nrow = 3, ncol = 6)
  par_sel[1, ] <- c(0, 0, 0, 0, 0, 0)
  par_sel[2, ] <- c(0, 0, 1, 1, -9, -9)
  par_sel[3, ] <- c(0.5, -2, 0.5, 0.5, -5, -5)
  sel <- get_selectivity(data, par_sel, mu_a, sd_a, len_lower, len_upper, len_mid)
  # Values should be in [0, 1]
  expect_true(all(sel >= 0 & sel <= 1))
})

test_that("get_selectivity is time-invariant", {
  mu_a <- 30 + (180 - 30) * (1 - exp(-0.2 * (0:39)))
  sd_a <- 0.1 * mu_a
  len_lower <- seq(10, 198, by = 2)
  len_upper <- seq(12, 200, by = 2)
  len_mid   <- seq(11, 199, by = 2)
  data <- list(
    n_fishery = 2, n_year = 10, n_age = 40,
    sel_type_f = c(1L, 2L)
  )
  par_sel <- matrix(0, nrow = 2, ncol = 6)
  par_sel[1, ] <- c(0, 0, 0, 0, 0, 0)
  par_sel[2, ] <- c(0, 0, 1, 1, -9, -9)
  sel <- get_selectivity(data, par_sel, mu_a, sd_a, len_lower, len_upper, len_mid)
  for (f in 1:2) {
    expect_equal(sel[f, 1, ], sel[f, 5, ])
    expect_equal(sel[f, 1, ], sel[f, 10, ])
  }
})

test_that("logistic selectivity-at-age is monotonically increasing", {
  mu_a <- 30 + (180 - 30) * (1 - exp(-0.2 * (0:39)))
  sd_a <- 0.1 * mu_a
  len_lower <- seq(10, 198, by = 2)
  len_upper <- seq(12, 200, by = 2)
  len_mid   <- seq(11, 199, by = 2)
  data <- list(
    n_fishery = 1, n_year = 5, n_age = 40,
    sel_type_f = c(1L)
  )
  par_sel <- matrix(c(0, 0, 0, 0, 0, 0), nrow = 1)
  sel <- get_selectivity(data, par_sel, mu_a, sd_a, len_lower, len_upper, len_mid)
  sel_a <- sel[1, 1, ]
  expect_true(all(diff(sel_a) >= -1e-10))
})

# Note: AD gradient test removed. Future capability when growth parameters are estimated.
# Currently mu_a and sd_a are fixed, so only par_sel flows through the AD tape.

# Test convert_ss3_selex_to_rtmb ----

test_that("SS3 to RTMB conversion round-trips correctly for double-normal", {
  sel_lengths <- seq(11, 199, by = 2)

  # SS3 parameters: peak=120, top_logit=-5, ascend_se=3, descend_se=4, start=-9, end=9
  ss3_pars <- matrix(c(120, -5, 3, 4, -9, 9), nrow = 1)
  sel_type_f <- 2L

  rtmb_pars <- convert_ss3_selex_to_rtmb(ss3_pars, sel_type_f, sel_lengths)
  ss3_back <- convert_rtmb_selex_to_ss3(rtmb_pars, sel_type_f, sel_lengths)

  expect_equal(ss3_back$peak_or_inflection, 120, tolerance = 1e-6)
  expect_equal(ss3_back$top_logit_or_width, -5, tolerance = 1e-6)
  expect_equal(ss3_back$ascend_se, 3, tolerance = 1e-6)
  expect_equal(ss3_back$descend_se, 4, tolerance = 1e-6)
  expect_equal(ss3_back$start_logit, -9, tolerance = 1e-6)
  expect_equal(ss3_back$end_logit, 9, tolerance = 1e-6)
})

test_that("SS3 to RTMB conversion round-trips correctly for logistic", {
  sel_lengths <- seq(11, 199, by = 2)

  # SS3 logistic: inflection=100, width=10
  ss3_pars <- matrix(c(100, 10), nrow = 1)
  sel_type_f <- 1L

  rtmb_pars <- convert_ss3_selex_to_rtmb(ss3_pars, sel_type_f, sel_lengths)
  ss3_back <- convert_rtmb_selex_to_ss3(rtmb_pars, sel_type_f, sel_lengths)

  expect_equal(ss3_back$peak_or_inflection, 100, tolerance = 1e-6)
  expect_equal(ss3_back$top_logit_or_width, 10, tolerance = 1e-6)
})

test_that("SS3 start_logit = -999 converts to RTMB e = -9", {
  sel_lengths <- seq(11, 199, by = 2)
  ss3_pars <- matrix(c(100, -5, 3, 4, -999, -9), nrow = 1)
  sel_type_f <- 2L

  rtmb_pars <- convert_ss3_selex_to_rtmb(ss3_pars, sel_type_f, sel_lengths)
  expect_equal(rtmb_pars[1, 5], -9)
})

test_that("Converted SS3 parameters produce valid selectivity curves", {
  sel_lengths <- seq(11, 199, by = 2)
  len_lower <- sel_lengths - 1
  len_upper <- sel_lengths + 1
  mu_a <- 30 + (180 - 30) * (1 - exp(-0.2 * (0:39)))
  sd_a <- 0.1 * mu_a

  # Typical SS3 longline parameters
  ss3_pars <- matrix(c(100, -5, 0, 0, -999, 9), nrow = 1)
  sel_type_f <- 2L

  rtmb_pars <- convert_ss3_selex_to_rtmb(ss3_pars, sel_type_f, sel_lengths)

  data <- list(
    n_fishery = 1, n_year = 5, n_age = 40,
    sel_type_f = sel_type_f
  )
  sel <- get_selectivity(data, rtmb_pars, mu_a, sd_a, len_lower, len_upper, sel_lengths)
  # Values should be in [0, 1] (no longer normalized to max=1 for AD compatibility)
  expect_true(all(sel >= 0 & sel <= 1))
})

test_that("Mixed logistic and double-normal fisheries convert correctly", {
  sel_lengths <- seq(11, 199, by = 2)

  # 3 fisheries: double-normal, logistic, double-normal
  ss3_pars <- matrix(NA, nrow = 3, ncol = 6)
  ss3_pars[1, ] <- c(100, -5, 3, 4, -999, 9)
  ss3_pars[2, ] <- c(100, 10, NA, NA, NA, NA)   # logistic: only first 2 used
  ss3_pars[3, ] <- c(40, -5, 0, 0, -999, -9)
  sel_type_f <- c(2L, 1L, 2L)

  rtmb_pars <- convert_ss3_selex_to_rtmb(ss3_pars, sel_type_f, sel_lengths)

  expect_equal(nrow(rtmb_pars), 3)
  expect_equal(ncol(rtmb_pars), 6)
  # Logistic fleet should have zeros in columns 3:6
  expect_equal(rtmb_pars[2, 3:6], c(0, 0, 0, 0))
})
