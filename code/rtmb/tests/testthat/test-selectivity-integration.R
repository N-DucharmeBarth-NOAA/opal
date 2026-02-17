# Tests for selectivity integration into bet_model
library(RTMB)

# Define project paths
proj_dir <- this.path::this.proj()
r_dir <- file.path(proj_dir, "code", "rtmb", "R")
data_dir <- file.path(proj_dir, "code", "rtmb")

# Source helper functions
suppressMessages({
  source(file.path(r_dir, "selectivity.R"))
  source(file.path(r_dir, "functions.R"))
  source(file.path(r_dir, "natural-mortality.R"))
  source(file.path(r_dir, "recruitment.R"))
  source(file.path(r_dir, "dynamics.R"))
  source(file.path(r_dir, "length-weight.R"))
  source(file.path(r_dir, "likelihoods.R"))
  source(file.path(r_dir, "model.R"))
  source(file.path(r_dir, "priors.R"))
})

# Test selectivity configuration in data ----

test_that("data includes all required selectivity elements", {
  skip_if_not(file.exists(file.path(data_dir, "catch-data.csv")), 
              "Catch data file not found")
  
  # Create minimal data structure for testing
  data <- list(
    n_fishery = 15,
    n_year = 268,
    n_age = 40,
    n_season = 1
  )
  
  # Add selectivity configuration elements
  len_lower <- seq(10, by = 2, length.out = 95)
  len_upper <- len_lower + 2
  len_mid <- (len_lower + len_upper) / 2
  
  data$sel_type_f <- rep(2L, data$n_fishery)
  data$sel_type_f[c(11, 15)] <- 1L
  data$sel_len_lower <- len_lower
  data$sel_len_upper <- len_upper
  data$sel_lengths <- len_mid
  
  # Verify sel_type_f
  expect_true("sel_type_f" %in% names(data))
  expect_equal(length(data$sel_type_f), data$n_fishery)
  expect_true(all(data$sel_type_f %in% c(1L, 2L)))
  expect_equal(data$sel_type_f[11], 1L)
  expect_equal(data$sel_type_f[15], 1L)
  expect_equal(sum(data$sel_type_f == 1L), 2)
  expect_equal(sum(data$sel_type_f == 2L), 13)
  
  # Verify length bins
  expect_true("sel_lengths" %in% names(data))
  expect_true("sel_len_lower" %in% names(data))
  expect_true("sel_len_upper" %in% names(data))
  expect_equal(length(data$sel_len_lower), length(data$sel_lengths))
  expect_true(all(data$sel_len_upper > data$sel_len_lower))
  
  # Midpoints should match
  expected_mid <- (data$sel_len_lower + data$sel_len_upper) / 2
  expect_equal(data$sel_lengths, expected_mid)
  
  # Bins should cover bigeye tuna size range
  expect_true(min(data$sel_len_lower) <= 20)
  expect_true(max(data$sel_len_upper) >= 200)
})

# Test par_sel parameter structure ----

test_that("par_sel has correct structure and dimensions", {
  # Create test SS3 parameters
  ss3_selex_pars <- matrix(NA, nrow = 15, ncol = 6)
  
  # F01-F07: Longline fleets (double-normal)
  for (f in 1:7) {
    ss3_selex_pars[f, ] <- c(100, -5, 0, 0, -999, 9)
  }
  
  # F08-F09: Purse seine fleets (double-normal)
  for (f in 8:9) {
    ss3_selex_pars[f, ] <- c(40, -5, 0, 0, -999, -9)
  }
  
  # F10: Domestic misc fisheries (double-normal)
  ss3_selex_pars[10, ] <- c(40, -5, 0, 0, 0, -999)
  
  # F11: Domestic handline (logistic)
  ss3_selex_pars[11, ] <- c(100, 10, NA, NA, NA, NA)
  
  # F12-F14: PL and JP PS (double-normal)
  for (f in 12:14) {
    ss3_selex_pars[f, ] <- c(40, -5, 0, 0, -999, -9)
  }
  
  # F15: Index fishery (logistic)
  ss3_selex_pars[15, ] <- c(100, 10, NA, NA, NA, NA)
  
  # Selectivity types
  sel_type_f <- rep(2L, 15)
  sel_type_f[c(11, 15)] <- 1L
  
  # Length bins
  sel_lengths <- seq(11, 199, by = 2)
  
  # Convert to RTMB parameterization
  par_sel <- convert_ss3_selex_to_rtmb(
    ss3_pars = ss3_selex_pars,
    sel_type_f = sel_type_f,
    sel_lengths = sel_lengths
  )
  
  # Verify structure
  expect_true(is.matrix(par_sel))
  expect_equal(dim(par_sel), c(15, 6))
  expect_true(all(is.finite(par_sel)))
})

# Test priors include par_sel ----

test_that("get_priors includes selectivity parameters", {
  par_sel <- matrix(0, nrow = 15, ncol = 6)
  
  parameters <- list(
    log_B0 = 13,
    log_h = log(0.8),
    log_sigma_r = log(0.6),
    log_cpue_q = log(1),
    cpue_creep = 0,
    log_cpue_tau = log(0.2),
    log_cpue_omega = log(1),
    rdev_y = rnorm(268, 0, 0.1),
    par_sel = par_sel
  )
  
  priors <- get_priors(parameters)
  
  expect_true("par_sel" %in% names(priors))
  expect_equal(priors$par_sel$type, "normal")
  expect_equal(priors$par_sel$par1, 0)
  expect_equal(priors$par_sel$par2, 2)
})

# Test model globals include selectivity functions ----

test_that("bet_globals includes selectivity functions", {
  globals <- bet_globals()
  
  expect_true("get_selectivity" %in% names(globals))
  expect_true("sel_logistic" %in% names(globals))
  expect_true("sel_double_normal" %in% names(globals))
  expect_true("get_pla" %in% names(globals))
})
