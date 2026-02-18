# Tests for growth module functions
library(RTMB)

# Define project paths
proj_dir <- this.path::this.proj()
r_dir <- file.path(proj_dir, "code", "rtmb", "R")

# Source helper functions
suppressMessages({
  source(file.path(r_dir, "growth.R"))
  source(file.path(r_dir, "selectivity.R"))
  source(file.path(r_dir, "functions.R"))
})

# Shared test fixtures
len_lower <- seq(10, by = 2, length.out = 95)
len_upper <- len_lower + 2
len_mid   <- (len_lower + len_upper) / 2
n_age <- 40L
A1 <- 1L
A2 <- 40L
L1 <- 30.0
L2 <- 180.0
k  <- 0.2
CV1 <- 0.15
CV2 <- 0.08

# Test get_growth ----

test_that("get_growth returns vector of length n_age", {
  mu_a <- get_growth(n_age, A1, A2, L1, L2, k)
  expect_equal(length(mu_a), n_age)
})

test_that("get_growth returns L1 at age A1 and L2 at age A2", {
  mu_a <- get_growth(n_age, A1, A2, L1, L2, k)
  expect_equal(mu_a[A1], L1, tolerance = 1e-10)
  expect_equal(mu_a[A2], L2, tolerance = 1e-10)
})

test_that("get_growth is monotonically increasing for k > 0 and L2 > L1", {
  mu_a <- get_growth(n_age, A1, A2, L1, L2, k)
  expect_true(all(diff(mu_a) > 0))
})

test_that("get_growth output values are between L1 and L2", {
  mu_a <- get_growth(n_age, A1, A2, L1, L2, k)
  expect_true(all(mu_a >= L1 - 1e-10))
  expect_true(all(mu_a <= L2 + 1e-10))
})

# Test get_sd_at_age ----

test_that("get_sd_at_age returns correct SD at reference ages", {
  mu_a <- get_growth(n_age, A1, A2, L1, L2, k)
  sd_a <- get_sd_at_age(mu_a, L1, L2, CV1, CV2)
  # At age A1: mu = L1, so SD = L1 * CV1
  expect_equal(sd_a[A1], L1 * CV1, tolerance = 1e-10)
  # At age A2: mu = L2, so SD = L2 * CV2
  expect_equal(sd_a[A2], L2 * CV2, tolerance = 1e-10)
})

test_that("get_sd_at_age returns vector of length n_age", {
  mu_a <- get_growth(n_age, A1, A2, L1, L2, k)
  sd_a <- get_sd_at_age(mu_a, L1, L2, CV1, CV2)
  expect_equal(length(sd_a), n_age)
})

test_that("get_sd_at_age returns positive values", {
  mu_a <- get_growth(n_age, A1, A2, L1, L2, k)
  sd_a <- get_sd_at_age(mu_a, L1, L2, CV1, CV2)
  expect_true(all(sd_a > 0))
})

# Test get_weight_at_length ----

test_that("get_weight_at_length returns positive values for positive lengths", {
  lw_a <- 1.15 * 2.942e-06
  lw_b <- 3.13088
  wt <- get_weight_at_length(len_mid, lw_a, lw_b)
  expect_true(all(wt > 0))
})

test_that("get_weight_at_length returns vector of same length as len_mid", {
  lw_a <- 1.15 * 2.942e-06
  lw_b <- 3.13088
  wt <- get_weight_at_length(len_mid, lw_a, lw_b)
  expect_equal(length(wt), length(len_mid))
})

test_that("get_weight_at_length increases with length (b > 0)", {
  lw_a <- 1.15 * 2.942e-06
  lw_b <- 3.13088
  wt <- get_weight_at_length(len_mid, lw_a, lw_b)
  expect_true(all(diff(wt) > 0))
})

# Test get_maturity_at_age ----

test_that("get_maturity_at_age returns values in [0, 1]", {
  mu_a <- get_growth(n_age, A1, A2, L1, L2, k)
  sd_a <- get_sd_at_age(mu_a, L1, L2, CV1, CV2)
  pla  <- get_pla(len_lower, len_upper, mu_a, sd_a)
  # Simple knife-edge maturity at 100 cm
  mat_l <- as.numeric(len_mid >= 100)
  mat_a <- get_maturity_at_age(pla, mat_l)
  expect_true(all(mat_a >= 0 - 1e-10))
  expect_true(all(mat_a <= 1 + 1e-10))
})

test_that("get_maturity_at_age returns vector of length n_age", {
  mu_a <- get_growth(n_age, A1, A2, L1, L2, k)
  sd_a <- get_sd_at_age(mu_a, L1, L2, CV1, CV2)
  pla  <- get_pla(len_lower, len_upper, mu_a, sd_a)
  mat_l <- as.numeric(len_mid >= 100)
  mat_a <- get_maturity_at_age(pla, mat_l)
  expect_equal(length(mat_a), n_age)
})

test_that("get_maturity_at_age is monotonically non-decreasing for knife-edge maturity", {
  mu_a <- get_growth(n_age, A1, A2, L1, L2, k)
  sd_a <- get_sd_at_age(mu_a, L1, L2, CV1, CV2)
  pla  <- get_pla(len_lower, len_upper, mu_a, sd_a)
  mat_l <- as.numeric(len_mid >= 100)
  mat_a <- get_maturity_at_age(pla, mat_l)
  expect_true(all(diff(mat_a) >= -1e-10))
})

# Test full pipeline ----

test_that("get_growth at A1 and A2 returns exactly L1 and L2 (boundary conditions)", {
  # Explicit boundary test with different parameter values
  mu_a2 <- get_growth(50L, 5L, 45L, 40.0, 200.0, 0.15)
  expect_equal(mu_a2[5],  40.0,  tolerance = 1e-10)
  expect_equal(mu_a2[45], 200.0, tolerance = 1e-10)
})

test_that("full pipeline: weight-at-age is positive and increasing over young ages", {
  mu_a <- get_growth(n_age, A1, A2, L1, L2, k)
  sd_a <- get_sd_at_age(mu_a, L1, L2, CV1, CV2)
  pla  <- get_pla(len_lower, len_upper, mu_a, sd_a)
  lw_a_par <- 1.15 * 2.942e-06
  lw_b_par <- 3.13088
  wt_at_len <- get_weight_at_length(len_mid, lw_a_par, lw_b_par)
  weight_a  <- as.vector(t(pla) %*% wt_at_len)
  expect_equal(length(weight_a), n_age)
  expect_true(all(weight_a > 0))
  # Weight should be increasing over first half of ages (young fish growing fast)
  expect_true(all(diff(weight_a[1:20]) > 0))
})

test_that("full pipeline: maturity-at-age is consistent with direct pla conversion", {
  mu_a <- get_growth(n_age, A1, A2, L1, L2, k)
  sd_a <- get_sd_at_age(mu_a, L1, L2, CV1, CV2)
  pla  <- get_pla(len_lower, len_upper, mu_a, sd_a)
  mat_l <- pmin(pmax(len_mid / 150, 0), 1)  # simple ramp maturity
  # Direct conversion
  mat_a_direct <- as.vector(t(pla) %*% mat_l)
  # Via get_maturity_at_age
  mat_a_module <- get_maturity_at_age(pla, mat_l)
  expect_equal(mat_a_module, mat_a_direct, tolerance = 1e-12)
})
