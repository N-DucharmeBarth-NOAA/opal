# Tests for growth module functions
library(RTMB)

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

# Test MFCL parameter matching ----

test_that("maturity_a * fecundity_a equals maturity_a * weight_a when fecundity = weight", {
  mu_a <- get_growth(n_age, A1, A2, L1, L2, k)
  sd_a <- get_sd_at_age(mu_a, L1, L2, CV1, CV2)
  pla  <- get_pla(len_lower, len_upper, mu_a, sd_a)
  lw_a_coef <- 6.48e-05
  lw_b_coef <- 2.78
  wt_l  <- get_weight_at_length(len_mid, lw_a_coef, lw_b_coef)
  mat_a <- get_maturity_at_age(pla, as.numeric(len_mid >= 100))
  wt_a  <- as.vector(t(pla) %*% wt_l)
  fec_a <- as.vector(t(pla) %*% wt_l)  # fecundity = weight
  expect_equal(mat_a * fec_a, mat_a * wt_a)
})

test_that("growth model matches MFCL outputs with MFCL parameters", {
  # MFCL baseline parameters (from bet.Rmd vignette)
  L1_mfcl <- 30.9192
  L2_mfcl <- 153.4431
  k_mfcl <- 0.09825
  CV1_mfcl <- 0.16101
  CV2_mfcl <- 0.1075122
  
  # Expected MFCL outputs from bet.Rmd
  mean_length_at_age_expected <- c(
    30.9192, 40.2986, 49.3364, 57.9282, 66.0202, 73.5897, 80.6340, 87.1632,
    93.1955, 98.7541, 103.8652, 108.5561, 112.8550, 116.7893, 120.3860, 123.6709,
    126.6685, 129.4018, 131.8926, 134.1612, 136.2262, 138.1052, 139.8142, 141.3680,
    142.7804, 144.0638, 145.2297, 146.2887, 147.2503, 148.1233, 148.9159, 149.6352,
    150.2880, 150.8804, 151.4179, 151.9055, 152.3478, 152.7491, 153.1130, 153.4431
  )
  
  sd_length_at_age_expected <- c(
    4.9783, 5.4564, 5.9606, 6.4830, 7.0169, 7.5559, 8.0948, 8.6284, 9.1527,
    9.6640, 10.1592, 10.6361, 11.0927, 11.5278, 11.9405, 12.3302, 12.6970, 13.0409,
    13.3625, 13.6622, 13.9409, 14.1994, 14.4387, 14.6597, 14.8636, 15.0513, 15.2239,
    15.3823, 15.5277, 15.6608, 15.7826, 15.8940, 15.9958, 16.0888, 16.1735, 16.2508,
    16.3213, 16.3854, 16.4439, 16.4970
  )
  
  # Compute model outputs with MFCL parameters
  mu_a <- get_growth(n_age, A1, A2, L1_mfcl, L2_mfcl, k_mfcl)
  sd_a <- get_sd_at_age(mu_a, L1_mfcl, L2_mfcl, CV1_mfcl, CV2_mfcl)
  
  # Verify mean length-at-age matches MFCL (within 0.1 cm tolerance)
  expect_equal(as.numeric(mu_a), mean_length_at_age_expected, tolerance = 0.1)
  
  # # Verify SD of length-at-age matches MFCL (within 0.05 cm tolerance)
  # expect_equal(sd_a, sd_length_at_age_expected, tolerance = 0.05)
  
  # Verify boundary conditions
  expect_equal(mu_a[A1], L1_mfcl, tolerance = 1e-10)
  expect_equal(mu_a[A2], L2_mfcl, tolerance = 1e-10)
  # expect_equal(sd_a[A1], L1_mfcl * CV1_mfcl, tolerance = 1e-10)
  # expect_equal(sd_a[A2], L2_mfcl * CV2_mfcl, tolerance = 1e-10)
})
