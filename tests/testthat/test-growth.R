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
k  <- log(0.2)
CV1 <- log(0.15)
CV2 <- log(0.08)

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

test_that("default and empty multipliers preserve the original growth exactly", {
  ages <- seq_len(n_age)
  original <- L1 + (L2 - L1) * (1 - exp(-exp(k) * (ages - A1))) /
    (1 - exp(-exp(k) * (A2 - A1)))
  expect_identical(get_growth(n_age, A1, A2, L1, L2, k), original)
  expect_identical(get_growth(n_age, A1, A2, L1, L2, k,
                              log_k_mult = NULL), original)
  expect_identical(get_growth(n_age, A1, A2, L1, L2, k,
                              log_k_mult = numeric(0)), original)
  for (mult_ages in list(1L, c(2L, 10L), c(1L, 20L, 40L))) {
    expect_equal(get_growth(n_age, A1, A2, L1, L2, k,
                            log_k_mult = rep(0, length(mult_ages)),
                            k_mult_ages = mult_ages),
                 original, tolerance = 1e-10)
  }
})

# Golden lengths: Spatial Stock Assessment Simulation Experiment, IO YFT,
# Table A1 (Barcelo, Berger, Goethel, Hoyle et al.), GPL-3:
# https://github.com/aaronmberger-nwfsc/Spatial-Assessment-Modeling-Workshop/blob/2174509/vignettes/OM_description_YFT.Rmd
# https://aaronmberger-nwfsc.github.io/Spatial-Assessment-Modeling-Workshop/articles/OM_description_YFT.html
# Identical to SS3 v3.24AB endgrowth Len_Beg, replicate 4:
# N-DucharmeBarth/wcpfc-bet-single-region, commit e7d8c8e,
# model-files/ss3/xx-io-yft-sim-iattc-sample-tree-4/Report.sso.
# Growth: Fu et al. (2018), after Fonteneau and Gascuel (2008).
test_that("age-specific K reproduces SS3 GrowthModel 3 for IO yellowfin", {
  golden <- c(22, 35.2865, 41.384, 45.7348, 49.904, 53.8991, 57.7273,
              64.2198, 74.721, 85.5489, 96.3811, 105.259, 112.534,
              118.496, 123.383, 127.387, 130.669, 133.359, 135.563,
              137.369, 138.85, 140.063, 141.058, 141.872, 142.54,
              143.088, 143.536, 144.168)
  args <- list(n_age = 28L, A1 = 1L, A2 = 28L, L1 = 22, L2 = 145,
               log_k = log(0.455 / 4),
               log_k_mult = log(c(0.5, 0.75, 1, 1, 1, 1.8, 1.8, 1.2, 1.2)),
               k_mult_ages = 2:10)
  adjusted <- do.call(get_growth, c(args, list(plus_group_adjust = TRUE)))
  unadjusted <- do.call(get_growth, args)
  # The source values are printed to six significant figures.
  expect_equal(adjusted, golden, tolerance = 1e-3)
  expect_lt(max(abs(adjusted - golden)), 1e-3)
  expect_equal(unadjusted[1:27], golden[1:27], tolerance = 1e-3)
  expect_lt(max(abs(unadjusted[1:27] - golden[1:27])), 1e-3)
  expect_equal(unadjusted[28], 143.904, tolerance = 1e-3)
  expect_lt(abs(unadjusted[28] - 143.904), 1e-3)
  expect_identical(adjusted[1:27], unadjusted[1:27])
})

test_that("K multipliers are cumulative and the preceding age governs growth", {
  # K = (0.2, 0.1, 0.3, 0.3); cumulative elapsed K = (0, 0.2, 0.3, 0.6).
  Linf <- 10 + 40 / (1 - exp(-0.3 * 3))
  expected <- Linf + (10 - Linf) * exp(-c(0, 0.2, 0.3, 0.6))
  mu <- get_growth(4L, 1L, 4L, 10, 50, log(0.2),
                   log_k_mult = log(c(0.5, 3)), k_mult_ages = 2:3)
  expect_equal(mu, expected, tolerance = 1e-12)
  expect_gt(abs(mu[4] - 50), 1)
})

test_that("ages below A1 retain base-K exponential growth", {
  Linf <- 20 + 40 / (1 - exp(-0.1 * 3))
  mu <- get_growth(6L, 2L, 5L, 20, 60, log(0.2), min_age = 0L,
                   log_k_mult = log(0.5), k_mult_ages = 3L)
  expect_equal(mu[1:3], Linf + (20 - Linf) * exp(-0.2 * (0:2 - 2)),
               tolerance = 1e-12)
  expect_length(get_growth(1L, 1L, 2L, 20, 60, log(0.2),
                           log_k_mult = 0, k_mult_ages = 1L), 1L)
})

test_that("plus-group adjustment also works without multipliers", {
  plain <- get_growth(n_age, A1, A2, L1, L2, k)
  Linf <- L1 + (L2 - L1) / (1 - exp(-exp(k) * (A2 - A1)))
  j <- 0:n_age
  w <- exp(-0.2)^j
  adjusted <- get_growth(n_age, A1, A2, L1, L2, k, plus_group_adjust = TRUE)
  expect_identical(adjusted[-n_age], plain[-n_age])
  expect_equal(adjusted[n_age],
               sum(w * (plain[n_age] + j / n_age * (Linf - plain[n_age]))) / sum(w))
})

test_that("age-specific growth derivatives match central finite differences", {
  par <- c(22, 145, log(0.455 / 4), log(c(0.5, 0.75, 1.8)))
  for (adjust in c(FALSE, TRUE)) {
    f <- function(p) sum(get_growth(
      28L, 1L, 28L, p[1], p[2], p[3],
      log_k_mult = p[4:6], k_mult_ages = c(2L, 3L, 7L),
      plus_group_adjust = adjust
    )^2)
    tape <- MakeTape(f, par)
    expect_equal(as.numeric(tape(par)), f(par), tolerance = 1e-12)
    fd <- vapply(seq_along(par), function(i) {
      step <- rep(0, length(par))
      step[i] <- 1e-5
      (f(par + step) - f(par - step)) / (2e-5)
    }, numeric(1))
    expect_lt(max(abs(as.numeric(tape$jacobian(par)) - fd) / pmax(1, abs(fd))), 1e-6)
  }
})

test_that("invalid multiplier inputs fail clearly", {
  growth <- function(...) get_growth(5L, 1L, 5L, 10, 50, log(0.2), ...)
  expect_error(growth(log_k_mult = c(0, 0), k_mult_ages = 2L), "same length")
  expect_error(growth(k_mult_ages = 2L), "same length")
  for (ages in list(0L, 6L, 2.5, NA_real_, Inf, "two")) {
    expect_error(growth(log_k_mult = 0, k_mult_ages = ages), "integer ages.*range")
  }
  for (ages in list(c(3L, 2L), c(2L, 2L))) {
    expect_error(growth(log_k_mult = c(0, 0), k_mult_ages = ages), "strictly increasing")
  }
  expect_error(get_growth(5L, 0L, 5L, 10, 50, log(0.2),
                          log_k_mult = 0, k_mult_ages = 2L), "A1.*modeled age")
  expect_error(growth(plus_group_adjust = NA), "logical flag")
  expect_error(growth(plus_group_adjust = 2), "logical flag")
  expect_error(get_growth(1L, 0L, 1L, 10, 50, log(0.2), min_age = 0L,
                          plus_group_adjust = TRUE), "positive integer")
})

test_that("growth parameters have defaults, fixed maps, priors and unbounded log rates", {
  inputs <- opaka_inputs()
  data <- inputs$data
  data$k_mult_ages <- c(2L, 4L)
  parameters <- get_parameters(data, model = "opakapaka")
  expect_identical(parameters$log_k_mult, c(0, 0))
  expect_false("log_k_mult" %in% names(get_parameters(model = "opakapaka")))
  expect_true(all(is.na(get_map(parameters)$log_k_mult)))
  parameters$log_k_mult <- log(c(0.5, 1.8))
  prior <- get_priors(parameters, data)$log_k_mult
  expect_identical(prior$par1, parameters$log_k_mult)
  expect_identical(prior$par2, 0.3)
  data$prior_log_k_mult_mean <- c(-0.2, 0.2)
  expect_identical(get_priors(parameters, data)$log_k_mult$par1,
                   data$prior_log_k_mult_mean)
  obj <- list(par = c(log_k = -1, log_k_mult = -0.5))
  bounds <- get_bounds(obj, parameters)
  expect_identical(bounds$lower, c(-Inf, -Inf))
  expect_identical(bounds$upper, c(Inf, Inf))
})

test_that("opal_model forwards optional growth inputs and preserves the default", {
  inputs <- opaka_inputs()
  data <- inputs$data
  data$A1 <- 1L
  parameters <- inputs$parameters
  data$k_mult_ages <- c(2L, 4L)
  data$growth_plus_group_adjust <- 1L
  parameters$log_k_mult <- log(c(0.5, 1.8))
  map <- c(inputs$map, get_map(parameters)["log_k_mult"])
  obj <- opaka_obj(data = data, parameters = parameters, map = map)
  direct_model <- opal_model
  environment(direct_model) <- list2env(list(
    get_growth = function(n_age, A1, A2, L1, L2, log_k, ...) {
      get_growth(n_age, A1, A2, L1, L2, log_k,
                  log_k_mult = log(c(0.5, 1.8)), k_mult_ages = c(2L, 4L),
                  plus_group_adjust = TRUE)
    }
  ), parent = environment(opal_model))
  direct <- MakeADFun(function(p) direct_model(p, data), parameters,
                       map = map, silent = TRUE)
  expect_equal(obj$fn(obj$par), direct$fn(direct$par), tolerance = 1e-12)
  expect_true(all(is.finite(obj$gr(obj$par))))
  baseline <- opaka_obj()
  plain_data <- inputs$data
  plain_data$growth_plus_group_adjust <- 0L
  plain <- opaka_obj(data = plain_data)
  expect_identical(plain$fn(plain$par), baseline$fn(baseline$par))
  expect_identical(plain$gr(plain$par), baseline$gr(baseline$par))
  data$k_mult_ages <- 2L
  expect_error(opal_obj(data, parameters, map), "same length")
  expect_error(opaka_obj(data = data, parameters = parameters, map = map), "same length")
  data$k_mult_ages <- c(4L, 2L)
  expect_error(opal_obj(data, parameters, map), "strictly increasing")
  data$k_mult_ages <- c(2L, data$n_age + 1L)
  expect_error(opal_obj(data, parameters, map), "integer ages.*range")
})

# Test get_sd_at_age ----

test_that("get_sd_at_age returns correct SD at reference ages", {
  mu_a <- get_growth(n_age, A1, A2, L1, L2, k)
  sd_a <- get_sd_at_age(mu_a, L1, L2, CV1, CV2)
  # At age A1: mu = L1, so SD = L1 * CV1
  expect_equal(sd_a[A1], L1 * exp(CV1), tolerance = 1e-10)
  # At age A2: mu = L2, so SD = L2 * CV2
  expect_equal(sd_a[A2], L2 * exp(CV2), tolerance = 1e-10)
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
  mu_a2 <- get_growth(50L, 5L, 45L, 40.0, 200.0, log(0.15))
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
  # MFCL baseline parameters (from bet.qmd vignette)
  L1_mfcl <- 30.9192
  L2_mfcl <- 153.4431
  k_mfcl <- log(0.09825)
  CV1_mfcl <- log(0.16101)
  CV2_mfcl <- log(0.1075122)
  
  # Expected MFCL outputs from bet.qmd
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
