# Tests for the joint age-length dynamics engine (do_dynamics_length)
library(RTMB)

# Shared synthetic setup (single fishery, one season)
n_age <- 25L; A1 <- 1L; A2 <- 25L; L1 <- 20; L2 <- 90; log_k <- log(0.18)
len_lower <- seq(5, by = 5, length.out = 20); len_upper <- len_lower + 5
len_mid <- (len_lower + len_upper) / 2
n_len <- length(len_lower); n_year <- 40L; n_season <- 1L; n_fishery <- 1L
mu_a <- get_growth(n_age, A1, A2, L1, L2, log_k)
B0 <- 5000; h <- 0.8; M_a <- rep(0.2, n_age)
mat_l <- 1 / (1 + exp(-0.3 * (len_mid - 45))); fec_l <- len_mid^3 / 1e5
sp_l <- mat_l * fec_l
wt_l <- 1e-5 * len_mid^3
sel_l <- 1 / (1 + exp(-0.4 * (len_mid - 40)))
catch <- c(rep(0, 5), seq(20, 120, length.out = 20), rep(120, 15))

mk_data <- function(catch_vec) {
  list(n_age = n_age, n_len = n_len, n_year = n_year, n_season = n_season,
       n_fishery = n_fishery, first_yr = 1, first_yr_catch = 1, catch_units_f = 1L,
       catch_obs_ysf = array(catch_vec, dim = c(n_year, 1, 1)))
}
params <- list(rdev_y = rep(0, n_year))

run_both <- function(log_CV1, log_CV2, catch_vec) {
  data <- mk_data(catch_vec)
  sd_a <- get_sd_at_age(mu_a, L1, L2, log_CV1, log_CV2)
  pla  <- get_pla(len_lower, len_upper, mu_a, sd_a)
  G    <- get_growth_matrix(len_lower, len_upper, len_mid, sd_a, L1, L2, log_k, A1, A2)
  rec  <- get_recruit_length_dist(len_lower, len_upper, mu_a[1], sd_a[1])
  sp_a <- as.vector(t(pla) %*% sp_l); wt_a <- as.vector(t(pla) %*% wt_l)
  sel_a <- as.vector(t(pla) %*% sel_l)
  wt_fya <- array(0, c(1, n_year, n_age)); sel_fya <- array(0, c(1, n_year, n_age))
  for (y in 1:n_year) { wt_fya[1, y, ] <- wt_a; sel_fya[1, y, ] <- sel_a }
  iA <- get_initial_numbers(B0, h, M_a, sp_a, sel_fa = matrix(sel_a, 1))
  dA <- do_dynamics(data, params, B0 = B0, R0 = iA$R0, alpha = iA$alpha, beta = iA$beta,
                    h = h, sigma_r = 0.6, M_a = M_a, spawning_potential_a = sp_a,
                    weight_fya = wt_fya, init_number_a = iA$Ninit,
                    init_number0_a = iA$Ninit0, sel_fya = sel_fya)
  wt_fyl <- array(0, c(1, n_year, n_len)); sel_fyl <- array(0, c(1, n_year, n_len))
  for (y in 1:n_year) { wt_fyl[1, y, ] <- wt_l; sel_fyl[1, y, ] <- sel_l }
  iL <- get_initial_numbers_length(B0, h, M_a, sp_l, pla, G, rec, init_F_f = 0,
                                   sel_fl = matrix(sel_l, 1))
  dL <- do_dynamics_length(data, params, B0 = B0, R0 = iL$R0, alpha = iL$alpha, beta = iL$beta,
                           h = h, sigma_r = 0.6, M_a = M_a, spawning_potential_l = sp_l,
                           weight_fyl = wt_fyl, recruit_dist_l = rec, G = G,
                           init_number_al = iL$Ninit_al, init_number0_al = iL$Ninit0_al,
                           sel_fyl = sel_fyl)
  list(A = dA, L = dL)
}

test_that("output arrays have the expected dimensions", {
  r <- run_both(log(0.12), log(0.08), catch)
  expect_equal(dim(r$L$number_ysal), c(n_year + 1, n_season, n_age, n_len))
  expect_equal(dim(r$L$number_ysa), c(n_year + 1, n_season, n_age))
  expect_equal(dim(r$L$NL_yl), c(n_year + 1, n_len))
  expect_equal(dim(r$L$catch_pred_fyl), c(n_fishery, n_year, n_len))
})

test_that("length marginal NL_yl is consistent with number_ysal", {
  r <- run_both(log(0.12), log(0.08), catch)
  for (y in c(1, 10, 25, n_year + 1)) {
    expect_equal(r$L$NL_yl[y, ], colSums(r$L$number_ysal[y, 1, , ]), tolerance = 1e-10)
  }
})

test_that("age marginal number_ysa is consistent with number_ysal", {
  r <- run_both(log(0.12), log(0.08), catch)
  for (y in c(1, 20, n_year + 1)) {
    expect_equal(r$L$number_ysa[y, 1, ], rowSums(r$L$number_ysal[y, 1, , ]), tolerance = 1e-10)
  }
})

test_that("predicted catch (weight) matches observed catch when feasible", {
  r <- run_both(log(0.12), log(0.08), catch)
  # catch_pred_ysf should track observed catch closely for the fished years
  yrs <- 10:n_year
  obs <- catch[yrs]
  pred <- r$L$catch_pred_ysf[yrs, 1, 1]
  expect_lt(max(abs(pred - obs) / obs), 0.05)
})

test_that("degenerate SD: length engine reduces to age-only dynamics", {
  r <- run_both(log(1e-3), log(1e-3), catch)
  # spawning biomass trajectory
  expect_lt(max(abs(r$L$spawning_biomass_y - r$A$spawning_biomass_y) /
                  r$A$spawning_biomass_y), 0.02)
  # terminal numbers-at-age (marginal)
  na <- r$A$number_ysa[n_year + 1, 1, ]
  nl <- r$L$number_ysa[n_year + 1, 1, ]
  expect_lt(max(abs(nl - na) / (na + 1e-8)), 0.03)
})

test_that("no-fishing degenerate run tracks age-only", {
  r <- run_both(log(1e-3), log(1e-3), rep(0, n_year))
  expect_lt(max(abs(r$L$spawning_biomass_y - r$A$spawning_biomass_y) /
                  r$A$spawning_biomass_y), 0.02)
})

test_that("M-at-length dynamics run and conserve sensibly", {
  data <- mk_data(catch)
  sd_a <- get_sd_at_age(mu_a, L1, L2, log(0.12), log(0.08))
  pla  <- get_pla(len_lower, len_upper, mu_a, sd_a)
  G    <- get_growth_matrix(len_lower, len_upper, len_mid, sd_a, L1, L2, log_k, A1, A2)
  rec  <- get_recruit_length_dist(len_lower, len_upper, mu_a[1], sd_a[1])
  M_l  <- rep(0.2, n_len)
  wt_fyl <- array(0, c(1, n_year, n_len)); sel_fyl <- array(0, c(1, n_year, n_len))
  for (y in 1:n_year) { wt_fyl[1, y, ] <- wt_l; sel_fyl[1, y, ] <- sel_l }
  iL <- get_initial_numbers_length(B0, h, M_l, sp_l, pla, G, rec, init_F_f = 0,
                                   sel_fl = matrix(sel_l, 1), M_basis = "length")
  dL <- do_dynamics_length(data, params, B0 = B0, R0 = iL$R0, alpha = iL$alpha, beta = iL$beta,
                           h = h, sigma_r = 0.6, M_a = M_l, spawning_potential_l = sp_l,
                           weight_fyl = wt_fyl, recruit_dist_l = rec, G = G,
                           init_number_al = iL$Ninit_al, init_number0_al = iL$Ninit0_al,
                           sel_fyl = sel_fyl, M_basis = "length")
  expect_equal(dL$spawning_biomass_y[1], B0, tolerance = 1e-6)
  expect_true(all(dL$spawning_biomass_y > 0))
})
