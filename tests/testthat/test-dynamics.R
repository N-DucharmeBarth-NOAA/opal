# Tests for dynamics module functions
library(RTMB)

# Test get_initial_numbers ----

test_that("get_initial_numbers with spawning_potential_a recovers B0 at equilibrium", {
  B0 <- 1e6
  h <- 0.8
  M_a <- rep(0.1, 10)
  maturity_a <- c(0, 0, 0, 0.5, 0.8, 1, 1, 1, 1, 1)
  weight_a <- seq(1, 10, length.out = 10)
  spawning_potential_a <- maturity_a * weight_a

  init <- get_initial_numbers(B0, h, M_a, spawning_potential_a)
  recovered_B0 <- sum(init$Ninit * spawning_potential_a)
  expect_equal(recovered_B0, B0, tolerance = 1e-6)
})

test_that("get_initial_numbers returns correct list elements", {
  B0 <- 1e6
  h <- 0.75
  M_a <- rep(0.2, 5)
  spawning_potential_a <- c(0, 0.2, 0.6, 1.0, 1.0)

  init <- get_initial_numbers(B0, h, M_a, spawning_potential_a)
  expect_named(init, c("Ninit", "R0", "alpha", "beta"))
  expect_equal(length(init$Ninit), length(M_a))
  expect_true(init$R0 > 0)
  expect_true(init$alpha > 0)
  expect_true(init$beta > 0)
})

test_that("get_initial_numbers BH parameters are internally consistent", {
  B0 <- 500000
  h <- 0.9
  M_a <- rep(0.15, 8)
  spawning_potential_a <- c(0, 0, 0.1, 0.4, 0.7, 1, 1, 1)

  init <- get_initial_numbers(B0, h, M_a, spawning_potential_a)
  # Check BH parameter derivation: alpha = 4*h*R0 / (5*h - 1)
  expect_equal(init$alpha, (4 * h * init$R0) / (5 * h - 1), tolerance = 1e-10)
  # Check beta: beta = B0 * (1 - h) / (5*h - 1)
  expect_equal(init$beta, (B0 * (1 - h)) / (5 * h - 1), tolerance = 1e-10)
})

test_that("do_dynamics returns expected dimensions with no fishing", {
  data <- list(
    first_yr = 1, first_yr_catch = 1,
    n_year = 2, n_season = 1, n_fishery = 1, n_age = 3,
    catch_obs_ysf = array(0, dim = c(2, 1, 1)),
    catch_units_f = 1
  )
  parameters <- list(rdev_y = rep(0, data$n_year))
  M_a <- rep(0.2, data$n_age)
  spawning_potential_a <- c(0, 1, 2)
  init <- get_initial_numbers(B0 = 1000, h = 0.75, M_a = M_a, spawning_potential_a = spawning_potential_a)

  dyn <- do_dynamics(
    data = data,
    parameters = parameters,
    B0 = 1000, R0 = init$R0, alpha = init$alpha, beta = init$beta, sigma_r = 0,
    M_a = M_a, spawning_potential_a = spawning_potential_a,
    weight_fya = array(1, dim = c(1, data$n_year, data$n_age)),
    init_number_a = init$Ninit,
    sel_fya = array(1, dim = c(1, data$n_year, data$n_age))
  )

  expect_equal(dim(dyn$number_ysa), c(data$n_year + 1, data$n_season, data$n_age))
  expect_equal(dim(dyn$catch_pred_fya), c(data$n_fishery, data$n_year, data$n_age))
  expect_equal(dyn$lp_penalty, 0, tolerance = 1e-12)
  expect_equal(sum(dyn$catch_pred_fya), 0, tolerance = 1e-12)
})

test_that("do_dynamics catch-at-age sums to observed catch in numbers", {
  data <- list(
    first_yr = 1, first_yr_catch = 1,
    n_year = 1, n_season = 1, n_fishery = 1, n_age = 3,
    catch_obs_ysf = array(10, dim = c(1, 1, 1)),
    catch_units_f = 2
  )
  parameters <- list(rdev_y = 0)
  M_a <- rep(0.2, data$n_age)
  spawning_potential_a <- c(0, 1, 2)
  init <- get_initial_numbers(B0 = 1000, h = 0.75, M_a = M_a, spawning_potential_a = spawning_potential_a)

  dyn <- do_dynamics(
    data = data,
    parameters = parameters,
    B0 = 1000, R0 = init$R0, alpha = init$alpha, beta = init$beta, sigma_r = 0,
    M_a = M_a, spawning_potential_a = spawning_potential_a,
    weight_fya = array(1, dim = c(1, data$n_year, data$n_age)),
    init_number_a = init$Ninit,
    sel_fya = array(1, dim = c(1, data$n_year, data$n_age))
  )

  expect_equal(sum(dyn$catch_pred_fya[1, 1, ]), data$catch_obs_ysf[1, 1, 1], tolerance = 1e-6)
})

test_that("get_harvest_rate reproduces observed catch in number and weight units", {
  number_ysa <- array(c(10, 20, 30), dim = c(2, 1, 3))
  sel_fya <- array(1, dim = c(1, 1, 3))

  data_numbers <- list(
    n_fishery = 1,
    n_age = 3,
    catch_obs_ysf = array(12, dim = c(1, 1, 1)),
    catch_units_f = 2
  )
  hr_numbers <- get_harvest_rate(data_numbers, 1, 1, number_ysa, sel_fya, array(1, dim = c(1, 1, 3)))
  expect_equal(sum(hr_numbers$h_rate_fa[1, ] * number_ysa[1, 1, ]), data_numbers$catch_obs_ysf[1, 1, 1], tolerance = 1e-7)

  weight_fya <- array(c(1, 2, 3), dim = c(1, 1, 3))
  data_weight <- list(
    n_fishery = 1,
    n_age = 3,
    catch_obs_ysf = array(14, dim = c(1, 1, 1)),
    catch_units_f = 1
  )
  hr_weight <- get_harvest_rate(data_weight, 1, 1, number_ysa, sel_fya, weight_fya)
  expect_equal(
    sum(hr_weight$h_rate_fa[1, ] * number_ysa[1, 1, ] * weight_fya[1, 1, ]),
    data_weight$catch_obs_ysf[1, 1, 1],
    tolerance = 1e-7
  )
})
