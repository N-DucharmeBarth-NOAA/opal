test_that("Dynamic B0 equals fished B0 when catch is zero", {
  mock_data <- list(
    first_yr = 1, first_yr_catch = 1, n_year = 5, n_season = 1,
    n_fishery = 1, n_age = 5, catch_units_f = c(1)
  )
  mock_data$catch_obs_ysf <- array(0, dim = c(5, 1, 1))

  mock_params <- list(rdev_y = rep(0.1, 5))
  bias_adj_y <- rep(1, 5)

  M_a <- rep(0.2, 5)
  spawning_potential_a <- c(0, 0, 10, 20, 30)
  weight_fya <- array(1, dim = c(1, 5, 5))
  sel_fya <- array(1, dim = c(1, 5, 5))

  init <- get_initial_numbers(B0 = 1000, h = 0.8, M_a = M_a,
                              spawning_potential_a = spawning_potential_a,
                              init_F_f = c(0), sel_fa = sel_fya[, 1, ],
                              init_rdev_a = rep(0, 5), sigma_r = 0.6,
                              init_bias_adj_a = rep(1, 5))

  dyn <- do_dynamics(mock_data, mock_params, B0 = 1000, R0 = init$R0,
                     alpha = init$alpha, beta = init$beta, h = 0.8, sigma_r = 0.6,
                     M_a = M_a, spawning_potential_a = spawning_potential_a,
                     weight_fya = weight_fya, init_number_a = init$Ninit,
                     init_number0_a = init$Ninit0, sel_fya = sel_fya,
                     bias_adj_y = bias_adj_y)

  expect_equal(dyn$spawning_biomass0_y, dyn$spawning_biomass_y)
  expect_equal(dyn$static_depletion_y, dyn$spawning_biomass_y / 1000)
  expect_equal(dyn$dynamic_depletion_y, rep(1, 6))
})
