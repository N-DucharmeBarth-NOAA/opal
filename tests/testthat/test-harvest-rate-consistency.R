test_that("get_harvest_rate reproduces do_dynamics in a mixed-unit case", {
  n_year <- 5L
  n_season <- 2L
  n_fishery <- 3L
  n_age <- 8L
  ages <- seq_len(n_age)
  catch <- array(0, c(n_year, n_season, n_fishery))
  catch[, , 1] <- 40
  catch[, , 2] <- 25
  catch[, 2, 3] <- 10
  data <- list(
    first_yr = 1L, first_yr_catch = 3L, n_year = n_year,
    n_season = n_season, n_fishery = n_fishery, n_age = n_age,
    catch_obs_ysf = catch, catch_units_f = c(1L, 2L, 1L)
  )

  sel_fya <- array(0, c(n_fishery, n_year, n_age))
  for (y in seq_len(n_year)) {
    sel_fya[1, y, ] <- plogis(ages - 3 - 0.1 * y)
    sel_fya[2, y, ] <- exp(-0.5 * ((ages - 4) / 1.5)^2)
    sel_fya[3, y, ] <- plogis(ages - 5)
  }
  weight_fya <- array(rep(0.2 * ages^1.5, each = n_fishery * n_year),
                      c(n_fishery, n_year, n_age))
  M_a <- rep(0.3, n_age)
  sp_a <- c(0, 0, 0.5, rep(1, 5)) * 0.2 * ages^1.5
  init <- get_initial_numbers(B0 = 1e4, h = 0.8, M_a = M_a,
                              spawning_potential_a = sp_a)

  dyn <- do_dynamics(
    data, list(rdev_y = c(0.2, -0.1, 0, 0.3, -0.2)), B0 = 1e4,
    R0 = init$R0, alpha = init$alpha, beta = init$beta, sigma_r = 0.6,
    M_a = M_a, spawning_potential_a = sp_a, weight_fya = weight_fya,
    init_number_a = init$Ninit, init_number0_a = init$Ninit0,
    sel_fya = sel_fya, bias_adj_y = rep(1, n_year)
  )

  fy <- data$first_yr_catch - data$first_yr + 1L
  penalty <- 0
  for (y in seq_len(n_year)) for (s in seq_len(n_season)) {
    lab <- sprintf("y=%d s=%d", y, s)
    if (y < fy) {
      expect_true(all(dyn$hrate_ysfa[y, s, , ] == 0), label = lab)
      next
    }
    hr <- get_harvest_rate(data, y, s, dyn$number_ysa, sel_fya, weight_fya)
    expect_equal(hr$h_rate_fa, dyn$hrate_ysfa[y, s, , ], tolerance = 1e-12,
                 label = lab)
    expect_equal(hr$h_rate_a, dyn$hrate_ysa[y, s, ], tolerance = 1e-12,
                 label = lab)
    penalty <- penalty + hr$penalty
  }
  expect_equal(dyn$lp_penalty, penalty, tolerance = 1e-12)
  expect_equal(dyn$catch_pred_ysf[fy:n_year, , ], catch[fy:n_year, , ],
               tolerance = 1e-8)
})

test_that("get_harvest_rate and do_dynamics agree when posfun is active", {
  n_age <- 3L
  data <- list(
    first_yr = 1L, first_yr_catch = 1L, n_year = 1L, n_season = 1L,
    n_fishery = 1L, n_age = n_age,
    catch_obs_ysf = array(1e6, c(1, 1, 1)), catch_units_f = 2L
  )
  parameters <- list(rdev_y = 0)
  M_a <- rep(0.2, n_age)
  spawning_potential_a <- c(0, 1, 2)
  init <- get_initial_numbers(1000, 0.75, M_a, spawning_potential_a)
  sel_fya <- array(1, c(1, 1, n_age))
  weight_fya <- array(1, c(1, 1, n_age))

  dyn <- do_dynamics(
    data, parameters, B0 = 1000, R0 = init$R0, alpha = init$alpha,
    beta = init$beta, sigma_r = 0, M_a = M_a,
    spawning_potential_a = spawning_potential_a, weight_fya = weight_fya,
    init_number_a = init$Ninit, init_number0_a = init$Ninit0,
    sel_fya = sel_fya
  )
  hr <- get_harvest_rate(data, 1, 1, dyn$number_ysa, sel_fya, weight_fya)

  expect_equal(as.numeric(hr$h_rate_fa),
               as.numeric(dyn$hrate_ysfa[1, 1, , ]), tolerance = 1e-12)
  expect_equal(hr$penalty, dyn$lp_penalty, tolerance = 1e-12)
  expect_gt(hr$penalty, 0)
})

test_that("get_harvest_rate reproduces opal_model harvest rates", {
  inputs <- opaka_inputs()
  data <- inputs$data
  obj <- opaka_obj(inputs)
  report <- obj$report(obj$par)
  fy <- data$first_yr_catch - data$first_yr + 1L
  penalty <- 0

  for (y in fy:data$n_year) for (s in seq_len(data$n_season)) {
    hr <- get_harvest_rate(data, y, s, report$number_ysa,
                           report$sel_fya, report$weight_fya_mod)
    expect_equal(hr$h_rate_fa, report$hrate_ysfa[y, s, , ], tolerance = 1e-12)
    penalty <- penalty + hr$penalty
  }
  expect_equal(report$lp_penalty, penalty, tolerance = 1e-12)
})

test_that("harvest-rate functions reject invalid catch unit codes", {
  data <- list(
    first_yr = 1L, first_yr_catch = 1L, n_year = 1L, n_season = 1L,
    n_fishery = 3L, n_age = 2L,
    catch_obs_ysf = array(0, c(1, 1, 3)), catch_units_f = c(1L, 3L, 1L)
  )
  expect_error(
    get_harvest_rate(data, 1, 1, array(1, c(2, 1, 2)),
                     array(1, c(3, 1, 2)), array(1, c(3, 1, 2))),
    "catch_units_f"
  )

  init <- get_initial_numbers(1000, 0.75, rep(0.2, 2), c(0, 1))
  expect_error(
    do_dynamics(
      data, list(rdev_y = 0), B0 = 1000, R0 = init$R0,
      alpha = init$alpha, beta = init$beta, sigma_r = 0,
      M_a = rep(0.2, 2), spawning_potential_a = c(0, 1),
      weight_fya = array(1, c(3, 1, 2)), init_number_a = init$Ninit,
      init_number0_a = init$Ninit0, sel_fya = array(1, c(3, 1, 2))
    ),
    "catch_units_f"
  )
})