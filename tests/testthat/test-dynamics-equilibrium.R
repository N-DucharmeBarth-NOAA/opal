eq_data <- function(n_year, n_season, n_fishery, n_age,
                    catch_ysf = NULL, catch_units_f = rep(2L, n_fishery)) {
  if (is.null(catch_ysf)) catch_ysf <- array(0, c(n_year, n_season, n_fishery))
  list(first_yr = 1L, first_yr_catch = 1L, n_year = n_year,
       n_season = n_season, n_fishery = n_fishery, n_age = n_age,
       catch_obs_ysf = catch_ysf, catch_units_f = catch_units_f)
}

run_eq <- function(data, init, B0, h, M_a, sp_a, sel_fa,
                   weight_a = rep(1, data$n_age), rdev = 0, bias = 0,
                   sigma_r = 0.6) {
  nF <- data$n_fishery
  nY <- data$n_year
  nA <- data$n_age
  do_dynamics(
    data, list(rdev_y = rep(rdev, nY)),
    B0 = B0, R0 = init$R0, alpha = init$alpha, beta = init$beta, h = h,
    sigma_r = sigma_r, M_a = M_a, spawning_potential_a = sp_a,
    weight_fya = array(rep(weight_a, each = nF * nY), c(nF, nY, nA)),
    init_number_a = init$Ninit, init_number0_a = init$Ninit0,
    sel_fya = aperm(array(sel_fa, c(nF, nA, nY)), c(1, 3, 2)),
    bias_adj_y = rep(bias, nY)
  )
}

equilibrium_catch_ysf <- function(init, sel_a, F0, M_a, n_year, n_season) {
  u <- 1 - exp(-F0 / n_season)
  seasonal_number_a <- init$Ninit
  catch_s <- numeric(n_season)
  for (s in seq_len(n_season)) {
    catch_s[s] <- u * sum(seasonal_number_a * sel_a)
    seasonal_number_a <- seasonal_number_a * (1 - u * sel_a) *
      exp(-M_a / n_season)
  }
  array(rep(catch_s, each = n_year), c(n_year, n_season, 1))
}

test_that("unfished zero-catch dynamics remain at equilibrium", {
  n_age <- 12L
  n_year <- 40L
  M_a <- c(seq(0.6, 0.25, length.out = 6), rep(0.2, 6))
  sp_a <- c(0, 0, 0.1, 0.4, 0.8, rep(1, 7)) * seq_len(n_age)

  for (n_season in c(1L, 4L)) for (h in c(0.5, 0.95)) {
    lab <- sprintf("n_season=%d, h=%.2f", n_season, h)
    data <- eq_data(n_year, n_season, 2L, n_age)
    init <- get_initial_numbers(B0 = 1e6, h = h, M_a = M_a,
                                spawning_potential_a = sp_a)
    dyn <- run_eq(data, init, 1e6, h, M_a, sp_a,
                  sel_fa = matrix(1, 2, n_age))

    expect_equal(dyn$spawning_biomass_y, rep(1e6, n_year + 1),
                 tolerance = 1e-10, label = lab)
    expect_equal(dyn$spawning_biomass0_y, rep(1e6, n_year + 1),
                 tolerance = 1e-10, label = lab)
    expect_equal(dyn$dynamic_depletion_y, rep(1, n_year + 1),
                 tolerance = 1e-10, label = lab)
    for (y in c(2L, 20L, n_year + 1L)) {
      expect_equal(dyn$number_ysa[y, 1, ], init$Ninit,
                   tolerance = 1e-10, label = lab)
    }
  }
})

test_that("init_F_f equilibrium is stationary under matching catch (knife-edge selectivity)", {
  n_age <- 12L
  n_year <- 40L
  M_a <- rep(0.25, n_age)
  sp_a <- c(0, 0, 0.2, 0.6, rep(1, 8))
  sel_a <- as.numeric(seq_len(n_age) >= 4)
  F0 <- 0.3
  u <- 1 - exp(-F0)
  init <- get_initial_numbers(
    B0 = 1e6, h = 0.8, M_a = M_a, spawning_potential_a = sp_a,
    init_F_f = F0, sel_fa = matrix(sel_a, 1)
  )
  catch <- array(u * sum(init$Ninit * sel_a), c(n_year, 1, 1))
  data <- eq_data(n_year, 1L, 1L, n_age, catch_ysf = catch,
                  catch_units_f = 2L)
  dyn <- run_eq(data, init, 1e6, 0.8, M_a, sp_a,
                sel_fa = matrix(sel_a, 1))

  expect_equal(dyn$spawning_biomass_y,
               rep(dyn$spawning_biomass_y[1], n_year + 1),
               tolerance = 1e-8)
  expect_equal(dyn$number_ysa[n_year + 1, 1, ], init$Ninit,
               tolerance = 1e-8)
  expect_equal(dyn$spawning_biomass0_y, rep(1e6, n_year + 1),
               tolerance = 1e-10)
  expect_lt(dyn$lp_penalty, 1e-12)
})

test_that("init_F_f equilibrium with logistic selectivity is stationary", {
  n_age <- 12L
  n_year <- 40L
  M_a <- rep(0.25, n_age)
  sp_a <- c(0, 0, 0.2, 0.6, rep(1, 8))
  sel_a <- plogis(seq_len(n_age) - 4)
  F0 <- 0.3

  for (n_season in c(1L, 4L)) {
    init <- get_initial_numbers(
      B0 = 1e6, h = 0.8, M_a = M_a, spawning_potential_a = sp_a,
      init_F_f = F0, sel_fa = matrix(sel_a, 1), n_season = n_season
    )
    catch <- equilibrium_catch_ysf(init, sel_a, F0, M_a, n_year, n_season)
    data <- eq_data(n_year, n_season, 1L, n_age, catch_ysf = catch,
                    catch_units_f = 2L)
    dyn <- run_eq(data, init, 1e6, 0.8, M_a, sp_a,
                  sel_fa = matrix(sel_a, 1))

    expect_equal(dyn$spawning_biomass_y,
                 rep(dyn$spawning_biomass_y[1], n_year + 1),
                 tolerance = 1e-8, label = sprintf("n_season=%d", n_season))
    expect_equal(dyn$number_ysa[n_year + 1, 1, ], init$Ninit,
                 tolerance = 1e-8, label = sprintf("n_season=%d", n_season))
  }
})

test_that("the unfished trajectory ignores catch", {
  n_age <- 12L
  n_year <- 40L
  M_a <- rep(0.25, n_age)
  sp_a <- c(0, 0, 0.2, 0.6, rep(1, 8))
  sel_a <- as.numeric(seq_len(n_age) >= 4)
  F0 <- 0.3
  u <- 1 - exp(-F0)
  init <- get_initial_numbers(
    B0 = 1e6, h = 0.8, M_a = M_a, spawning_potential_a = sp_a,
    init_F_f = F0, sel_fa = matrix(sel_a, 1)
  )
  catch <- array(0.25 * u * sum(init$Ninit * sel_a), c(n_year, 1, 1))
  data <- eq_data(n_year, 1L, 1L, n_age, catch_ysf = catch,
                  catch_units_f = 2L)
  dyn <- run_eq(data, init, 1e6, 0.8, M_a, sp_a,
                sel_fa = matrix(sel_a, 1))

  expect_equal(dyn$spawning_biomass0_y, rep(1e6, n_year + 1),
               tolerance = 1e-10)
  expect_equal(dyn$dynamic_depletion_y,
               dyn$spawning_biomass_y / dyn$spawning_biomass0_y,
               tolerance = 1e-10)
})