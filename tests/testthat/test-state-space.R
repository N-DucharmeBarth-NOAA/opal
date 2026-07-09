library(RTMB)

state_space_fixture <- function() {
  data <- list(
    first_yr = 1L,
    first_yr_catch = 1L,
    n_year = 3L,
    n_season = 1L,
    n_fishery = 1L,
    n_age = 3L,
    catch_obs_ysf = array(0, dim = c(3L, 1L, 1L)),
    catch_units_f = 2L
  )
  M_a <- c(0.3, 0.2, 0.2)
  spawning_potential_a <- c(0, 0.5, 1)
  init <- get_initial_numbers(
    B0 = 1000,
    h = 0.75,
    M_a = M_a,
    spawning_potential_a = spawning_potential_a
  )
  list(
    data = data,
    parameters = list(rdev_y = rep(0, data$n_year)),
    M_a = M_a,
    spawning_potential_a = spawning_potential_a,
    init = init,
    weight_fya = array(1, dim = c(1L, data$n_year, data$n_age)),
    sel_fya = array(1, dim = c(1L, data$n_year, data$n_age))
  )
}

run_state_space_fixture <- function(x, data = x$data, parameters = x$parameters) {
  do_dynamics(
    data = data,
    parameters = parameters,
    B0 = 1000,
    R0 = x$init$R0,
    alpha = x$init$alpha,
    beta = x$init$beta,
    h = 0.75,
    sigma_r = 0.2,
    M_a = x$M_a,
    spawning_potential_a = x$spawning_potential_a,
    weight_fya = x$weight_fya,
    init_number_a = x$init$Ninit,
    init_number0_a = x$init$Ninit0,
    sel_fya = x$sel_fya
  )
}

test_that("state process density is zero-residual at its transition mean", {
  pred <- c(100, 80, 50)
  sigma <- c(0.1, 0.2, 0.3)
  out <- get_state_process_nll(
    log_number_state_a = log(pred),
    number_pred_a = pred,
    log_sigma_state = log(sigma),
    bias_correct = FALSE
  )

  expect_equal(out$residual_a, rep(0, 3), tolerance = 1e-12)
  expect_equal(
    out$nll,
    -sum(dnorm(log(pred), log(pred), sigma, log = TRUE)),
    tolerance = 1e-12
  )
})

test_that("state-space parameters initialise from a deterministic trajectory", {
  x <- state_space_fixture()
  deterministic <- run_state_space_fixture(x)
  parameters <- initialize_state_space_parameters(
    x$parameters,
    deterministic$number_ysa,
    process_sigma = 0.15
  )

  expect_equal(dim(parameters$log_number_state_ya), c(x$data$n_year, x$data$n_age))
  expect_equal(
    exp(parameters$log_number_state_ya),
    deterministic$number_ysa[2:(x$data$n_year + 1L), 1L, ]
  )
  expect_equal(parameters$log_sigma_state, log(0.15))
})

test_that("state-space switch preserves deterministic dynamics at matching states", {
  x <- state_space_fixture()
  deterministic <- run_state_space_fixture(x)
  parameters <- initialize_state_space_parameters(
    x$parameters,
    deterministic$number_ysa,
    process_sigma = 0.1
  )
  data <- x$data
  data$state_space_switch <- 1L
  data$state_process_bias_correct <- 0L

  state_space <- run_state_space_fixture(x, data = data, parameters = parameters)

  expect_equal(state_space$number_ysa, deterministic$number_ysa, tolerance = 1e-10)
  expect_equal(state_space$number_pred_ya, exp(parameters$log_number_state_ya),
               tolerance = 1e-10)
  expect_equal(state_space$state_residual_ya, matrix(0, data$n_year, data$n_age),
               tolerance = 1e-10)
  expect_equal(
    state_space$lp_state,
    -data$n_year * data$n_age * dnorm(0, 0, 0.1, log = TRUE),
    tolerance = 1e-8
  )
})

test_that("latent states feed forward into later transition predictions", {
  x <- state_space_fixture()
  deterministic <- run_state_space_fixture(x)
  parameters <- initialize_state_space_parameters(x$parameters, deterministic$number_ysa)
  parameters$log_number_state_ya[1L, 2L] <-
    parameters$log_number_state_ya[1L, 2L] + log(1.2)
  data <- x$data
  data$state_space_switch <- 1L
  data$state_process_bias_correct <- 0L

  state_space <- run_state_space_fixture(x, data = data, parameters = parameters)

  expect_equal(
    state_space$number_ysa[2L, 1L, 2L],
    deterministic$number_ysa[2L, 1L, 2L] * 1.2,
    tolerance = 1e-10
  )
  expect_false(isTRUE(all.equal(
    state_space$number_pred_ya[2L, ],
    deterministic$number_ysa[3L, 1L, ]
  )))
})

test_that("latent log states can be registered as RTMB random effects", {
  objective <- function(parameters, data) {
    get_state_process_nll(
      log_number_state_a = parameters$log_number_state_a,
      number_pred_a = data$number_pred_a,
      log_sigma_state = parameters$log_sigma_state,
      bias_correct = FALSE
    )$nll
  }
  parameters <- list(
    log_number_state_a = log(c(100, 80, 50)),
    log_sigma_state = log(0.2)
  )
  obj <- MakeADFun(
    func = cmb(objective, list(number_pred_a = c(100, 80, 50))),
    parameters = parameters,
    random = "log_number_state_a",
    silent = TRUE
  )

  expect_length(obj$env$random, 3L)
  expect_true(is.finite(obj$fn()))
})
