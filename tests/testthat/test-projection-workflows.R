test_that("MVN projections join the historical trajectory for the same draw", {
  x <- opal_fit(small_opal_object())
  obj <- opal_rtmb(x, fresh = TRUE)
  n <- 2L
  args <- list(n_proj = 2L, n_iter = n, rdev_y = matrix(0, n, 2),
               sel_fya = project_selectivity(x, n_proj = 2, n_iter = n),
               catch_ysf = array(0, c(2, 1, 2)))
  # This fixture estimates only log_B0: its scalar Gaussian draws are an
  # independent reference for the full historical/projection bridge.
  set.seed(72)
  draws <- as.numeric(obj$par) + rnorm(n) / sqrt(as.numeric(obj$he()))
  expected <- lapply(draws, function(value) obj$report(value))
  set.seed(72)
  result <- do.call(project_dynamics, c(list(data = x, uncertainty = "mvn",
                                            return_hist = TRUE), args))
  expect_equal(dim(result$hist_sbio), c(n, x$data$n_year + 1L))
  for (i in seq_len(n)) {
    expect_equal(result$hist_sbio[i, ], as.numeric(expected[[i]]$spawning_biomass_y))
    expect_equal(result$dyn[[i]]$spawning_biomass_y[1],
                  tail(result$hist_sbio[i, ], 1))
    expect_true(all(is.finite(result$dyn[[i]]$number_ysa)))
    expect_true(all(result$dyn[[i]]$number_ysa >= 0))
  }
  set.seed(72)
  without_history <- do.call(project_dynamics, c(list(data = x), args))
  expect_equal(without_history, result$dyn)
  obj <- opal_rtmb(x, fresh = TRUE)
  # The finite-difference fallback has its own approximate curvature.
  # Verify its draws against that curvature, rather than equating it to AD.
  numerical_hessian <- stats::optimHess(obj$par, obj$fn, obj$gr)
  set.seed(72)
  numerical_draws <- as.numeric(obj$par) + rnorm(n) / sqrt(as.numeric(numerical_hessian))
  numerical_bridge <- vapply(numerical_draws, function(value) {
    tail(as.numeric(obj$report(value)$spawning_biomass_y), 1)
  }, numeric(1))
  obj <- opal_rtmb(x, fresh = TRUE)
  obj["he"] <- list(NULL)
  set.seed(72)
  numerical <- do.call(project_dynamics, c(list(data = x$data, object = obj), args))
  expect_equal(vapply(numerical, function(draw) draw$spawning_biomass_y[1], numeric(1)),
                numerical_bridge)
  obj$he <- function() matrix(0, 1, 1)
  expect_error(do.call(project_dynamics, c(list(data = x$data, object = obj), args)),
                "Hessian is singular")
})

test_that("projection helpers reject mismatched inputs and missing fit state", {
  x <- opal_fit(small_opal_object())
  args <- list(n_proj = 2L, n_iter = 2L, rdev_y = matrix(0, 2, 2),
               sel_fya = project_selectivity(x, n_proj = 2, n_iter = 2),
               catch_ysf = array(0, c(2, 1, 2)))
  run <- function(changes) do.call(project_dynamics, c(list(data = x),
                                                         utils::modifyList(args, changes)))
  expect_error(run(list(rdev_y = c(0, 0))), "2-D matrix")
  expect_error(run(list(rdev_y = matrix(0, 1, 2))), "number of rows")
  expect_error(run(list(sel_fya = matrix(1, 2, 2))), "4-D array")
  expect_error(run(list(sel_fya = array(1, c(1, 2, 2, 5)))), "first dimension")
  expect_error(run(list(catch_ysf = array(0, c(3, 1, 2)))), "3-D array")
  expect_error(project_rec_devs(small_opal_object()), "requires a fitted model")
  expect_error(project_rec_devs(x, obj = opal_rtmb(x)), "stored in the opal_obj")
  expect_error(project_rec_devs(x$data, obj = opal_rtmb(x), uncertainty = "fit"),
                "Use uncertainty with an opal_obj")
  expect_error(opal_project(x, name = ""), "one result name")
  expect_error(opal_project(x, data = x$data), "inputs come from x")
  sampled <- opal_attach_mcmc(x, mock_opal_sampler(opal_rtmb(x), num_samples = 2),
                              check = FALSE)
  expect_error(project_rec_devs(sampled, uncertainty = "mcmc", n_iter = 5),
                "available posterior draws")
  args$n_iter <- 5L
  args$rdev_y <- matrix(0, 5, 2)
  args$sel_fya <- array(1, c(5, 2, 2, 5))
  expect_error(do.call(project_dynamics,
                       c(list(data = sampled, uncertainty = "mcmc"), args)),
                "available MCMC draws")
})

test_that("recruitment forecasts preserve horizons and selected historical windows", {
  data <- make_opal_fit_data()
  data$n_year <- data$last_yr <- 20L
  data$catch_obs_ysf <- array(100, c(20, 1, 2))
  parameters <- make_opal_fit_parameters(data)
  parameters$log_h <- log(0.75)
  parameters$rdev_y <- sin(seq_len(20) * 1.3) / 5 + seq_len(20) / 100
  x <- opal_obj(data, parameters, make_opal_fit_map(parameters))
  obj <- opal_rtmb(x)
  x <- opal_attach_fit(x, list(par = obj$par, objective = obj$fn(obj$par),
                              convergence = 0L), check = FALSE)
  for (use_arima in c(TRUE, FALSE)) {
    set.seed(123)
    result <- project_rec_devs(x, first_yr = 5, last_yr = 18,
      n_proj = 3, n_iter = 2, arima = use_arima, max.p = 1, max.d = 0, max.q = 0)
    expect_equal(dim(result$rdev_y), c(2L, 3L))
    expect_equal(colnames(result$rdev_y), as.character(21:23))
    expect_true(all(is.finite(result$rdev_y)))
    expect_true(all(result$arima_pars[, 1] %in% 0:1))
    # Reproduce the first forecast independently from the selected history.
    history <- parameters$rdev_y[5:18]
    set.seed(123)
    if (use_arima) {
      model <- forecast::auto.arima(history, max.p = 1, max.d = 0, max.q = 0,
        approximation = FALSE, stepwise = FALSE, ic = "bic")
      expected <- stats::simulate(model, nsim = 3, future = TRUE, bootstrap = TRUE)
    } else {
      model <- stats::ar(history, order.max = 1)
      expected <- stats::arima.sim(list(ar = model$ar), n = 3, sd = sqrt(model$var.pred))
    }
    expect_equal(as.numeric(result$rdev_y[1, ]), as.numeric(expected))
  }
  constant <- opal_fit(small_opal_object())
  result <- project_rec_devs(constant, arima = FALSE, n_iter = 2, n_proj = 3)
  expect_equal(unname(result$rdev_y), matrix(0, 2, 3))
  expect_equal(result$arima_pars, matrix(0, 2, 3))
})

test_that("selectivity forecasts respect variable, constant, and inactive series", {
  data <- list(first_yr = 2000L, last_yr = 2003L, age_a = 1:2,
               removal_switch_f = c(0, 1))
  sel <- array(0.5, c(2, 4, 2))
  sel[1, , 1] <- c(0.1, 0.2, 0.4, 0.8)
  object <- list(report = function() list(sel_fya = sel))
  set.seed(44)
  result <- project_selectivity(data, object, first_yr = 2001,
                                  n_iter = 2, n_proj = 3)
  expect_equal(dim(result), c(2L, 2L, 3L, 2L))
  expect_equal(dimnames(result)$year, as.character(2004:2006))
  expect_equal(unname(result[, 1, , 2]), matrix(0.5, 2, 3))
  expect_true(all(result[, 2, , ] == 0))
  set.seed(44)
  expected <- exp(rnorm(3, mean(log(sel[1, 2:4, 1])), sd(log(sel[1, 2:4, 1]))))
  expect_equal(as.numeric(result[1, 1, , 1]), expected)
})
