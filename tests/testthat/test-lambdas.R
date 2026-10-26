test_that("likelihood lambdas default to one and weighted terms sum to the objective", {
  data <- synth_full_data(wf_switch = 1L, lf_switch = 1L)
  parameters <- synth_parameters(data)
  map <- synth_map(parameters)
  map$log_cpue_q <- NULL
  map$log_cpue_tau <- NULL
  map$log_lf_tau <- NULL
  map$log_wf_tau <- NULL
  baseline <- synth_obj(data, parameters, map)
  data$lambdas <- list(cpue = 1, lf = 1, wf = 1, rec = 1,
                       init_rec = 1, prior = 1, penalty = 1)
  defaults <- synth_obj(data, parameters, map)
  expect_identical(defaults$fn(), baseline$fn())
  expect_equal(defaults$report()$lp_lf, baseline$report()$lp_lf)

  data$lambdas <- list(cpue = 0.5, lf = c(2, 0), wf = 3, rec = 2)
  weighted <- synth_obj(data, parameters, map)
  report <- weighted$report()
  expect_equal(weighted$fn(), report$lp_prior_w + report$lp_penalty_w +
    report$lp_rec_w + report$lp_init_rec_w + sum(report$lp_cpue_w) +
    sum(report$lp_lf_w) + sum(report$lp_wf_w), tolerance = 1e-8)
  expect_true(all(is.finite(weighted$gr())))
})

# Doubling one component's lambda adds exactly its raw term to the objective,
# leaves the raw report unchanged, and doubles the weighted report.
expect_lambda_doubles <- function(build, data, name) {
  term <- paste0("lp_", name)
  base <- build(data)
  raw <- base$report()
  expect_gt(abs(sum(raw[[term]])), 0, label = paste("raw", term))
  data$lambdas <- stats::setNames(list(2), name)
  doubled <- build(data)
  report <- doubled$report()
  expect_equal(doubled$fn() - base$fn(), sum(raw[[term]]), tolerance = 1e-8,
               label = paste("objective change for", name))
  expect_equal(report[[term]], raw[[term]], label = paste("raw", term))
  expect_equal(report[[paste0(term, "_w")]], 2 * raw[[term]],
               label = paste("weighted", term))
}

test_that("a lambda of 2 adds exactly the raw term for each component", {
  inputs <- opaka_inputs()
  opaka <- function(data) opaka_obj(inputs, data = data)
  for (name in c("prior", "rec", "init_rec", "cpue")) {
    expect_lambda_doubles(opaka, inputs$data, name)
  }
  # Opaka's LF counts are fractional, and its multinomial tape rounds them, so
  # its reported lp_lf differs slightly from the taped term; use integer counts.
  for (name in c("lf", "wf")) {
    expect_lambda_doubles(function(data) synth_obj(data), synth_full_data(), name)
  }
  # A low B0 makes the catch exceed the population, so the penalty is active.
  penalised <- synth_full_data(lf_switch = 0L, wf_switch = 0L)
  parameters <- synth_parameters(penalised)
  parameters$log_B0 <- 14
  expect_lambda_doubles(function(data) synth_obj(data, parameters), penalised, "penalty")
})

test_that("all-one lambdas preserve bundled quickstart and IO yellowfin objectives", {
  quickstart <- opaka_inputs()
  quickstart_obj <- opaka_obj(quickstart)
  quickstart$data$lambdas <- list(cpue = 1, lf = 1, wf = 1, rec = 1,
    init_rec = 1, prior = 1, penalty = 1)
  quickstart_all_one <- opaka_obj(quickstart, data = quickstart$data)
  expect_identical(quickstart_all_one$fn(), quickstart_obj$fn())

  utils::data("io_yft_data", package = "opal", envir = environment())
  utils::data("io_yft_parameters", package = "opal", envir = environment())
  utils::data("io_yft_map", package = "opal", envir = environment())
  io_obj <- RTMB::MakeADFun(cmb(opal_model, io_yft_data),
    io_yft_parameters, map = io_yft_map, silent = TRUE)
  io_yft_data$lambdas <- list(cpue = 1, lf = 1, wf = 1, rec = 1,
    init_rec = 1, prior = 1, penalty = 1)
  io_all_one <- RTMB::MakeADFun(cmb(opal_model, io_yft_data),
    io_yft_parameters, map = io_yft_map, silent = TRUE)
  expect_identical(io_all_one$fn(), io_obj$fn())
})

test_that("composition lambdas align with fishery-grouped observation order", {
  data <- synth_full_data(wf_switch = 1L, lf_switch = 1L)
  data$lf_n_f <- data$wf_n_f <- c(1L, 1L)
  data$lf_fishery_f <- data$wf_fishery_f <- c(1L, 2L)
  data$lf_year_fi <- data$wf_year_fi <- list(1L, 2L)
  data$lf_n_fi <- data$wf_n_fi <- list(15, 15)
  data$lf_n_int_fi <- data$wf_n_int_fi <- list(15L, 15L)
  data$lambdas <- list(lf = c(0, 1), wf = c(1, 0))
  parameters <- synth_parameters(data)
  map <- synth_map(parameters)
  map$log_lf_tau <- NULL
  map$log_wf_tau <- NULL
  obj <- synth_obj(data, parameters, map)
  report <- obj$report()
  expect_equal(report$lp_lf_w, report$lp_lf * c(0, 1))
  expect_equal(report$lp_wf_w, report$lp_wf * c(1, 0))
  expect_equal(sum(report$lp_lf_w), sum(report$lp_lf[2]))
  expect_equal(sum(report$lp_wf_w), sum(report$lp_wf[1]))
  expect_true(all(is.finite(obj$gr())))
})

test_that("CPUE lambdas map by index and default missing indices to one", {
  data <- synth_data()
  data$n_index <- 2L
  data$cpue_data$index <- c(1L, 2L)
  data$lambdas <- list(cpue = c(0, 2))
  report <- synth_obj(data)$report()
  expect_equal(report$lp_cpue_w, report$lp_cpue * c(0, 2))
  expect_equal(get_lambdas(synth_data())$cpue, 1)
})

test_that("CPUE index checks match get_cpue_like() and only apply when CPUE is on", {
  data <- synth_data()
  data$cpue_data$index <- c(1L, 2L)
  data$n_index <- 1L
  # A stale n_index is extended, as get_cpue_like() does.
  expect_equal(get_lambdas(data)$cpue, c(1, 1))
  expect_true(is.finite(synth_obj(data)$fn()))
  data$cpue_data$index <- c(1, 1.5)
  expect_error(get_lambdas(data), "positive integers")
  data$cpue_switch <- 0L
  expect_equal(get_lambdas(data)$cpue, 1)
})

test_that("lambda values are validated and scalar values expand", {
  data <- synth_full_data()
  expect_equal(get_lambdas(data)$lf, c(1, 1))
  expect_equal(get_lambdas(c(data, list(lambdas = list(lf = 2))))$lf, c(2, 2))
  expect_error(get_lambdas(c(data, list(lambdas = list(lff = 1)))), "Unknown lambda")
  expect_error(get_lambdas(c(data, list(lambdas = list(lff = NULL)))), "Unknown lambda")
  expect_error(get_lambdas(c(data, list(lambdas = list(lf = -1)))), "non-negative")
  expect_error(get_lambdas(c(data, list(lambdas = list(lf = NA_real_)))), "finite")
  expect_error(get_lambdas(c(data, list(lambdas = list(lf = c(1, 2, 3))))),
                           "length one or 2")
  expect_equal(get_lambdas(c(data, list(lambdas = list(lf = NULL))))$lf,
               c(1, 1))
})

test_that("lambda weighting refuses mismatched lengths", {
  weight <- opal:::.opal_weight_likelihood
  expect_equal(weight(c(1, 2, 3), c(0, 1, 2)), c(0, 2, 6))
  expect_equal(weight(c(1, 2), 1), c(1, 2))
  expect_error(weight(c(1, 2, 3), c(0, 1)), "does not match")
  expect_error(weight(c(1, 2, 3), numeric()), "does not match")
})

test_that("switched-off components ignore per-observation lambdas", {
  data <- synth_full_data(wf_switch = 0L, lf_switch = 0L)
  data$lambdas <- list(lf = c(0, 1), wf = c(1, 0))
  report <- synth_obj(data)$report()
  expect_identical(report$lp_lf_w, 0)
  expect_identical(report$lp_wf_w, 0)
})

test_that("mapped-only priors are resolved before staged runtime construction", {
  data <- make_opal_fit_data()
  parameters <- make_opal_fit_parameters(data)
  map <- make_opal_fit_map(parameters)
  data$priors <- list(log_cpue_q = list(type = "normal", par1 = 2,
    par2 = 1, index = which(names(parameters) == "log_cpue_q")))
  ordinary <- opal_build(opal_obj(data, parameters, map))
  data$lambdas <- list(prior_mapped_only = TRUE)
  mapped_only <- opal_build(opal_obj(data, parameters, map))
  expect_gt(opal_report(ordinary)$lp_prior, 0)
  expect_no_warning(report <- opal_report(mapped_only))
  expect_equal(report$lp_prior, 0)
})

test_that("direct MakeADFun warns that prior_mapped_only was not applied", {
  data <- make_opal_fit_data()
  parameters <- make_opal_fit_parameters(data)
  map <- make_opal_fit_map(parameters)
  data$priors <- list(log_cpue_q = list(type = "normal", par1 = 2,
    par2 = 1, index = which(names(parameters) == "log_cpue_q")))
  data$lambdas <- list(prior_mapped_only = TRUE)
  expect_warning(obj <- RTMB::MakeADFun(cmb(opal_model, data), parameters,
    map = map, silent = TRUE), "not resolved against the parameter map")
  expect_gt(suppressWarnings(obj$report())$lp_prior, 0)
})

test_that("legacy conversion refuses unresolved mapped-only priors", {
  data <- make_opal_fit_data()
  parameters <- make_opal_fit_parameters(data)
  map <- make_opal_fit_map(parameters)
  data$priors <- list(log_cpue_q = list(type = "normal", par1 = 2,
    par2 = 1, index = which(names(parameters) == "log_cpue_q")))
  data$lambdas <- list(prior_mapped_only = TRUE)
  obj <- suppressWarnings(RTMB::MakeADFun(cmb(opal_model, data), parameters,
    map = map, silent = TRUE))
  opt <- list(par = obj$par, objective = suppressWarnings(obj$fn(obj$par)),
              convergence = 0L)
  legacy <- suppressWarnings(opal_fit(data, obj, opt, check = FALSE))
  expect_error(opal_from_fit(legacy), "legacy fits do not resolve")
})

test_that("MCMC joint sampling resolves mapped-only priors", {
  x <- small_opal_object(random = TRUE)
  data <- x$data
  data$priors <- list(log_cpue_q = list(type = "normal", par1 = 2, par2 = 1,
    index = which(names(x$parameters) == "log_cpue_q")))
  data$lambdas <- list(prior_mapped_only = TRUE)
  x <- opal_update(x, data = data)
  captured <- NULL
  capture <- function(obj, ...) {
    captured <<- list(random = obj$env$random, lp_prior = obj$report()$lp_prior)
    stop("captured")
  }
  expect_warning(opal_mcmc(x, sampler = capture), "captured")
  expect_length(captured$random, 0L)
  expect_equal(captured$lp_prior, 0)
})

test_that("lambda updates invalidate staged fits and warn about unidentified parameters", {
  x <- small_opal_object()
  x <- opal_build(x)
  updated <- opal_update(x, lambdas = list(cpue = 0))
  expect_null(updated$build)
  expect_null(updated$fit$opt)
  map <- updated$map
  map$log_cpue_q <- NULL
  expect_warning(opal:::.opal_warn_lambda_identifiability(
    updated$data, updated$parameters, map),
    "map the corresponding parameters off")
})

test_that("fit checks warn once and record parameters a zero lambda leaves unidentified", {
  collect <- function(expr) {
    messages <- character()
    value <- withCallingHandlers(expr, warning = function(w) {
      messages <<- c(messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
    list(value = value, n = sum(grepl("unidentified", messages)))
  }
  x <- opal_build(small_opal_object())
  map <- x$map
  map$log_cpue_q <- NULL
  x <- opal_update(x, map = map, bounds = NULL, lambdas = list(cpue = 0))
  fitted <- collect(opal_fit(x))
  expect_identical(fitted$n, 1L)
  expect_identical(fitted$value$validation$fit$metrics$lambda_identifiability,
                   "CPUE index 1")
  checked <- collect(opal_check(fitted$value))
  expect_identical(checked$n, 1L)
  clean <- collect(opal_check(opal_fit(small_opal_object(), check = FALSE)))
  expect_identical(clean$n, 0L)
  expect_identical(clean$value$validation$fit$metrics$lambda_identifiability,
                   character())
})

test_that("prior lambda helper updates an opal object and penalty lambdas warn", {
  data <- make_opal_fit_data()
  parameters <- make_opal_fit_parameters(data)
  map <- make_opal_fit_map(parameters)
  data$priors <- list(log_cpue_q = list(type = "normal", par1 = 0,
    par2 = 1, index = which(names(parameters) == "log_cpue_q")))
  x <- opal_obj(data, parameters, map)
  updated <- opal_set_prior_lambda(x, "log_cpue_q", 0)
  expect_equal(updated$data$priors$log_cpue_q$lambda, 0)
  unresolved_prior <- list(log_lf_tau = list(type = "normal",
    par1 = c(0, 0), par2 = 1,
    index = which(names(parameters) == "log_lf_tau")))
  data$priors <- unresolved_prior
  unresolved <- opal_obj(data)
  unresolved <- opal_set_prior_lambda(unresolved, "log_lf_tau", 0, element = 2)
  expect_equal(unresolved$data$priors$log_lf_tau$lambda, c(1, 0))
  expect_error(opal_set_prior_lambda(unresolved, "log_lf_tau", 0, element = 3),
               "has 2 elements")
  data$priors$log_lf_tau$par1 <- 0
  scalar <- opal_obj(data)
  expect_error(opal_set_prior_lambda(scalar, "log_lf_tau", 0, element = 2),
               "length of `log_lf_tau` is unknown")
  expect_equal(opal_set_prior_lambda(scalar, "log_lf_tau", 0)$data$priors$log_lf_tau$lambda, 0)
  resolved <- opal_obj(data, parameters, map)
  expect_equal(opal_set_prior_lambda(resolved, "log_lf_tau", 0, element = 2)$data$priors$log_lf_tau$lambda,
               c(1, 0))
  expect_error(opal_set_prior_lambda(resolved, "log_lf_tau", 0, element = 3), "has 2 elements")
  data$lambdas <- list(penalty = 0.5)
  expect_warning(opal_obj(data, parameters, map), "penalty lambda below 1")
  expect_warning(opal:::.opal_warn_tempered(
    list(n_fishery = 2L, n_index = 1L, lambdas = list(lf = 0.5)),
    context = "posterior"), "temper")
})
