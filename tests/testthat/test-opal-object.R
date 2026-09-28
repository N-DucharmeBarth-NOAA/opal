test_that("construction, resolution, and fitting have distinct stages", {
  x <- small_opal_object()
  expect_s3_class(x, "opal_obj")
  expect_identical(summary(x)$stage, "configured")
  expect_output(print(x), "configured")
  expect_null(x$fit$opt)
  b <- opal_build(x)
  expect_identical(summary(b)$stage, "built")
  expect_identical(x$parameters, b$parameters)
  expect_identical(x$map, b$map)
  expect_false(opal:::.opal_contains_nonportable(b))
  expect_identical(opal_rtmb(b), opal_rtmb(b))
  expect_false(identical(opal_rtmb(b), opal_rtmb(b, fresh = TRUE)))
  fit <- opal_fit(b, check = FALSE)
  expect_identical(summary(fit)$stage, "fitted")
  expect_lt(fit$fit$opt$objective, b$build$objective)
  expect_length(fit$fit$diagnostics$optimisation$passes, 2L)
  expect_equal(opal_rtmb(fit, fresh = TRUE)$fn(fit$fit$opt$par), fit$fit$opt$objective)
  checked <- suppressWarnings(opal_check(fit))
  expect_identical(checked$validation$fit$identity, fit$identity)
  expect_type(checked$validation$fit$passes, "logical")
  expect_identical(checked$validation$fit$error, NULL)
  expect_error(opal_check(x), "No fitted")
  expect_error(opal_fit(x, n_passes = 0), "positive integer")
  expect_error(opal_obj(x$data, x$parameters, list(unknown = factor(NA))), "unknown")
  expect_error(opal_obj(x$data, x$parameters, list(rdev_y = factor(1))), "Invalid map")
  expect_error(opal_obj(x$data, x$parameters, x$map, random = "absent"), "random")
  expect_error(opal_build(opal_obj(x$data)), "explicit bundled")
  expect_identical(summary(opal_obj(x$data))$stage, "data")
  expect_identical(opal_obj(x$data, x$parameters, list())$map, list())
})

test_that("configuration changes invalidate dependent state but annotations do not", {
  x <- opal_fit(small_opal_object(), check = FALSE)
  x <- opal_mcmc(x, sampler = mock_opal_sampler, chains = 2, check = FALSE)
  obj <- opal_rtmb(x)
  tagged <- opal_update(x, metadata = list(stock = "test"), control = list(iter.max = 50))
  expect_identical(tagged$identity, x$identity)
  expect_identical(tagged$mcmc, x$mcmc)
  expect_identical(opal_rtmb(tagged), obj)
  data <- x$data
  data$cpue_data$value[1] <- 0.6
  y <- opal_update(x, data = data)
  expect_null(y$fit$opt)
  expect_null(y$mcmc)
  expect_null(y$build)
  expect_identical(y$bounds, x$bounds)
  expect_identical(y$map, x$map)
  expect_identical(y$validation, list())
  expect_null(opal_update(x, priors = list(test = 1))$fit$opt)
  bounds <- x$bounds
  bounds$upper[] <- 21
  expect_null(opal_update(x, bounds = bounds)$fit$opt)
  parameters <- x$parameters
  parameters$log_B0 <- 19
  expect_null(opal_update(x, parameters = parameters)$mcmc)
  altered <- x
  altered$data <- data
  expect_error(opal_rtmb(altered), "modified directly")
  expect_identical(opal_update(x, data = x$data)$identity, x$identity)
})

test_that("fresh runtime protects fitted state and random-effect modes", {
  for (random in c(FALSE, TRUE)) {
    x <- opal_fit(small_opal_object(random), check = FALSE)
    expected <- opal_report(x)
    cache <- opal_rtmb(x)
    invisible(cache$fn(cache$par + 0.1))
    expect_equal(opal_report(x), expected)
    fresh <- opal_rtmb(x, fresh = TRUE)
    expect_equal(fresh$par, x$fit$opt$par)
    expect_equal(fresh$env$last.par.best, x$fit$last_par_best)
    expect_equal(as.numeric(fresh$fn(fresh$par)), x$fit$opt$objective, tolerance = 1e-6)
    expect_lte(opal_fit(x, check = FALSE)$fit$opt$objective, x$fit$opt$objective + 1e-6)
    expect_false(opal:::.opal_contains_nonportable(x))
  }
})

test_that("legacy conversion and every stage survive saved payload verification", {
  legacy <- make_opal_fit_fixture()
  converted <- opal_from_fit(legacy)
  expect_equal(converted$fit$opt, legacy$fit$opt)
  expect_equal(converted$mcmc$samples, legacy$mcmc$samples)
  expect_equal(opal_report(converted), opal_fit_report(legacy))
  stages <- list(opal_obj(legacy$data), small_opal_object(),
    opal_build(small_opal_object()), opal_fit(small_opal_object(), check = FALSE), converted)
  for (x in stages) {
    file <- tempfile(fileext = ".rds")
    opal_save(x, file)
    y <- opal_read(file)
    expect_identical(y, x)
    expect_false(opal:::.opal_contains_nonportable(y))
    expect_error(opal_save(x, file), "already exists")
    opal_save(x, file, overwrite = TRUE)
    expect_identical(opal_read(file), x)
    saved <- readRDS(file)
    saved$object$data$n_year <- 99L
    saveRDS(saved, file)
    expect_error(opal_read(file), "modified Opal file")
    unlink(file)
  }
  path <- tempfile(fileext = ".rds")
  save_opal_fit(legacy, path)
  expect_s3_class(opal_read(path), "opal_obj")
  unlink(path)
})

test_that("MCMC uses an isolated runtime and preserves attempts and provenance", {
  x <- opal_fit(small_opal_object(), check = FALSE)
  cached <- opal_rtmb(x)
  sampled <- opal_mcmc(x, sampler = mock_opal_sampler, chains = 2, seed = 123,
    check_args = list(max_rhat = 1.1, min_ess = 10))
  expect_true(sampled$validation$mcmc$passes)
  expect_identical(summary(sampled)$stage, "sampled")
  expect_identical(opal_rtmb(sampled), cached)
  expect_identical(sampled$fit, x$fit)
  expect_identical(sampled$mcmc$settings$seed, 123)
  again <- opal_fit(sampled, check = FALSE)
  expect_identical(again$mcmc, sampled$mcmc)
  expect_identical(again$mcmc$fit_id, sampled$identity$fit)
  expect_warning(failed <- opal_mcmc(sampled, sampler = function(...) stop("mock failure")), "mock failure")
  expect_identical(failed$mcmc, sampled$mcmc)
  expect_identical(tail(failed$mcmc_history, 1)[[1]]$status, "error")
  unchecked <- opal_mcmc(sampled, sampler = mock_opal_sampler, chains = 2, seed = 2, check = FALSE)
  expect_identical(unchecked$mcmc, sampled$mcmc)
  expect_identical(tail(unchecked$mcmc_history, 1)[[1]]$status, "not_selected")
  bad <- sampled
  bad$mcmc$samples[1] <- 999
  expect_error(validate_opal_obj(bad, results = TRUE), "modified")
  expect_error(opal_save(bad, tempfile()), "modified")
  expect_error(opal_attach_mcmc(x, matrix(1, 2, 1, dimnames = list(NULL, "bad"))), "layout")
  configured <- opal_mcmc(small_opal_object(), sampler = mock_opal_sampler, chains = 2, check = FALSE)
  expect_null(configured$fit$opt)
  expect_s3_class(opal_as_tmbfit(configured), "tmbfit")
})

test_that("object plots resolve data and fitted predictions without cache mutation", {
  x <- opal_fit(small_opal_object(), check = FALSE)
  for (fun in list(plot_catch, plot_cpue, plot_initial_numbers, plot_hrate, plot_recruitment)) {
    plot <- suppressMessages(fun(x))
    expect_s3_class(plot, "ggplot")
    expect_s3_class(ggplot2::ggplot_build(plot), "ggplot_built")
  }
  expect_s3_class(plot_biomass_spawning(x), "ggplot")
  expect_s3_class(plot_biomass_spawning(list(x, x), labels = c("A", "B")), "ggplot")
  expect_error(plot_cpue(x, opal_rtmb(x)), "separate objective")
})

test_that("the new interface agrees with the unchanged Opakapaka reference", {
  inputs <- opaka_inputs()
  x <- opal_build(opal_obj(inputs$data, inputs$parameters, inputs$map))
  reference <- readRDS(test_path("_reference", "opaka-quickstart.rds"))
  for (point in c("start", "mle")) {
    ref <- reference[[point]]
    obj <- opal_rtmb(x, fresh = TRUE)
    expect_equal(obj$fn(ref$par), ref$nll, tolerance = reference_tolerance())
    if (!is.null(ref$gr)) expect_equal(as.vector(obj$gr(ref$par)), ref$gr, tolerance = reference_tolerance())
    report <- obj$report(ref$par)
    for (name in names(ref$report)) expect_equal(report[[name]], ref$report[[name]], tolerance = reference_tolerance())
  }
  fit <- suppressWarnings(opal_fit(x))
  expect_equal(fit$fit$opt$objective, reference$mle$nll, tolerance = reference_tolerance())
  expect_identical(fit$model$scientific_version, reference$meta$scientific_version)
  # Default map resolution is separate from an explicitly empty map.
  defaults <- opal_build(opal_obj(inputs$data, inputs$parameters))
  expect_identical(defaults$configuration$origins$map, "default")
  parameters <- defaults$parameters
  parameters$log_B0 <- parameters$log_B0 + 0.01
  changed <- opal_update(defaults, parameters = parameters)
  expect_null(changed$map)
  expect_null(changed$bounds)
})

test_that("projections require an explicit posterior choice and retain provenance", {
  x <- opal_fit(small_opal_object(), check = FALSE)
  selectivity <- project_selectivity(x, n_proj = 2, n_iter = 2)
  args <- list(n_proj = 2L, n_iter = 2L, rdev_y = matrix(0, 2, 2),
    sel_fya = selectivity, catch_ysf = array(0, c(2, 1, 2)), return_hist = TRUE)
  x <- opal_mcmc(x, sampler = mock_opal_sampler, chains = 2, check = FALSE)
  expect_error(do.call(project_dynamics, c(list(data = x), args)), "Choose uncertainty")
  result <- suppressMessages(do.call(opal_project, c(list(x = x, uncertainty = "mcmc", seed = 23), args)))
  expect_length(result$derived$projection$result$dyn, 2)
  expect_identical(result$derived$projection$identity, opal:::.opal_check_identity(x, "mcmc"))
  expect_identical(result$derived$projection$settings$seed, 23)
  expect_identical(opal_update(result, metadata = list(label = "test"))$derived, result$derived)
  expect_identical(opal_update(result, priors = list())$derived, list())
  random <- opal_fit(small_opal_object(TRUE), check = FALSE)
  expect_error(opal:::.opal_projection_input(random, uncertainty = "mvn"), "random effects")
})

test_that("short real SparseNUTS runs work before and after fitting", {
  skip_on_cran()
  for (mode in c("unfitted", "fitted", "random")) {
    x <- small_opal_object(random = mode == "random")
    if (mode != "unfitted") x <- opal_fit(x, check = FALSE)
    sampled <- suppressWarnings(opal_mcmc(x, seed = 510, chains = 2L, cores = 1L,
      num_samples = 20L, num_warmup = 20L, refresh = 0, print = FALSE, check = FALSE))
    expect_false(is.null(sampled$mcmc), info = mode)
    expect_identical(sampled$fit, x$fit)
    expect_identical(sampled$mcmc$iter, 20L)
    expect_identical(sampled$mcmc$chains, 2L)
    expect_identical(sampled$mcmc$parameter_scope, if (mode == "random") "complete" else "active")
    expect_silent(validate_opal_obj(sampled, results = TRUE))
  }
})

test_that("failed or incomplete diagnostics never certify a posterior", {
  x <- opal_build(small_opal_object())
  m <- mock_opal_sampler(opal_rtmb(x))
  m$sampler_params[[1]][20, "divergent__"] <- 1
  expect_warning(y <- opal_attach_mcmc(x, m), "check failed")
  expect_false(y$validation$mcmc$passes)
  expect_identical(y$validation$mcmc$metrics$divergences, 1)
  expect_false(is.null(y$mcmc))
  m$sampler_params <- NULL
  expect_warning(y <- opal_attach_mcmc(x, m), "check failed")
  expect_false(y$validation$mcmc$metrics$sampler_diagnostics_known)
  expect_error(opal_check(y, scope = "mcmc", stop_on_failure = TRUE), "check failed")
})

test_that("all stages can be loaded in a fresh installed-package R process", {
  package_path <- system.file(package = "opal")
  skip_if_not(file.exists(file.path(package_path, "Meta", "package.rds")),
              "Fresh-process loading runs against the installed package during R CMD check")
  path <- tempfile("opal-stages-")
  dir.create(path)
  on.exit(unlink(path, recursive = TRUE), add = TRUE)
  x <- small_opal_object()
  fitted <- opal_fit(x, check = FALSE)
  stages <- list(opal_obj(x$data), x, opal_build(x), fitted,
    opal_mcmc(fitted, sampler = mock_opal_sampler, check = FALSE),
    opal_fit(small_opal_object(TRUE), check = FALSE))
  files <- file.path(path, paste0(seq_along(stages), ".rds"))
  for (i in seq_along(stages)) opal_save(stages[[i]], files[i])
  script <- file.path(path, "read.R")
  writeLines(c(
    sprintf(".libPaths(c(%s, .libPaths()))", paste(deparse(dirname(package_path)), collapse = "")),
    "library(opal)",
    "for (file in commandArgs(TRUE)) {",
    "  x <- opal_read(file)",
    "  validate_opal_obj(x, results = TRUE)",
    "  if (!is.null(x$parameters) && !is.null(x$map)) {",
    "    x <- opal_build(x)",
    "    stopifnot(all(is.finite(opal_report(x)$number_ysa)))",
    "  }",
    "}"), script)
  output <- system2(file.path(R.home("bin"), "Rscript"),
    c("--vanilla", shQuote(script), shQuote(files)), stdout = TRUE, stderr = TRUE)
  expect_null(attr(output, "status"), info = paste(output, collapse = "\n"))
})

test_that("external attachment retains the supplied point, not a cached better point", {
  fitted <- opal_fit(small_opal_object(), check = FALSE)
  object <- opal_rtmb(fitted, fresh = TRUE)
  par <- fitted$fit$opt$par + 0.2
  objective <- as.numeric(object$fn(par))
  expect_gt(objective, fitted$fit$opt$objective)
  attached <- opal_attach_fit(fitted, list(par = par, objective = objective,
    convergence = 1L, message = "intentional non-optimum"), check = FALSE)
  expect_equal(opal_rtmb(attached)$par, par)
  expect_equal(attached$fit$last_par_best, par)
  expect_equal(opal_report(attached), object$report(par))
  expect_warning(checked <- opal_check(attached), "check failed")
  expect_false(checked$validation$fit$passes)
  expect_identical(checked$fit$opt, attached$fit$opt)
})

test_that("recruitment projection reconstructs fixed and random-effect layouts", {
  fixed <- opal_fit(small_opal_object(), check = FALSE)
  rec <- project_rec_devs(fixed, n_proj = 2, n_iter = 2)
  expect_equal(unname(rec$rdev_y), matrix(0, 2, 2))
  sampled <- opal_mcmc(fixed, sampler = mock_opal_sampler, check = FALSE)
  expect_error(project_rec_devs(sampled), "Choose uncertainty")
  rec <- project_rec_devs(sampled, uncertainty = "mcmc", n_proj = 2, n_iter = 2)
  expect_equal(unname(rec$rdev_y), matrix(0, 2, 2))
  random <- opal_fit(small_opal_object(TRUE), check = FALSE)
  # The random-effect modes must be reconstructed even though they are absent
  # from the marginal optimiser vector.
  rec <- project_rec_devs(random, uncertainty = "fit", n_proj = 2, n_iter = 1)
  expect_equal(dim(rec$rdev_y), c(1L, 2L))
  expect_true(all(is.finite(rec$rdev_y)))
})

test_that("parameter-only posterior imports do not lose their final variable", {
  x <- opal_build(small_opal_object())
  draws <- matrix(c(15, 15.1), ncol = 1, dimnames = list(NULL, "log_B0"))
  x <- opal_attach_mcmc(x, draws, check = FALSE)
  expect_equal(as.matrix(SparseNUTS::extract_samples(opal_as_tmbfit(x))), draws,
               ignore_attr = TRUE)
  expect_true(all(is.na(opal_as_tmbfit(x)$samples[, , "lp__"])))
})

test_that("the R 4.6 object fixture remains readable across R versions", {
  x <- opal_read(test_path("_fixtures", "opal-obj-v1.rds"), rebuild = TRUE)
  expect_identical(summary(x)$stage, "sampled")
  expect_silent(validate_opal_obj(x, results = TRUE))
  expect_equal(as.numeric(opal_rtmb(x, fresh = TRUE)$fn(x$fit$opt$par)),
               x$fit$opt$objective, tolerance = 1e-6)
  expect_identical(x$mcmc$par_names, "log_B0")
})
