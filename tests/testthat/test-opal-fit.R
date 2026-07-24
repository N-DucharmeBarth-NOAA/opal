test_that("opal_fit stores portable fitted state and posterior draws", {
  fit <- make_opal_fit_fixture()

  expect_s3_class(fit, "opal_fit")
  expect_identical(fit$schema_version, 1L)
  expect_identical(fit$model$name, "opal_model")
  expect_named(fit$bounds, c("lower", "upper"))
  expect_s3_class(fit$mcmc, "opal_mcmc")
  expect_equal(dim(fit$mcmc$samples), c(2L, 1L, 2L))
  expect_identical(fit$mcmc$parameter_scope, "active")
  expect_named(fit$derived, "projection")
  expect_false(opal:::.opal_contains_nonportable(
    opal:::.opal_fit_runtime_payload(fit)
  ))
  expect_invisible(validate_opal_fit(fit))
})

test_that("opal_fit rebuilds the objective and report", {
  fit <- make_opal_fit_fixture()
  obj <- opal_fit_object(fit, fresh = TRUE)

  expect_true(is.function(obj$fn))
  expect_equal(obj$fn(obj$par), fit$fit$opt$objective, tolerance = 1e-8)
  expect_named(opal_fit_report(fit))
  expect_true(opal_fit_compatibility(fit)$compatible)
})

test_that("opal_fit has concise print and summary methods", {
  fit <- make_opal_fit_fixture()

  expect_output(print(fit), "<opal_fit>")
  result <- summary(fit)
  expect_s3_class(result, "summary.opal_fit")
  expect_equal(result$optimization$n_parameters, length(fit$fit$opt$par))
  expect_output(print(result), "opal fitted-model summary")
  expect_output(print(fit$mcmc), "<opal_mcmc>")
})

test_that("opal_fit round trips without serializing the RTMB objective", {
  fit <- make_opal_fit_fixture()
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)

  expect_invisible(save_opal_fit(fit, path))
  expect_error(save_opal_fit(fit, path), "already exists")
  expect_invisible(save_opal_fit(fit, path, overwrite = TRUE))

  restored <- read_opal_fit(path, strict = TRUE)
  expect_s3_class(restored, "opal_fit")
  expect_equal(restored$fit$opt, fit$fit$opt)
  expect_equal(restored$mcmc$samples, fit$mcmc$samples)
  expect_true(is.function(opal_fit_object(restored)$fn))
})

test_that("opal_fit detects changed serialized payloads", {
  fit <- make_opal_fit_fixture()
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  save_opal_fit(fit, path)

  changed <- readRDS(path)
  changed$data$n_year <- 3L
  saveRDS(changed, path)

  expect_error(read_opal_fit(path), "runtime identity")
  expect_error(opal_fit_object(changed), "modified")
})

test_that("update_opal_fit manages attached results", {
  fit <- make_opal_fit_fixture(mcmc = FALSE)
  old_runtime_id <- fit$runtime_id
  active_names <- opal:::.opal_expand_parameter_names(names(fit$fit$opt$par))
  posterior <- matrix(
    c(fit$fit$opt$par, -fit$fit$opt$objective),
    nrow = 1L,
    dimnames = list(NULL, c(active_names, "lp__"))
  )

  updated <- update_opal_fit(
    fit,
    mcmc = posterior,
    diagnostics = list(max_gradient = 1e-7),
    derived = list(retrospective = data.frame(peel = 1L, rho = 0.02)),
    metadata = list(analyst = "test")
  )

  expect_s3_class(updated$mcmc, "opal_mcmc")
  expect_equal(updated$fit$diagnostics$max_gradient, 1e-7)
  expect_named(updated$derived, c("projection", "retrospective"))
  expect_identical(updated$provenance$metadata$analyst, "test")
  expect_false(identical(updated$runtime_id, old_runtime_id))
  expect_null(update_opal_fit(updated, mcmc = NULL)$mcmc)
})

test_that("opal_as_tmbfit provides SparseNUTS-compatible stored output", {
  fit <- make_opal_fit_fixture()
  converted <- opal_as_tmbfit(fit)

  expect_s3_class(converted, "tmbfit")
  expect_identical(converted$par_names, fit$mcmc$par_names)
  expect_equal(converted$warmup, 1L)
  expect_equal(dim(converted$samples)[1L], 3L)

  no_mcmc <- make_opal_fit_fixture(mcmc = FALSE)
  expect_error(opal_as_tmbfit(no_mcmc), "no stored MCMC")
})

test_that("opal_fit rejects malformed optimization and posterior state", {
  fit <- make_opal_fit_fixture(mcmc = FALSE)
  obj <- opal_fit_object(fit)
  bad_opt <- fit$fit$opt
  bad_opt$objective <- bad_opt$objective + 1

  expect_error(
    opal_fit(fit$data, obj, bad_opt),
    "does not reproduce"
  )

  bad_samples <- matrix(
    c(fit$fit$opt$par, -fit$fit$opt$objective),
    nrow = 1L
  )
  colnames(bad_samples) <- c("wrong_parameter", "lp__")
  expect_error(
    update_opal_fit(fit, mcmc = bad_samples),
    "do not match"
  )
})
