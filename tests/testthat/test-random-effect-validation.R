test_that("random-effect fits pass checks using the marginal Hessian", {
  x <- opal_fit(small_opal_object(TRUE), check = FALSE)
  original_fit <- x$fit
  original_report <- opal_report(x)
  object <- opal_rtmb(x, fresh = TRUE)
  hessian <- stats::optimHess(x$fit$opt$par, object$fn, object$gr)
  expect_true(all(eigen(hessian, symmetric = TRUE)$values > 0))

  checked <- opal_check(x, stop_on_failure = TRUE)
  expect_true(checked$validation$fit$passes)
  expect_null(checked$validation$fit$error)
  expect_true(checked$validation$fit$metrics$positive_hessian)
  expect_true(checked$validation$fit$metrics$biology$passes)
  expect_lt(checked$validation$fit$metrics$max_gradient, 1e-3)
  expect_identical(checked$fit$opt, original_fit$opt)
  expect_identical(checked$fit$parameters, original_fit$parameters)
  expect_identical(checked$fit$last_par_best, original_fit$last_par_best)
  expect_equal(opal_report(checked), original_report)
  expect_true(opal_fit(small_opal_object(TRUE))$validation$fit$passes)
})

test_that("random-effect checking still rejects an unsuccessful optimiser", {
  x <- opal_fit(small_opal_object(TRUE), check = FALSE)
  opt <- x$fit$opt
  opt$convergence <- 1L
  x <- opal_attach_fit(x, opt, check = FALSE)
  expect_warning(x <- opal_check(x), "fit check failed")
  expect_false(x$validation$fit$passes)
  expect_null(x$validation$fit$error)
  expect_identical(x$validation$fit$metrics$convergence, 1L)
  expect_true(x$validation$fit$metrics$positive_hessian)
  expect_true(x$validation$fit$metrics$biology$passes)
})

test_that("accepted random-effect grid fits proceed to joint sampling", {
  grid <- opal_grid(small_opal_object(TRUE), list(base = list()))
  expect_true(grid$summary$passes)
  calls <- 0L
  sampler <- function(obj, ...) {
    calls <<- calls + 1L
    mock_opal_sampler(obj, ...)
  }
  sampled <- opal_grid_mcmc(grid, sampler = sampler, chains = 2L,
                           num_samples = 20L, seed = 41L, check = FALSE)
  expect_identical(calls, 1L)
  expect_s3_class(sampled$models$base$mcmc, "opal_mcmc")
  expect_identical(sampled$models$base$mcmc$parameter_scope, "complete")
  expect_identical(sampled$models$base$fit, grid$models$base$fit)
})

test_that("older validation checkpoints are not reused after the Hessian fix", {
  directory <- tempfile()
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE))
  x <- small_opal_object(TRUE)
  old <- testthat::with_mocked_bindings(
    opal_grid(x, list(base = list()), directory),
    .opal_validation_version = "opal_validation_v2", .package = "opal")
  expect_identical(summary(old$models$base)$checks$fit, "stale")
  current <- opal_grid(x, list(base = list()), directory)
  expect_false(current$summary$reused)
  expect_true(current$summary$passes)
  expect_true(opal_grid(x, list(base = list()), directory)$summary$reused)
})
