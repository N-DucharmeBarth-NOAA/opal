test_that("converged but constrained initial equilibrium fails biological acceptance", {
  d <- make_opal_fit_data()
  d$catch_obs_ysf[] <- 0
  d$sel_fa_external <- matrix(1, d$n_fishery, d$n_age)
  p <- make_opal_fit_parameters(d)
  p$log_init_F_f <- rep(log(0.4), d$n_fishery)
  m <- make_opal_fit_map(p)
  m$log_B0 <- factor(NA)
  m$log_init_F_f <- factor(rep(NA, d$n_fishery))
  m$log_cpue_q <- NULL
  x <- opal_fit(opal_obj(d, p, m), check = FALSE)
  expect_warning(checked <- opal_check(x), "check failed")
  expect_lt(checked$validation$fit$metrics$max_gradient, 1e-3)
  expect_false(checked$validation$fit$passes)
  expect_false(checked$validation$fit$metrics$biology$checks[["initial_equilibrium"]])
  expect_gt(opal_report(x)$lp_init_penalty, 1)
})

test_that("biological validation checks catch, harvest, mortality, and steepness", {
  x <- opal_fit(small_opal_object())
  expect_true(opal_diagnose(x)$passes)
  r <- opal_report(x)
  diagnose <- function(r, p = x$fit$parameters) opal:::.opal_diagnose_report(x$data, p, r)
  bad <- r; bad$catch_pred_ysf[1] <- 0
  expect_false(diagnose(bad)$checks[["catch_reconstruction"]])
  bad <- r; bad$hrate_ysa[1] <- 1.01
  expect_false(diagnose(bad)$checks[["harvest"]])
  bad <- r; bad$M_a[1] <- -0.1
  expect_false(diagnose(bad)$checks[["mortality"]])
  p <- x$fit$parameters; p$log_h <- 0.7
  expect_false(diagnose(r, p)$checks[["steepness"]])
  expect_true(opal:::.opal_validation_current(x, "fit"))
  stale <- x; stale$validation$fit$version <- "old"
  expect_identical(summary(stale)$checks$fit, "stale")
  stale <- x; stale$validation$fit$settings$catch_tolerance <- 1
  expect_identical(summary(stale)$checks$fit, "stale")
})

test_that("posterior checks inspect every retained joint draw and report identities", {
  x <- opal_fit(small_opal_object())
  m <- mock_opal_sampler(opal_rtmb(x), num_samples = 120)
  m$samples[12, 2, 1] <- 100
  y <- opal_attach_mcmc(x, m, check = FALSE)
  expect_warning(y <- opal_check(y, "mcmc", min_ess = 1), "check failed")
  biology <- y$validation$mcmc$metrics$biology
  expect_equal(biology$checked, 240L)
  expect_true(any(biology$failures$iteration == 12 & biology$failures$chain == 2))
  expect_match(biology$failures$reason[biology$failures$iteration == 12], "bounds")
  expect_false(y$validation$mcmc$passes)
  random <- opal_fit(small_opal_object(TRUE), check = FALSE)
  random <- opal_attach_mcmc(random, mock_opal_sampler(opal_rtmb(random)), check = FALSE)
  expect_error(opal:::.opal_posterior_context(random), "complete joint")
})
