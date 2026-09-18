# Integration tests for weight likelihood in opal_model() AD tape
library(RTMB)

# Tests: opal_globals ----------------------------------------------------------

test_that("opal_globals includes get_weight_like, rebin_counts, rebin_matrix", {
  g <- opal_globals()
  expect_true("get_weight_like" %in% names(g))
  expect_true("rebin_counts"    %in% names(g))
  expect_true("rebin_matrix"    %in% names(g))
})

test_that("opal_model uses external selectivity-at-age when supplied", {
  d <- synth_data()
  external_sel <- matrix(
    c(0.1, 0.3, 0.6, 0.9, 1.0,
      1.0, 0.8, 0.5, 0.2, 0.1),
    nrow = d$n_fishery,
    ncol = d$n_age,
    byrow = TRUE
  )
  d$sel_fa_external <- external_sel

  obj <- synth_obj(d)
  sel_fya <- obj$report()$sel_fya

  expect_equal(dim(sel_fya), c(d$n_fishery, d$n_year, d$n_age))
  for (y in seq_len(d$n_year)) {
    expect_equal(sel_fya[, y, ], external_sel)
  }
})

test_that("opal_model uses year-specific external selectivity array when supplied", {
  d <- synth_data()
  external_sel <- array(0, dim = c(d$n_fishery, d$n_year, d$n_age))
  external_sel[, 1, ] <- matrix(
    c(0.1, 0.3, 0.6, 0.9, 1.0,
      1.0, 0.8, 0.5, 0.2, 0.1),
    nrow = d$n_fishery,
    byrow = TRUE
  )
  external_sel[, 2, ] <- matrix(
    c(0.2, 0.4, 0.7, 0.95, 1.0,
      0.9, 0.7, 0.4, 0.15, 0.05),
    nrow = d$n_fishery,
    byrow = TRUE
  )
  d$sel_fa_external <- external_sel

  obj <- synth_obj(d)

  expect_equal(obj$report()$sel_fya, external_sel)
})

# Tests: full model with WF data and gradient check ---------------------------

test_that("obj$fn() is finite with LF and WF data active", {
  obj_full <- synth_obj(synth_full_data(wf_switch = 1L, lf_switch = 1L))
  expect_true(is.finite(obj_full$fn()))
})

test_that("opal_model objective includes all reported likelihood components", {
  obj_full <- synth_obj(synth_full_data(wf_switch = 1L, lf_switch = 1L))
  nll <- obj_full$fn()
  rep <- obj_full$report()
  expected <- rep$lp_prior + rep$lp_penalty + rep$lp_rec + rep$lp_init_rec +
    sum(rep$lp_cpue) + sum(rep$lp_lf) + sum(rep$lp_wf)
  expect_equal(nll, expected, tolerance = 1e-8)
})

test_that("initial recruitment-deviation prior is omitted when init_rdev_a is absent", {
  d <- synth_data()
  parameters <- synth_parameters(d)
  map <- synth_map(parameters)
  obj <- synth_obj(d, parameters, map)
  rpt <- obj$report()

  expect_equal(rpt$lp_init_rec, 0)
  expect_equal(rpt$init_rdev_a, rep(0, d$n_age))

  parameters$init_rdev_a <- rep(0.1, d$n_age)
  map <- synth_map(parameters)
  map$init_rdev_a <- factor(rep(NA, d$n_age))
  obj <- synth_obj(d, parameters, map)
  rpt <- obj$report()

  expect_gt(rpt$lp_init_rec, 0)
  expect_equal(rpt$init_rdev_a, parameters$init_rdev_a)
})

test_that("obj$gr() is finite with LF and WF data active", {
  obj_full <- synth_obj(synth_full_data(wf_switch = 1L, lf_switch = 1L))
  gr <- obj_full$gr()
  expect_true(all(is.finite(gr)))
  # Only log_B0 is free in the map above
  expect_equal(length(gr), 1L)
})

# Tests: WF-only model (lf_switch = 0) ----------------------------------------

test_that("model runs with wf_switch = 1 and lf_switch = 0", {
  d          <- synth_full_data(wf_switch = 1L, lf_switch = 0L)
  parameters <- synth_parameters(d)
  map        <- synth_map(parameters)
  obj        <- synth_obj(d, parameters, map)
  expect_true(is.finite(obj$fn()))
  rpt <- obj$report()
  expect_true(sum(rpt$lp_wf) > 0)
})

test_that("lp_lf and lp_wf both contribute to NLL, and lp_wf is reported", {
  obj_full <- synth_obj(synth_full_data(wf_switch = 1L, lf_switch = 1L))
  obj_full$fn()
  rpt <- obj_full$report()
  expect_true(sum(rpt$lp_lf) > 0)
  expect_true(sum(rpt$lp_wf) > 0)
  # Also check that obj$report()$lp_wf is returned and is a numeric vector
  expect_true("lp_wf" %in% names(rpt))
  expect_true(is.numeric(rpt$lp_wf))
  expect_true(length(rpt$lp_wf) > 0)
})

test_that("WF compositions for no-catch fleets use selected abundance", {
  d <- synth_full_data(wf_switch = 1L, lf_switch = 0L)
  d$catch_obs_ysf[, , 2] <- 0
  d$wf_fishery <- c(2L, 2L)
  d$wf_fishery_f <- 2L
  d$wf_n_f <- 2L
  d$sel_fa_external <- matrix(
    c(1, 1, 1, 1, 1,
      0, 0, 0, 0, 1),
    nrow = d$n_fishery,
    byrow = TRUE
  )

  parameters <- synth_parameters(d)
  map <- synth_map(parameters)
  obj <- synth_obj(d, parameters, map)
  obj$fn()
  pred <- as.numeric(obj$report()$wf_pred[[1]][1, ])

  expect_equal(sum(d$catch_obs_ysf[, , 2]), 0)
  expect_gt(max(pred) - min(pred), 1e-4)
  expect_gt(sum(pred[8:d$n_wt]), sum(pred[seq_len(7)]))
})

# Tests: WF disabled -----------------------------------------------------------

test_that("lp_wf is 0 when wf_switch = 0 (data prepared with switch off)", {
  d          <- synth_full_data(wf_switch = 0L, lf_switch = 1L)
  parameters <- synth_parameters(d)
  map        <- synth_map(parameters)
  obj        <- synth_obj(d, parameters, map)
  obj$fn()
  rpt <- obj$report()
  expect_equal(sum(rpt$lp_wf), 0)
})

test_that("lp_wf is 0 when no WF data prepared (backward-compat default)", {
  d          <- synth_data()
  d$lf_switch <- 1L
  d$n_lf      <- 0L
  parameters <- synth_parameters(d)
  # Remove log_wf_tau from parameters (not needed without WF data)
  parameters$log_wf_tau <- NULL
  map        <- synth_map(parameters)
  obj        <- synth_obj(d, parameters, map)
  obj$fn()
  rpt <- obj$report()
  expect_equal(sum(rpt$lp_wf), 0)
})

test_that("lp_lf is 0 when no LF data prepared (backward-compat default)", {
  # Build data with length structure but no LF composition fields
  d <- synth_data()
  d$lf_switch <- 0L
  d$n_lf      <- 0L
  parameters <- synth_parameters(d)
  # Remove log_lf_tau since there is no LF data
  parameters$log_lf_tau <- NULL
  map        <- synth_map(parameters)
  obj        <- synth_obj(d, parameters, map)
  obj$fn()
  rpt <- obj$report()
  expect_equal(sum(rpt$lp_lf), 0)
})

test_that("lp_cpue is 0 when cpue_switch absent (backward-compat default)", {
  # Build data with length structure but no CPUE/LF composition fields
  d <- synth_data()
  d$lf_switch <- 0L
  d$n_lf      <- 0L
  # Remove cpue_switch to exercise the backward-compat guard
  d$cpue_switch <- NULL
  parameters <- synth_parameters(d)
  parameters$log_lf_tau <- NULL
  map        <- synth_map(parameters)
  obj        <- synth_obj(d, parameters, map)
  obj$fn()
  rpt <- obj$report()
  expect_equal(sum(rpt$lp_cpue), 0)
})

# Tests: optimization smoke test -----------------------------------------------

test_that("5 nlminb iterations complete without error", {
  obj_full <- synth_obj(synth_full_data(wf_switch = 1L, lf_switch = 1L))
  opt <- nlminb(start = obj_full$par, objective = obj_full$fn, gradient = obj_full$gr,
                control = list(iter.max = 5, eval.max = 5))
  expect_true(is.finite(opt$objective))
  expect_true(opt$convergence %in% c(0L, 1L))
})
