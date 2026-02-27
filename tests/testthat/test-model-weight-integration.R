# Integration tests for weight likelihood in opal_model() AD tape
library(RTMB)

# Shared fixture helpers -------------------------------------------------------

# Build the data list with both LF and WF composition data.
make_full_data <- function(wf_switch = 1L, lf_switch = 1L,
                           wf_keep_fisheries = c(6, 7),
                           lf_keep_fisheries = c(8, 9)) {
  data(wcpo_bet_data, package = "opal", envir = environment())
  data(wcpo_bet_lf,   package = "opal", envir = environment())
  data(wcpo_bet_wf,   package = "opal", envir = environment())
  d <- wcpo_bet_data

  # LF data
  lf_wide <- wcpo_bet_lf |>
    tidyr::pivot_wider(
      id_cols     = c(fishery, year, month, ts),
      names_from  = bin,
      values_from = value,
      values_fill = 0
    ) |>
    dplyr::arrange(fishery, ts)
  d <- prep_lf_data(d, lf_wide, lf_keep_fisheries = lf_keep_fisheries,
                    lf_switch = lf_switch)

  # Weight bin scalars (required by prep_wf_data)
  d$wt_bin_start <- 1
  d$wt_bin_width <- 1
  d$n_wt         <- 200L

  # WF data
  wf_wide <- wcpo_bet_wf |>
    tidyr::pivot_wider(
      id_cols     = c(fishery, year, month, ts),
      names_from  = bin,
      values_from = value,
      values_fill = 0
    ) |>
    dplyr::arrange(fishery, ts)
  d <- prep_wf_data(d, wf_wide, wf_keep_fisheries = wf_keep_fisheries,
                    wf_switch = wf_switch)

  d
}

# Build the data list with LF data only (no weight data prepared).
make_lf_only_data <- function(lf_switch = 1L,
                              lf_keep_fisheries = c(8, 9)) {
  data(wcpo_bet_data, package = "opal", envir = environment())
  data(wcpo_bet_lf,   package = "opal", envir = environment())
  d <- wcpo_bet_data

  lf_wide <- wcpo_bet_lf |>
    tidyr::pivot_wider(
      id_cols     = c(fishery, year, month, ts),
      names_from  = bin,
      values_from = value,
      values_fill = 0
    ) |>
    dplyr::arrange(fishery, ts)
  d <- prep_lf_data(d, lf_wide, lf_keep_fisheries = lf_keep_fisheries,
                    lf_switch = lf_switch)
  d
}

# Build the parameters list given a data object.
make_parameters <- function(d) {
  data(wcpo_bet_parameters, package = "opal", envir = environment())
  list(
    log_B0        = 20,
    log_h         = as.numeric(wcpo_bet_parameters$log_h),
    log_sigma_r   = as.numeric(wcpo_bet_parameters$log_sigma_r),
    log_cpue_q    = as.numeric(wcpo_bet_parameters$log_cpue_q),
    cpue_creep    = as.numeric(wcpo_bet_parameters$cpue_creep),
    log_cpue_tau  = as.numeric(wcpo_bet_parameters$log_cpue_tau),
    log_cpue_omega = as.numeric(wcpo_bet_parameters$log_cpue_omega),
    log_lf_tau    = rep(log(0.1), d$n_fishery),
    log_wf_tau    = rep(0, d$n_fishery),
    log_L1        = as.numeric(wcpo_bet_parameters$log_L1),
    log_L2        = as.numeric(wcpo_bet_parameters$log_L2),
    log_k         = as.numeric(wcpo_bet_parameters$log_k),
    log_CV1       = as.numeric(wcpo_bet_parameters$log_CV1),
    log_CV2       = as.numeric(wcpo_bet_parameters$log_CV2),
    par_sel       = as.matrix(wcpo_bet_parameters$par_sel),
    rdev_y        = as.numeric(wcpo_bet_parameters$rdev_y)
  )
}

# Build the parameter map (everything fixed except log_B0).
# log_wf_tau is only included when present in parameters.
make_map <- function(parameters) {
  m <- list(
    log_h          = factor(NA),
    log_sigma_r    = factor(NA),
    log_cpue_q     = factor(NA),
    cpue_creep     = factor(NA),
    log_cpue_tau   = factor(NA),
    log_cpue_omega = factor(NA),
    log_lf_tau     = factor(rep(NA, length(parameters$log_lf_tau))),
    log_L1         = factor(NA),
    log_L2         = factor(NA),
    log_k          = factor(NA),
    log_CV1        = factor(NA),
    log_CV2        = factor(NA),
    par_sel        = factor(matrix(NA, nrow(parameters$par_sel),
                                   ncol(parameters$par_sel))),
    rdev_y         = factor(rep(NA, length(parameters$rdev_y)))
  )
  if (!is.null(parameters$log_wf_tau)) {
    m$log_wf_tau <- factor(rep(NA, length(parameters$log_wf_tau)))
  }
  m
}

# Build the full AD object from data + parameters + map.
make_obj <- function(d, parameters = NULL, map = NULL) {
  if (is.null(parameters)) parameters <- make_parameters(d)
  if (is.null(map))        map        <- make_map(parameters)
  suppressWarnings(
    RTMB::MakeADFun(func = cmb(opal_model, d), parameters = parameters,
                    map = map, silent = TRUE)
  )
}

# Tests: bet_globals -----------------------------------------------------------

test_that("bet_globals includes get_weight_like, rebin_counts, rebin_matrix", {
  g <- bet_globals()
  expect_true("get_weight_like" %in% names(g))
  expect_true("rebin_counts"    %in% names(g))
  expect_true("rebin_matrix"    %in% names(g))
})

# Tests: full model with WF data -----------------------------------------------

test_that("obj$fn() is finite with LF and WF data active", {
  d          <- make_full_data(wf_switch = 1L, lf_switch = 1L)
  parameters <- make_parameters(d)
  map        <- make_map(parameters)
  obj        <- make_obj(d, parameters, map)
  expect_true(is.finite(obj$fn()))
})

# Tests: gradient check --------------------------------------------------------

test_that("obj$gr() is finite with LF and WF data active", {
  d          <- make_full_data(wf_switch = 1L, lf_switch = 1L)
  parameters <- make_parameters(d)
  map        <- make_map(parameters)
  obj        <- make_obj(d, parameters, map)
  gr <- obj$gr()
  expect_true(all(is.finite(gr)))
  # Only log_B0 is free in the map above
  expect_equal(length(gr), 1L)
})

# Tests: WF-only model (lf_switch = 0) ----------------------------------------

test_that("model runs with wf_switch = 1 and lf_switch = 0", {
  d          <- make_full_data(wf_switch = 1L, lf_switch = 0L)
  parameters <- make_parameters(d)
  map        <- make_map(parameters)
  obj        <- make_obj(d, parameters, map)
  expect_true(is.finite(obj$fn()))
  rpt <- obj$report()
  expect_true(sum(rpt$lp_wf) > 0)
})

# Tests: LF + WF both active ---------------------------------------------------

test_that("lp_lf and lp_wf both contribute to NLL", {
  d          <- make_full_data(wf_switch = 1L, lf_switch = 1L)
  parameters <- make_parameters(d)
  map        <- make_map(parameters)
  obj        <- make_obj(d, parameters, map)
  obj$fn()
  rpt <- obj$report()
  expect_true(sum(rpt$lp_lf) > 0)
  expect_true(sum(rpt$lp_wf) > 0)
})

# Tests: WF disabled -----------------------------------------------------------

test_that("lp_wf is 0 when wf_switch = 0 (data prepared with switch off)", {
  d          <- make_full_data(wf_switch = 0L, lf_switch = 1L)
  parameters <- make_parameters(d)
  map        <- make_map(parameters)
  obj        <- make_obj(d, parameters, map)
  obj$fn()
  rpt <- obj$report()
  expect_equal(sum(rpt$lp_wf), 0)
})

test_that("lp_wf is 0 when no WF data prepared (backward-compat default)", {
  d          <- make_lf_only_data(lf_switch = 1L)
  parameters <- make_parameters(d)
  # Remove log_wf_tau from parameters (not needed without WF data)
  parameters$log_wf_tau <- NULL
  map        <- make_map(parameters)
  obj        <- make_obj(d, parameters, map)
  obj$fn()
  rpt <- obj$report()
  expect_equal(sum(rpt$lp_wf), 0)
})

# Tests: REPORT(lp_wf) ---------------------------------------------------------

test_that("obj$report()$lp_wf is returned and is a numeric vector", {
  d          <- make_full_data(wf_switch = 1L, lf_switch = 1L)
  parameters <- make_parameters(d)
  map        <- make_map(parameters)
  obj        <- make_obj(d, parameters, map)
  obj$fn()
  rpt <- obj$report()
  expect_true("lp_wf" %in% names(rpt))
  expect_true(is.numeric(rpt$lp_wf))
  expect_true(length(rpt$lp_wf) > 0)
})

# Tests: optimization smoke test -----------------------------------------------

test_that("10 nlminb iterations complete without error", {
  d          <- make_full_data(wf_switch = 1L, lf_switch = 1L)
  parameters <- make_parameters(d)
  map        <- make_map(parameters)
  obj        <- make_obj(d, parameters, map)
  opt <- nlminb(start = obj$par, objective = obj$fn, gradient = obj$gr,
                control = list(iter.max = 10, eval.max = 20))
  expect_true(is.finite(opt$objective))
  expect_true(opt$convergence %in% c(0L, 1L))
})
