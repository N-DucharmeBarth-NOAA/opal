# Integration tests for weight likelihood in opal_model() AD tape
library(RTMB)

# Shared fixture helpers -------------------------------------------------------

# Build synthetic data with both LF and WF composition data (small, fast).
make_synthetic_full_data <- function(wf_switch = 1L, lf_switch = 1L) {
  d <- make_synthetic_data()
  
  # Disable removal filtering so composition likelihoods are actually used
  d$removal_switch_f <- rep(0L, d$n_fishery)
  
  # Add minimal LF composition for 1 fishery, 2 years
  if (lf_switch > 0L) {
    d$lf_switch <- lf_switch
    d$n_lf <- 2L  # 1 fishery × 2 years
    d$lf_year    <- c(1L, 2L)  # 1-based model timestep indices, not calendar years
    d$lf_season  <- c(1L, 1L)
    d$lf_fishery <- c(1L, 1L)
    d$lf_fishery_f <- 1L  # only fishery 1 has LF data
    d$lf_n_f <- 2L        # fishery 1 has 2 observations
    d$lf_minbin  <- c(1L, 1L)
    d$lf_maxbin  <- c(15L, 15L)
    # Flat vector of observations (2 obs × 15 bins)
    d$lf_obs     <- c(rep(1, 15), rep(1, 15))  # uniform counts
    d$lf_n       <- c(15, 15)  # total counts per observation
    d$lf_var_adj <- c(1.0, 1.0)
    # Prepare flattened obs vectors for different likelihood types
    d$lf_obs_flat <- c(rep(1, 15), rep(1, 15))
    d$lf_obs_ints <- c(rep(1L, 15), rep(1L, 15))
    d$lf_obs_prop <- d$lf_obs_flat / c(15, 15)
  }
  
  # Add minimal WF composition for 1 fishery, 2 years
  if (wf_switch > 0L) {
    d$wf_switch <- wf_switch
    d$n_wf <- 2L
    d$wf_year    <- c(1L, 2L)  # 1-based model timestep indices, not calendar years
    d$wf_season  <- c(1L, 1L)
    d$wf_fishery <- c(1L, 1L)
    d$wf_minbin  <- c(1L, 1L)
    d$wf_maxbin  <- c(15L, 15L)  # Match n_len = 15
    # Flat vector of observations (2 obs × 15 bins)
    d$wf_obs_flat <- c(rep(1, 15), rep(1, 15))
    d$wf_obs_ints <- c(rep(15L, 15), rep(15L, 15))  # denom for Dirichlet
    d$wf_obs_prop <- d$wf_obs_flat / c(15, 15)
    d$wf_n_f      <- 1L  # 1 fishery with data
    d$wf_fishery_f <- 1L
    d$wf_fishery  <- c(1L, 1L)
    d$wf_year     <- c(1L, 2L)  # 1-based model timestep indices
    d$wf_n        <- c(15L, 15L)
    # Weight binning
    d$wt_bin_start <- 1
    d$wt_bin_width <- 1
    d$n_wt <- 15L  # Match n_len for simplicity
    d$wf_rebin_matrix <- diag(15)  # Identity: 1-to-1 mapping from length to weight bins
  }
  
  d
}

# Build the full data list with both LF and WF composition data (uses wcpo_bet_data).
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

# Build minimal synthetic data for fast testing (avoids the slow wcpo_bet_data).
# Dimensions: 2 years, 1 season, 2 fisheries, 5 ages, 15 length bins.
make_synthetic_data <- function() {
  list(
    # Dimensions
    n_year    = 2L,
    n_season  = 1L,
    n_age     = 5L,
    n_fishery = 2L,
    n_len     = 15L,
    min_age   = 1L,
    max_age   = 5L,
    first_yr  = 1L,
    first_yr_catch = 1L,
    last_yr   = 2L,
    
    # Length structure (15 bins starting at 20 cm, width 2 cm)
    len_bin_start = 20,
    len_bin_width = 2,
    len_lower = seq(20, by = 2, length.out = 15),
    len_upper = seq(22, by = 2, length.out = 15),
    len_mid   = seq(21, by = 2, length.out = 15),
    
    # Age structure
    A1 = 1L,
    A2 = 5L,
    
    # Biology (scalar or age-based)
    M = rep(0.3, 5),
    maturity = c(0, 0.2, 0.5, 0.8, 1.0),
    fecundity = c(0, 100, 500, 1000, 1500),
    lw_a = 0.00001,
    lw_b = 3.0,
    
    # Catch observations (2 years × 1 season × 2 fisheries)
    catch_obs_ysf = array(c(100, 200, 150, 180), dim = c(2, 1, 2)),
    catch_units_f = c(1L, 1L),  # 1 = weight, 2 = numbers
    removal_switch_f = c(0L, 0L),  # 0 = use composition data, 1 = skip (removal only)
    
    # Selectivity structure
    sel_type_f = c(1L, 1L),  # 1 = logistic, 2 = double-normal
    
    # CPUE (minimal)
    cpue_switch = 1L,
    cpue_data = data.frame(
      ts = c(1L, 2L),
      fishery = c(1L, 1L),
      value = c(0.5, 0.48),
      se = c(0.1, 0.1),
      units = c(1L, 1L)  # 1 = weight, 2 = numbers
    ),
    
    # LF composition (will be added by tests if needed)
    lf_switch = 0L,
    n_lf = 0L,
    
    # WF composition (will be added by tests if needed)
    wf_switch = 0L,
    n_wf = 0L,
    
    # Growth (Schnute VB parameters, as log values)
    # These are dummy values; actual tests override them
    log_L1 = log(30),
    log_L2 = log(60),
    log_k  = log(0.2)
  )
}

# Build the parameters list given a data object.
make_parameters <- function(d) {
  # Create minimal parameters matching data dimensions
  list(
    log_B0        = 20,
    log_h         = 0.7,  # steepness ~ 2
    log_sigma_r   = log(0.6),
    log_cpue_q    = rep(0, 1),  # 1 CPUE index
    cpue_creep    = rep(0, 1),
    log_cpue_tau  = rep(log(0.2), 1),
    log_cpue_omega = rep(log(0.1), 1),
    log_lf_tau    = rep(log(0.1), d$n_fishery),
    log_wf_tau    = rep(0, d$n_fishery),
    log_L1        = log(30),
    log_L2        = log(60),
    log_k         = log(0.2),
    log_CV1       = log(0.1),
    log_CV2       = log(0.05),
    par_sel       = matrix(c(40, 5, 0, 0, 0, 0,
                             45, 4, 0, 0, 0, 0), 
                           nrow = d$n_fishery, ncol = 6, byrow = TRUE),
    rdev_y        = rep(0, d$n_year)
  )
}

# Build the parameter map (everything fixed except log_B0).
# log_lf_tau and log_wf_tau are only included when present in parameters.
make_map <- function(parameters) {
  m <- list(
    log_h          = factor(NA),
    log_sigma_r    = factor(NA),
    log_cpue_q     = factor(NA),
    cpue_creep     = factor(NA),
    log_cpue_tau   = factor(NA),
    log_cpue_omega = factor(NA),
    log_L1         = factor(NA),
    log_L2         = factor(NA),
    log_k          = factor(NA),
    log_CV1        = factor(NA),
    log_CV2        = factor(NA),
    par_sel        = factor(matrix(NA, nrow(parameters$par_sel),
                                   ncol(parameters$par_sel))),
    rdev_y         = factor(rep(NA, length(parameters$rdev_y)))
  )
  if (!is.null(parameters$log_lf_tau)) {
    m$log_lf_tau <- factor(rep(NA, length(parameters$log_lf_tau)))
  }
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

# Tests: full model with WF data and gradient check ---------------------------

# Build once, reuse across related tests (using fast synthetic data)
d_full     <- make_synthetic_full_data(wf_switch = 1L, lf_switch = 1L)
params_full <- make_parameters(d_full)
map_full   <- make_map(params_full)
obj_full   <- make_obj(d_full, params_full, map_full)

test_that("obj$fn() is finite with LF and WF data active", {
  expect_true(is.finite(obj_full$fn()))
})

test_that("obj$gr() is finite with LF and WF data active", {
  gr <- obj_full$gr()
  expect_true(all(is.finite(gr)))
  # Only log_B0 is free in the map above
  expect_equal(length(gr), 1L)
})

# Tests: WF-only model (lf_switch = 0) ----------------------------------------

test_that("model runs with wf_switch = 1 and lf_switch = 0", {
  d          <- make_synthetic_full_data(wf_switch = 1L, lf_switch = 0L)
  parameters <- make_parameters(d)
  map        <- make_map(parameters)
  obj        <- make_obj(d, parameters, map)
  expect_true(is.finite(obj$fn()))
  rpt <- obj$report()
  expect_true(sum(rpt$lp_wf) > 0)
})

test_that("lp_lf and lp_wf both contribute to NLL, and lp_wf is reported", {
  obj_full$fn()
  rpt <- obj_full$report()
  expect_true(sum(rpt$lp_lf) > 0)
  expect_true(sum(rpt$lp_wf) > 0)
  # Also check that obj$report()$lp_wf is returned and is a numeric vector
  expect_true("lp_wf" %in% names(rpt))
  expect_true(is.numeric(rpt$lp_wf))
  expect_true(length(rpt$lp_wf) > 0)
})

# Tests: WF disabled -----------------------------------------------------------

test_that("lp_wf is 0 when wf_switch = 0 (data prepared with switch off)", {
  d          <- make_synthetic_full_data(wf_switch = 0L, lf_switch = 1L)
  parameters <- make_parameters(d)
  map        <- make_map(parameters)
  obj        <- make_obj(d, parameters, map)
  obj$fn()
  rpt <- obj$report()
  expect_equal(sum(rpt$lp_wf), 0)
})

test_that("lp_wf is 0 when no WF data prepared (backward-compat default)", {
  d          <- make_synthetic_data()
  d$lf_switch <- 1L
  d$n_lf      <- 0L
  parameters <- make_parameters(d)
  # Remove log_wf_tau from parameters (not needed without WF data)
  parameters$log_wf_tau <- NULL
  map        <- make_map(parameters)
  obj        <- make_obj(d, parameters, map)
  obj$fn()
  rpt <- obj$report()
  expect_equal(sum(rpt$lp_wf), 0)
})

test_that("lp_lf is 0 when no LF data prepared (backward-compat default)", {
  # Build data with length structure but no LF composition fields
  d <- make_synthetic_data()
  d$lf_switch <- 0L
  d$n_lf      <- 0L
  parameters <- make_parameters(d)
  # Remove log_lf_tau since there is no LF data
  parameters$log_lf_tau <- NULL
  map        <- make_map(parameters)
  obj        <- make_obj(d, parameters, map)
  obj$fn()
  rpt <- obj$report()
  expect_equal(sum(rpt$lp_lf), 0)
})

test_that("lp_cpue is 0 when cpue_switch absent (backward-compat default)", {
  # Build data with length structure but no CPUE/LF composition fields
  d <- make_synthetic_data()
  d$lf_switch <- 0L
  d$n_lf      <- 0L
  # Remove cpue_switch to exercise the backward-compat guard
  d$cpue_switch <- NULL
  parameters <- make_parameters(d)
  parameters$log_lf_tau <- NULL
  map        <- make_map(parameters)
  obj        <- make_obj(d, parameters, map)
  obj$fn()
  rpt <- obj$report()
  expect_equal(sum(rpt$lp_cpue), 0)
})

# Tests: optimization smoke test -----------------------------------------------

test_that("5 nlminb iterations complete without error", {
  opt <- nlminb(start = obj_full$par, objective = obj_full$fn, gradient = obj_full$gr,
                control = list(iter.max = 5, eval.max = 5))
  expect_true(is.finite(opt$objective))
  expect_true(opt$convergence %in% c(0L, 1L))
})
