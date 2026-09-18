synth_data <- function() {
  list(
    n_year = 2L, n_season = 1L, n_age = 5L, n_fishery = 2L, n_len = 15L,
    min_age = 1L, max_age = 5L, first_yr = 1L, first_yr_catch = 1L, last_yr = 2L,
    len_bin_start = 20, len_bin_width = 2,
    len_lower = seq(20, by = 2, length.out = 15),
    len_upper = seq(22, by = 2, length.out = 15),
    len_mid = seq(21, by = 2, length.out = 15),
    A1 = 1L, A2 = 5L,
    M = rep(0.3, 5), maturity = c(0, 0.2, 0.5, 0.8, 1.0),
    fecundity = c(0, 100, 500, 1000, 1500), lw_a = 0.00001, lw_b = 3.0,
    catch_obs_ysf = array(c(100, 200, 150, 180), dim = c(2, 1, 2)),
    catch_units_f = c(1L, 1L), removal_switch_f = c(0L, 0L),
    sel_type_f = c(1L, 1L), cpue_switch = 1L,
    cpue_data = data.frame(ts = c(1L, 2L), fishery = c(1L, 1L),
                           value = c(0.5, 0.48), se = c(0.1, 0.1),
                           units = c(1L, 1L)),
    lf_switch = 0L, n_lf = 0L, wf_switch = 0L, n_wf = 0L,
    log_L1 = log(30), log_L2 = log(60), log_k = log(0.2)
  )
}

synth_full_data <- function(wf_switch = 1L, lf_switch = 1L) {
  d <- synth_data()
  d$removal_switch_f <- rep(0L, d$n_fishery)
  if (lf_switch > 0L) {
    d$lf_switch <- lf_switch
    d$n_lf <- 2L
    d$lf_year <- c(1L, 2L)
    d$lf_season <- c(1L, 1L)
    d$lf_fishery <- c(1L, 1L)
    d$lf_fishery_f <- 1L
    d$lf_n_f <- 2L
    d$lf_minbin <- c(1L, 1L)
    d$lf_maxbin <- c(15L, 15L)
    d$lf_obs <- c(rep(1, 15), rep(1, 15))
    d$lf_n <- c(15, 15)
    d$lf_var_adj <- c(1.0, 1.0)
    d$lf_obs_flat <- c(rep(1, 15), rep(1, 15))
    d$lf_obs_ints <- c(rep(1L, 15), rep(1L, 15))
    d$lf_obs_prop <- d$lf_obs_flat / c(15, 15)
  }
  if (wf_switch > 0L) {
    d$wf_switch <- wf_switch
    d$n_wf <- 2L
    d$wf_year <- c(1L, 2L)
    d$wf_season <- c(1L, 1L)
    d$wf_fishery <- c(1L, 1L)
    d$wf_minbin <- c(1L, 1L)
    d$wf_maxbin <- c(15L, 15L)
    d$wf_obs_flat <- c(rep(1, 15), rep(1, 15))
    d$wf_obs_ints <- c(rep(15L, 15), rep(15L, 15))
    d$wf_obs_prop <- d$wf_obs_flat / c(15, 15)
    d$wf_n_f <- 1L
    d$wf_fishery_f <- 1L
    d$wf_n <- c(15L, 15L)
    d$wt_bin_start <- 1
    d$wt_bin_width <- 1
    d$n_wt <- 15L
    d$wf_rebin_matrix <- diag(15)
  }
  d
}

synth_parameters <- function(d) {
  list(
    log_B0 = 20, log_h = 0.7, log_sigma_r = log(0.6),
    log_cpue_q = rep(0, 1), cpue_creep = rep(0, 1),
    log_cpue_tau = rep(log(0.2), 1), log_cpue_omega = rep(log(0.1), 1),
    log_lf_tau = rep(log(0.1), d$n_fishery),
    log_wf_tau = rep(0, d$n_fishery),
    log_L1 = log(30), log_L2 = log(60), log_k = log(0.2),
    log_CV1 = log(0.1), log_CV2 = log(0.05),
    par_sel = matrix(c(40, 5, 0, 0, 0, 0, 45, 4, 0, 0, 0, 0),
                     nrow = d$n_fishery, ncol = 6, byrow = TRUE),
    rdev_y = rep(0, d$n_year)
  )
}

synth_map <- function(parameters) {
  map <- list(
    log_h = factor(NA), log_sigma_r = factor(NA), log_cpue_q = factor(NA),
    cpue_creep = factor(NA), log_cpue_tau = factor(NA),
    log_cpue_omega = factor(NA), log_L1 = factor(NA), log_L2 = factor(NA),
    log_k = factor(NA), log_CV1 = factor(NA), log_CV2 = factor(NA),
    par_sel = factor(matrix(NA, nrow(parameters$par_sel), ncol(parameters$par_sel))),
    rdev_y = factor(rep(NA, length(parameters$rdev_y)))
  )
  if (!is.null(parameters$log_lf_tau)) map$log_lf_tau <- factor(rep(NA, length(parameters$log_lf_tau)))
  if (!is.null(parameters$log_wf_tau)) map$log_wf_tau <- factor(rep(NA, length(parameters$log_wf_tau)))
  map
}

synth_obj <- function(d, parameters = NULL, map = NULL) {
  if (is.null(parameters)) parameters <- synth_parameters(d)
  if (is.null(map)) map <- synth_map(parameters)
  suppressWarnings(RTMB::MakeADFun(
    func = cmb(opal_model, d), parameters = parameters, map = map, silent = TRUE
  ))
}