make_opal_fit_data <- function() {
  list(
    n_year = 2L,
    n_season = 1L,
    n_age = 5L,
    n_fishery = 2L,
    n_len = 15L,
    min_age = 1L,
    max_age = 5L,
    first_yr = 1L,
    first_yr_catch = 1L,
    last_yr = 2L,
    len_bin_start = 20,
    len_bin_width = 2,
    len_lower = seq(20, by = 2, length.out = 15),
    len_upper = seq(22, by = 2, length.out = 15),
    len_mid = seq(21, by = 2, length.out = 15),
    A1 = 1L,
    A2 = 5L,
    M = rep(0.3, 5),
    maturity = c(0, 0.2, 0.5, 0.8, 1),
    fecundity = c(0, 100, 500, 1000, 1500),
    lw_a = 0.00001,
    lw_b = 3,
    catch_obs_ysf = array(c(100, 200, 150, 180), dim = c(2, 1, 2)),
    catch_units_f = c(1L, 1L),
    removal_switch_f = c(0L, 0L),
    sel_type_f = c(1L, 1L),
    cpue_switch = 1L,
    cpue_data = data.frame(
      ts = c(1L, 2L),
      fishery = c(1L, 1L),
      value = c(0.5, 0.48),
      se = c(0.1, 0.1),
      units = c(1L, 1L)
    ),
    lf_switch = 0L,
    n_lf = 0L,
    wf_switch = 0L,
    n_wf = 0L,
    log_L1 = log(30),
    log_L2 = log(60),
    log_k = log(0.2)
  )
}

make_opal_fit_parameters <- function(data) {
  list(
    log_B0 = 20,
    log_h = 0.7,
    log_sigma_r = log(0.6),
    log_cpue_q = 0,
    cpue_creep = 0,
    log_cpue_tau = log(0.2),
    log_cpue_omega = log(0.1),
    log_lf_tau = rep(log(0.1), data$n_fishery),
    log_wf_tau = rep(0, data$n_fishery),
    log_L1 = log(30),
    log_L2 = log(60),
    log_k = log(0.2),
    log_CV1 = log(0.1),
    log_CV2 = log(0.05),
    par_sel = matrix(
      c(40, 5, 0, 0, 0, 0, 45, 4, 0, 0, 0, 0),
      nrow = data$n_fishery,
      byrow = TRUE
    ),
    rdev_y = rep(0, data$n_year)
  )
}

make_opal_fit_map <- function(parameters) {
  list(
    log_h = factor(NA),
    log_sigma_r = factor(NA),
    log_cpue_q = factor(NA),
    cpue_creep = factor(NA),
    log_cpue_tau = factor(NA),
    log_cpue_omega = factor(NA),
    log_lf_tau = factor(rep(NA, length(parameters$log_lf_tau))),
    log_wf_tau = factor(rep(NA, length(parameters$log_wf_tau))),
    log_L1 = factor(NA),
    log_L2 = factor(NA),
    log_k = factor(NA),
    log_CV1 = factor(NA),
    log_CV2 = factor(NA),
    par_sel = factor(matrix(
      NA,
      nrow = nrow(parameters$par_sel),
      ncol = ncol(parameters$par_sel)
    )),
    rdev_y = factor(rep(NA, length(parameters$rdev_y)))
  )
}

make_opal_fit_fixture <- local({
  cached <- NULL

  function(mcmc = TRUE) {
    if (is.null(cached)) {
      data <- make_opal_fit_data()
      parameters <- make_opal_fit_parameters(data)
      map <- make_opal_fit_map(parameters)
      obj <- suppressWarnings(RTMB::MakeADFun(
        func = cmb(opal_model, data),
        parameters = parameters,
        map = map,
        silent = TRUE
      ))
      opt <- list(
        par = obj$par,
        objective = obj$fn(obj$par),
        convergence = 0L,
        message = "test fixture"
      )
      cached <<- list(data = data, obj = obj, opt = opt)
    }

    fixture <- cached
    posterior <- NULL
    if (isTRUE(mcmc)) {
      active_names <- opal:::.opal_expand_parameter_names(
        names(fixture$obj$par)
      )
      posterior <- rbind(
        c(unname(fixture$obj$par), lp__ = -fixture$opt$objective),
        c(unname(fixture$obj$par) + 0.01, lp__ = -fixture$opt$objective)
      )
      colnames(posterior) <- c(active_names, "lp__")
    }

    opal_fit(
      data = fixture$data,
      obj = fixture$obj,
      opt = fixture$opt,
      bounds = list(
        lower = stats::setNames(-Inf, names(fixture$opt$par)),
        upper = stats::setNames(Inf, names(fixture$opt$par))
      ),
      diagnostics = list(max_gradient = 0),
      metadata = list(stock = "synthetic"),
      mcmc = posterior,
      derived = list(projection = data.frame(year = 2025L, biomass = 100))
    )
  }
})
