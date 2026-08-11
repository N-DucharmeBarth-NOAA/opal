#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(RTMB)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(purrr)
})

source_opal_r <- function(path = ".") {
  files <- list.files(file.path(path, "R"), pattern = "\\.R$", full.names = TRUE)
  for (file in files) sys.source(file, envir = .GlobalEnv)
  invisible(files)
}

normalize_scenario <- function(scenario) {
  scenario <- as.list(scenario)
  lapply(scenario, function(x) x[[1]])
}

add_length_bins <- function(data) {
  data$len_lower <- seq(
    from = data$len_bin_start,
    by = data$len_bin_width,
    length.out = data$n_len
  )
  data$len_upper <- data$len_lower + data$len_bin_width
  data$len_mid <- data$len_lower + data$len_bin_width / 2
  data
}

load_bet_inputs <- function() {
  load("data/wcpo_bet_data.rda")
  load("data/wcpo_bet_parameters.rda")
  data <- add_length_bins(wcpo_bet_data)
  list(data = data, parameters = wcpo_bet_parameters)
}

set_cpue_index <- function(data, cpue_index) {
  cpue_index <- match.arg(cpue_index, c("single", "quarter"))
  if (cpue_index == "single") {
    data$cpue_data$index <- rep(1L, nrow(data$cpue_data))
    data$n_index <- 1L
  } else {
    data$cpue_data$index <- as.integer(
      factor(data$cpue_data$month, levels = c(2, 5, 8, 11))
    )
    data$n_index <- 4L
  }
  data
}

make_par_sel <- function(data, wcpo_pars, selectivity_scale) {
  selectivity_scale <- match.arg(selectivity_scale, c("current_scaled", "raw_unconverted"))
  par_sel <- as.matrix(wcpo_pars$par_sel)
  if (selectivity_scale == "current_scaled") {
    double_normal_f <- data$sel_type_f == 2L
    par_sel[double_normal_f, 3:4] <- par_sel[double_normal_f, 3:4] - log(sd(data$len_mid))
  }
  par_sel
}

make_parameters <- function(data,
                            wcpo_pars,
                            selectivity_scale = "current_scaled",
                            tau_start = "current_0.1",
                            rdev_mode = "estimate",
                            init_dev_mode = "absent") {
  tau_start <- match.arg(tau_start, c("current_0.1", "published_zero"))
  rdev_mode <- match.arg(rdev_mode, c("estimate", "fixed_reference", "fixed_zero"))
  init_dev_mode <- match.arg(init_dev_mode, c("absent", "estimate"))

  log_cpue_tau <- if (tau_start == "current_0.1") {
    rep(log(0.1), data$n_index)
  } else {
    rep(as.numeric(wcpo_pars$log_cpue_tau), data$n_index)
  }

  rdev_y <- if (rdev_mode == "fixed_zero") {
    rep(0, data$n_year)
  } else {
    as.numeric(wcpo_pars$rdev_y)
  }

  parameters <- list(
    log_B0 = 20,
    log_h = as.numeric(wcpo_pars$log_h),
    log_sigma_r = as.numeric(wcpo_pars$log_sigma_r),
    log_cpue_q = rep(as.numeric(wcpo_pars$log_cpue_q), data$n_index),
    cpue_creep = rep(as.numeric(wcpo_pars$cpue_creep), data$n_index),
    log_cpue_tau = log_cpue_tau,
    log_cpue_omega = rep(as.numeric(wcpo_pars$log_cpue_omega), data$n_index),
    log_lf_tau = as.numeric(log(rep(0.1, data$n_fishery))),
    log_wf_tau = rep(0, data$n_fishery),
    log_L1 = as.numeric(wcpo_pars$log_L1),
    log_L2 = as.numeric(wcpo_pars$log_L2),
    log_k = as.numeric(wcpo_pars$log_k),
    log_CV1 = as.numeric(wcpo_pars$log_CV1),
    log_CV2 = as.numeric(wcpo_pars$log_CV2),
    par_sel = make_par_sel(data, wcpo_pars, selectivity_scale),
    rdev_y = rdev_y
  )

  if (init_dev_mode == "estimate") {
    parameters$init_rdev_a <- rep(0, data$n_age)
  }

  parameters
}

make_map <- function(data,
                     parameters,
                     rdev_mode = "estimate",
                     init_dev_mode = "absent",
                     estimate_q = TRUE) {
  rdev_mode <- match.arg(rdev_mode, c("estimate", "fixed_reference", "fixed_zero"))
  init_dev_mode <- match.arg(init_dev_mode, c("absent", "estimate"))

  map_sel <- matrix(NA, nrow(parameters$par_sel), ncol(parameters$par_sel))
  map_rdev <- if (rdev_mode == "estimate") {
    seq_along(parameters$rdev_y)
  } else {
    rep(NA, length(parameters$rdev_y))
  }

  map <- list(
    log_h = factor(NA),
    log_sigma_r = factor(NA),
    log_cpue_q = factor(if (estimate_q) seq_len(data$n_index) else rep(NA, data$n_index)),
    cpue_creep = factor(rep(NA, data$n_index)),
    log_cpue_tau = factor(rep(NA, data$n_index)),
    log_cpue_omega = factor(rep(NA, data$n_index)),
    log_lf_tau = factor(rep(NA, data$n_fishery)),
    log_wf_tau = factor(rep(NA, data$n_fishery)),
    log_L1 = factor(NA),
    log_L2 = factor(NA),
    log_k = factor(NA),
    log_CV1 = factor(NA),
    log_CV2 = factor(NA),
    par_sel = factor(map_sel),
    rdev_y = factor(map_rdev)
  )

  if ("init_rdev_a" %in% names(parameters)) {
    map$init_rdev_a <- factor(seq_len(data$n_age))
  }

  map
}

make_bounds <- function(obj) {
  lower <- rep(-Inf, length(obj$par))
  upper <- rep(Inf, length(obj$par))
  lower[names(obj$par) == "log_B0"] <- log(1)
  upper[names(obj$par) == "log_B0"] <- 22
  lower[names(obj$par) == "log_cpue_q"] <- log(0.1)
  upper[names(obj$par) == "log_cpue_q"] <- log(10)
  lower[names(obj$par) == "rdev_y"] <- -5
  upper[names(obj$par) == "rdev_y"] <- 5
  lower[names(obj$par) == "init_rdev_a"] <- -5
  upper[names(obj$par) == "init_rdev_a"] <- 5
  list(lower = lower, upper = upper)
}

prepare_data <- function(data,
                         use_priors = FALSE,
                         prior_log_B0_mean = NA_real_,
                         lf_switch = 0L,
                         wf_switch = 0L) {
  data <- add_length_bins(data)
  data$lf_switch <- as.integer(lf_switch)
  data$wf_switch <- as.integer(wf_switch)
  if (lf_switch == 0L) data$n_lf <- 0L
  if (wf_switch == 0L) data$n_wf <- 0L
  if (is.finite(prior_log_B0_mean)) data$prior_log_B0_mean <- prior_log_B0_mean
  data$priors <- if (use_priors) NULL else list()
  data
}

build_fit <- function(scenario, base_data, wcpo_pars) {
  scenario <- normalize_scenario(scenario)
  data <- base_data |>
    set_cpue_index(scenario$cpue_index) |>
    prepare_data(
      use_priors = scenario$use_priors,
      prior_log_B0_mean = scenario$prior_log_B0_mean,
      lf_switch = scenario$lf_switch,
      wf_switch = scenario$wf_switch
    )

  parameters <- make_parameters(
    data = data,
    wcpo_pars = wcpo_pars,
    selectivity_scale = scenario$selectivity_scale,
    tau_start = scenario$tau_start,
    rdev_mode = scenario$rdev_mode,
    init_dev_mode = scenario$init_dev_mode
  )

  if (isTRUE(scenario$use_priors)) {
    data$priors <- get_priors(parameters = parameters, data = data)
  }

  map <- make_map(
    data = data,
    parameters = parameters,
    rdev_mode = scenario$rdev_mode,
    init_dev_mode = scenario$init_dev_mode
  )

  obj <- MakeADFun(func = cmb(opal_model, data), parameters = parameters, map = map)
  obj$env$tracemgc <- FALSE
  bounds <- make_bounds(obj)
  list(data = data, parameters = parameters, map = map, obj = obj, bounds = bounds)
}

fit_scenario <- function(scenario, base_data, wcpo_pars,
                         eval_max = 3000L, iter_max = 3000L,
                         max_restarts = 3L) {
  scenario <- normalize_scenario(scenario)
  built <- build_fit(scenario, base_data, wcpo_pars)
  obj <- built$obj
  bounds <- built$bounds
  opt <- list(par = obj$par, convergence = NA_integer_, message = "not run")
  fit_grad <- Inf
  control <- list(eval.max = eval_max, iter.max = iter_max)

  for (i in seq_len(max_restarts)) {
    opt <- nlminb(
      start = opt$par,
      objective = obj$fn,
      gradient = obj$gr,
      hessian = obj$he,
      lower = bounds$lower,
      upper = bounds$upper,
      control = control
    )
    fit_grad <- max(abs(obj$gr(opt$par)))
    if (is.finite(fit_grad) && fit_grad < 1e-2) break
  }

  obj$env$last.par.best <- opt$par
  rep <- obj$report(opt$par)
  par_list <- obj$env$parList(opt$par)
  cpue_obs <- built$data$cpue_data$value
  cpue_pred <- as.numeric(rep$cpue_pred)
  log_resid <- log(cpue_obs) - log(cpue_pred)
  sb <- as.numeric(rep$spawning_biomass_y)
  early_n <- min(8L, length(cpue_obs))
  early_sb_end <- min(9L, length(sb))
  scenario_name <- scenario$name
  scenario_cpue_index <- scenario$cpue_index
  scenario_selectivity_scale <- scenario$selectivity_scale
  scenario_tau_start <- scenario$tau_start
  scenario_rdev_mode <- scenario$rdev_mode
  scenario_init_dev_mode <- scenario$init_dev_mode
  scenario_use_priors <- scenario$use_priors
  scenario_prior_log_B0_mean <- scenario$prior_log_B0_mean
  scenario_lf_switch <- scenario$lf_switch
  scenario_wf_switch <- scenario$wf_switch

  metrics <- tibble(
    scenario = scenario_name,
    cpue_index = scenario_cpue_index,
    selectivity_scale = scenario_selectivity_scale,
    tau_start = scenario_tau_start,
    rdev_mode = scenario_rdev_mode,
    init_dev_mode = scenario_init_dev_mode,
    use_priors = scenario_use_priors,
    prior_log_B0_mean = scenario_prior_log_B0_mean,
    lf_switch = scenario_lf_switch,
    wf_switch = scenario_wf_switch,
    n_par = length(obj$par),
    active_parameters = paste(unique(names(obj$par)), collapse = ","),
    convergence = opt$convergence,
    message = opt$message,
    nll = obj$fn(opt$par),
    max_gradient = fit_grad,
    lp_prior = as.numeric(rep$lp_prior),
    lp_penalty = as.numeric(rep$lp_penalty),
    lp_rec = as.numeric(rep$lp_rec),
    lp_init_rec = if ("lp_init_rec" %in% names(rep)) as.numeric(rep$lp_init_rec) else NA_real_,
    lp_cpue = sum(as.numeric(rep$lp_cpue)),
    lp_lf = sum(as.numeric(rep$lp_lf)),
    lp_wf = sum(as.numeric(rep$lp_wf)),
    log_cpue_rmse = sqrt(mean(log_resid^2)),
    early_log_cpue_rmse = sqrt(mean(log_resid[seq_len(early_n)]^2)),
    early_log_cpue_bias = mean(log_resid[seq_len(early_n)]),
    B0 = as.numeric(rep$B0),
    first_spawning_biomass = sb[1],
    early_spawning_biomass_ratio = sb[early_sb_end] / sb[1],
    early_spawning_biomass_max_ratio = max(sb[seq_len(early_sb_end)]) / sb[1],
    final_spawning_biomass_ratio = tail(sb, 1) / as.numeric(rep$B0),
    q_est = paste(round(as.numeric(par_list$log_cpue_q), 4), collapse = ","),
    first_8_rdev = paste(round(as.numeric(par_list$rdev_y[seq_len(early_n)]), 4), collapse = ",")
  )

  cpue_df <- built$data$cpue_data |>
    mutate(
      scenario = scenario_name,
      obs_id = row_number(),
      pred = cpue_pred,
      log_resid = log_resid
    )

  sb_df <- tibble(
    scenario = scenario_name,
    time_step = seq_along(sb),
    spawning_biomass = sb,
    relative_to_initial = sb / sb[1],
    relative_to_B0 = sb / as.numeric(rep$B0)
  )

  list(
    scenario = scenario,
    metrics = metrics,
    cpue = cpue_df,
    spawning_biomass = sb_df,
    fit = built,
    opt = opt,
    report = rep
  )
}

make_scenarios <- function() {
  tibble::tribble(
    ~name, ~cpue_index, ~selectivity_scale, ~tau_start, ~rdev_mode, ~init_dev_mode, ~use_priors, ~lf_switch, ~wf_switch,
    "published_single_rdev_est_tau0_noprior", "single", "current_scaled", "published_zero", "estimate", "absent", FALSE, 0L, 0L,
    "single_rdev_est_tau01_noprior", "single", "current_scaled", "current_0.1", "estimate", "absent", FALSE, 0L, 0L,
    "quarter_rdev_est_tau01_noprior", "quarter", "current_scaled", "current_0.1", "estimate", "absent", FALSE, 0L, 0L,
    "quarter_rdev_est_tau01_prior", "quarter", "current_scaled", "current_0.1", "estimate", "absent", TRUE, 0L, 0L,
    "quarter_rdev_fixed_ref_tau01_noprior", "quarter", "current_scaled", "current_0.1", "fixed_reference", "absent", FALSE, 0L, 0L,
    "quarter_rdev_fixed_zero_tau01_noprior", "quarter", "current_scaled", "current_0.1", "fixed_zero", "absent", FALSE, 0L, 0L,
    "quarter_initdev_est_rdev_fixed_zero", "quarter", "current_scaled", "current_0.1", "fixed_zero", "estimate", FALSE, 0L, 0L,
    "quarter_initdev_est_rdev_est", "quarter", "current_scaled", "current_0.1", "estimate", "estimate", FALSE, 0L, 0L,
    "quarter_raw_sel_rdev_est_tau01_noprior", "quarter", "raw_unconverted", "current_0.1", "estimate", "absent", FALSE, 0L, 0L
  ) |>
    mutate(prior_log_B0_mean = NA_real_) |>
    bind_rows(tibble::tibble(
      name = "quarter_rdev_est_tau01_prior_B0_14",
      cpue_index = "quarter",
      selectivity_scale = "current_scaled",
      tau_start = "current_0.1",
      rdev_mode = "estimate",
      init_dev_mode = "absent",
      use_priors = TRUE,
      lf_switch = 0L,
      wf_switch = 0L,
      prior_log_B0_mean = log(1.7e6)
    ))
}

write_outputs <- function(results, out_dir) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  metrics <- bind_rows(map(results, "metrics"))
  cpue <- bind_rows(map(results, "cpue"))
  sb <- bind_rows(map(results, "spawning_biomass"))

  write.csv(metrics, file.path(out_dir, "scenario_metrics.csv"), row.names = FALSE)
  write.csv(cpue, file.path(out_dir, "scenario_cpue_predictions.csv"), row.names = FALSE)
  write.csv(sb, file.path(out_dir, "scenario_spawning_biomass.csv"), row.names = FALSE)
  saveRDS(results, file.path(out_dir, "scenario_results.rds"))

  p_cpue <- ggplot(cpue, aes(x = obs_id)) +
    geom_point(aes(y = value), size = 0.8, alpha = 0.7) +
    geom_line(aes(y = pred), colour = "red3", linewidth = 0.5) +
    facet_wrap(~ scenario, scales = "free_y", ncol = 2) +
    labs(x = "CPUE observation", y = "CPUE")
  ggsave(file.path(out_dir, "cpue_fits.png"), p_cpue, width = 12, height = 14, dpi = 150)

  p_cpue_early <- cpue |>
    filter(obs_id <= 32) |>
    ggplot(aes(x = obs_id)) +
    geom_point(aes(y = value), size = 1.1, alpha = 0.8) +
    geom_line(aes(y = pred), colour = "red3", linewidth = 0.6) +
    facet_wrap(~ scenario, scales = "free_y", ncol = 2) +
    labs(x = "CPUE observation", y = "CPUE", title = "Early CPUE fit")
  ggsave(file.path(out_dir, "cpue_fits_first32.png"), p_cpue_early, width = 12, height = 14, dpi = 150)

  p_sb <- ggplot(sb, aes(x = time_step, y = relative_to_initial)) +
    geom_hline(yintercept = 1, colour = "grey60", linewidth = 0.3) +
    geom_line(colour = "steelblue4", linewidth = 0.6) +
    facet_wrap(~ scenario, scales = "free_y", ncol = 2) +
    labs(x = "Model time step", y = "Spawning biomass / initial")
  ggsave(file.path(out_dir, "spawning_biomass_relative_initial.png"), p_sb, width = 12, height = 14, dpi = 150)

  invisible(list(metrics = metrics, cpue = cpue, spawning_biomass = sb))
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out_dir <- if (length(args) >= 1L) args[[1]] else "inst/scripts/bet_cpue_diagnostics"
  eval_max <- as.integer(Sys.getenv("OPAL_BET_EVAL_MAX", "500"))
  iter_max <- as.integer(Sys.getenv("OPAL_BET_ITER_MAX", as.character(eval_max)))
  max_restarts <- as.integer(Sys.getenv("OPAL_BET_MAX_RESTARTS", "1"))
  scenario_filter <- Sys.getenv("OPAL_BET_SCENARIOS", "")
  source_opal_r()
  inputs <- load_bet_inputs()
  scenarios <- make_scenarios()
  if (nzchar(scenario_filter)) {
    keep <- trimws(strsplit(scenario_filter, ",", fixed = TRUE)[[1]])
    scenarios <- scenarios |>
      filter(name %in% keep | as.character(row_number()) %in% keep)
  }

  cat(
    "Running ", nrow(scenarios), " BET CPUE diagnostic scenarios",
    " (eval_max=", eval_max, ", iter_max=", iter_max,
    ", max_restarts=", max_restarts, ")\n",
    sep = ""
  )
  flush.console()
  results <- vector("list", nrow(scenarios))
  for (i in seq_len(nrow(scenarios))) {
    scenario <- scenarios[i, ]
    cat("[", i, "/", nrow(scenarios), "] ", scenario$name, "\n", sep = "")
    flush.console()
    results[[i]] <- fit_scenario(
      scenario,
      inputs$data,
      inputs$parameters,
      eval_max = eval_max,
      iter_max = iter_max,
      max_restarts = max_restarts
    )
  }

  outputs <- write_outputs(results, out_dir)
  print(outputs$metrics |>
    arrange(early_log_cpue_rmse) |>
    select(
      scenario, n_par, convergence, max_gradient, nll,
      lp_cpue, lp_prior, lp_penalty, log_cpue_rmse,
      early_log_cpue_rmse, early_spawning_biomass_ratio,
      early_spawning_biomass_max_ratio, q_est, first_8_rdev
    ))
}

invisible(NULL)
