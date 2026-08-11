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

num_env <- function(name, default) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || !nzchar(value)) return(default)
  as.numeric(value)
}

print_progress <- function(...) {
  cat(...)
  flush(stdout())
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
  load("data/wcpo_bet_lf.rda")
  list(
    data = add_length_bins(wcpo_bet_data),
    parameters = wcpo_bet_parameters,
    lf = wcpo_bet_lf
  )
}

set_quarter_cpue <- function(data) {
  data$cpue_data$index <- as.integer(
    factor(data$cpue_data$month, levels = c(2, 5, 8, 11))
  )
  data$n_index <- 4L
  data
}

make_lf_wide <- function(lf) {
  lf |>
    pivot_wider(
      id_cols = c(fishery, year, month, ts),
      names_from = bin,
      values_from = value,
      values_fill = 0
    ) |>
    arrange(fishery, ts)
}

prepare_lf_data <- function(data, lf_wide, target_fishery,
                            lf_var_adjust_scalar = 80) {
  lf_var_adjust <- rep(1, data$n_fishery)
  lf_var_adjust[] <- lf_var_adjust_scalar
  data <- prep_lf_data(
    data = data,
    lf_wide = lf_wide,
    lf_keep_fisheries = target_fishery,
    lf_var_adjust = lf_var_adjust,
    lf_switch = 1L
  )
  data$lf_switch <- 1L
  data$wf_switch <- 0L
  data$n_wf <- 0L
  data$priors <- list()
  data
}

make_par_sel <- function(data, wcpo_pars) {
  par_sel <- as.matrix(wcpo_pars$par_sel)
  double_normal_f <- data$sel_type_f == 2L
  par_sel[double_normal_f, 3:4] <- par_sel[double_normal_f, 3:4] - log(sd(data$len_mid))
  par_sel
}

make_parameters <- function(data, wcpo_pars,
                            log_B0_start = num_env("OPAL_BET_LF_LOG_B0_START", 15)) {
  list(
    log_B0 = log_B0_start,
    log_h = as.numeric(wcpo_pars$log_h),
    log_sigma_r = as.numeric(wcpo_pars$log_sigma_r),
    log_cpue_q = rep(0, data$n_index),
    cpue_creep = rep(as.numeric(wcpo_pars$cpue_creep), data$n_index),
    log_cpue_tau = rep(log(0.1), data$n_index),
    log_cpue_omega = rep(as.numeric(wcpo_pars$log_cpue_omega), data$n_index),
    log_lf_tau = as.numeric(log(rep(0.1, data$n_fishery))),
    log_wf_tau = rep(0, data$n_fishery),
    log_L1 = as.numeric(wcpo_pars$log_L1),
    log_L2 = as.numeric(wcpo_pars$log_L2),
    log_k = as.numeric(wcpo_pars$log_k),
    log_CV1 = as.numeric(wcpo_pars$log_CV1),
    log_CV2 = as.numeric(wcpo_pars$log_CV2),
    par_sel = make_par_sel(data, wcpo_pars),
    rdev_y = as.numeric(wcpo_pars$rdev_y)
  )
}

selectivity_columns <- function(data, fishery) {
  if (data$sel_type_f[fishery] == 1L) c(1L, 2L) else c(1L, 3L, 4L)
}

selectivity_natural <- function(data, par, fishery) {
  mu_len <- mean(data$len_mid)
  sd_len <- sd(data$len_mid)
  if (data$sel_type_f[fishery] == 1L) {
    values <- c(
      inflection_cm = mu_len + par[1] * sd_len,
      width95_cm = exp(par[2]) * sd_len
    )
  } else {
    values <- c(
      peak_cm = mu_len + par[1] * sd_len,
      asc_denom_cm2 = exp(par[3]) * sd_len^2,
      desc_denom_cm2 = exp(par[4]) * sd_len^2
    )
  }
  paste(paste(names(values), round(as.numeric(values), 4), sep = "="), collapse = ",")
}

make_map <- function(data, parameters, target_fishery) {
  map_sel <- matrix(NA_integer_, nrow(parameters$par_sel), ncol(parameters$par_sel))
  active_sel_cols <- selectivity_columns(data, target_fishery)
  map_sel[target_fishery, active_sel_cols] <- seq_along(active_sel_cols)

  map <- list(
    log_h = factor(NA),
    log_sigma_r = factor(NA),
    log_cpue_q = factor(seq_len(data$n_index)),
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
    rdev_y = factor(seq_along(parameters$rdev_y))
  )

  list(map = map, map_sel = map_sel, active_sel_cols = active_sel_cols)
}

make_bounds <- function(obj, data, map, map_sel,
                        log_B0_lower = num_env("OPAL_BET_LF_LOG_B0_LOWER", 14),
                        log_B0_upper = num_env("OPAL_BET_LF_LOG_B0_UPPER", 22),
                        logistic_width_lower = num_env("OPAL_BET_LF_LOGISTIC_WIDTH_LOWER", 1e-6),
                        logistic_width_upper = num_env("OPAL_BET_LF_LOGISTIC_WIDTH_UPPER", 500),
                        width_denom_lower = num_env("OPAL_BET_LF_WIDTH_DENOM_LOWER", 4),
                        width_denom_upper = num_env("OPAL_BET_LF_WIDTH_DENOM_UPPER", exp(8))) {
  lower <- rep(-Inf, length(obj$par))
  upper <- rep(Inf, length(obj$par))

  lower[names(obj$par) == "log_B0"] <- log_B0_lower
  upper[names(obj$par) == "log_B0"] <- log_B0_upper
  lower[names(obj$par) == "log_cpue_q"] <- log(0.1)
  upper[names(obj$par) == "log_cpue_q"] <- log(10)
  lower[names(obj$par) == "rdev_y"] <- -5
  upper[names(obj$par) == "rdev_y"] <- 5

  mu_len <- mean(data$len_mid)
  sd_len <- sd(data$len_mid)
  par_sel_pos <- which(names(obj$par) == "par_sel")
  par_sel_level_pos <- setNames(par_sel_pos, levels(map$par_sel))

  for (f in seq_len(data$n_fishery)) {
    idx <- rep(NA_integer_, ncol(map_sel))
    for (p in seq_len(ncol(map_sel))) {
      level <- as.character(map_sel[f, p])
      if (!is.na(level)) idx[p] <- par_sel_level_pos[[level]]
    }

    if (data$sel_type_f[f] == 1L) {
      if (!is.na(idx[1])) {
        lower[idx[1]] <- (5 - mu_len) / sd_len
        upper[idx[1]] <- (200 - mu_len) / sd_len
      }
      if (!is.na(idx[2])) {
        lower[idx[2]] <- log(logistic_width_lower / sd_len)
        upper[idx[2]] <- log(logistic_width_upper / sd_len)
      }
    } else {
      if (!is.na(idx[1])) {
        lower[idx[1]] <- (10.1 - mu_len) / sd_len
        upper[idx[1]] <- (200 - mu_len) / sd_len
      }
      if (!is.na(idx[2])) {
        lower[idx[2]] <- -7
        upper[idx[2]] <- 7
      }
      shift <- 2 * log(sd_len)
      width_log_lower <- log(width_denom_lower)
      width_log_upper <- log(width_denom_upper)
      if (!is.na(idx[3])) {
        lower[idx[3]] <- width_log_lower - shift
        upper[idx[3]] <- width_log_upper - shift
      }
      if (!is.na(idx[4])) {
        lower[idx[4]] <- width_log_lower - shift
        upper[idx[4]] <- width_log_upper - shift
      }
      if (!is.na(idx[5])) {
        lower[idx[5]] <- -9
        upper[idx[5]] <- 9
      }
      if (!is.na(idx[6])) {
        lower[idx[6]] <- -9
        upper[idx[6]] <- 9
      }
    }
  }

  list(lower = lower, upper = upper)
}

build_fit <- function(base_data, wcpo_pars, lf_wide, target_fishery) {
  data <- base_data |>
    add_length_bins() |>
    set_quarter_cpue() |>
    prepare_lf_data(lf_wide = lf_wide, target_fishery = target_fishery)

  parameters <- make_parameters(data, wcpo_pars)
  map_info <- make_map(data, parameters, target_fishery)

  obj <- MakeADFun(
    func = cmb(opal_model, data),
    parameters = parameters,
    map = map_info$map
  )
  obj$env$tracemgc <- FALSE
  bounds <- make_bounds(obj, data, map_info$map, map_info$map_sel)

  list(
    data = data,
    parameters = parameters,
    map = map_info$map,
    map_sel = map_info$map_sel,
    active_sel_cols = map_info$active_sel_cols,
    obj = obj,
    bounds = bounds
  )
}

fit_one_lf <- function(base_data, wcpo_pars, lf_wide, target_fishery,
                       eval_max = 1000L, iter_max = 1000L,
                       max_restarts = 2L, trace = FALSE) {
  built <- build_fit(base_data, wcpo_pars, lf_wide, target_fishery)
  obj <- built$obj
  opt <- list(par = obj$par, convergence = NA_integer_, message = "not run")
  fit_grad <- Inf
  control <- list(eval.max = eval_max, iter.max = iter_max)

  if (isTRUE(trace)) {
    start_nll <- tryCatch(obj$fn(obj$par), error = function(e) NaN)
    print_progress("  start nll=", signif(start_nll, 8), "\n")
  }

  for (i in seq_len(max_restarts)) {
    opt <- nlminb(
      start = opt$par,
      objective = obj$fn,
      gradient = obj$gr,
      lower = built$bounds$lower,
      upper = built$bounds$upper,
      control = control
    )
    fit_grad <- max(abs(obj$gr(opt$par)))
    if (isTRUE(trace)) {
      print_progress(
        "  restart ", i,
        ": convergence=", opt$convergence,
        ", nll=", signif(opt$objective, 8),
        ", max_grad=", signif(fit_grad, 6),
        "\n"
      )
    }
    if (is.finite(fit_grad) && fit_grad < 1e-2) break
  }

  obj$env$last.par.best <- opt$par
  rep <- obj$report(opt$par)
  par_list <- obj$env$parList(opt$par)
  metrics <- make_metrics(built, opt, fit_grad, rep, par_list, target_fishery)
  list(
    target_fishery = target_fishery,
    fit = built,
    opt = opt,
    report = rep,
    par_list = par_list,
    metrics = metrics,
    cpue = make_cpue_df(built$data, rep, target_fishery),
    lf = make_lf_df(built$data, rep, target_fishery),
    selectivity = make_selectivity_df(built$data, built$parameters, par_list, target_fishery),
    spawning_biomass = make_spawning_biomass_df(rep, target_fishery)
  )
}

make_cpue_df <- function(data, rep, target_fishery) {
  cpue_pred <- as.numeric(rep$cpue_pred)
  data$cpue_data |>
    mutate(
      target_lf_fishery = target_fishery,
      obs_id = row_number(),
      pred = cpue_pred,
      log_resid = log(value) - log(pred)
    )
}

make_lf_df <- function(data, rep, target_fishery) {
  lf_pred <- rep$lf_pred[[1]]
  obs_counts <- data$lf_obs_data$obs[[1]]
  obs_prop <- obs_counts / rowSums(obs_counts)
  bmin <- data$lf_minbin[target_fishery]
  bmax <- data$lf_maxbin[target_fishery]
  len <- data$len_mid[bmin:bmax]

  pred_df <- as.data.frame(lf_pred)
  names(pred_df) <- len
  obs_df <- as.data.frame(obs_prop)
  names(obs_df) <- len

  bind_rows(
    obs_df |>
      mutate(obs_id = row_number(), year = data$lf_year, type = "observed"),
    pred_df |>
      mutate(obs_id = row_number(), year = data$lf_year, type = "predicted")
  ) |>
    pivot_longer(
      cols = -c(obs_id, year, type),
      names_to = "length",
      values_to = "proportion"
    ) |>
    mutate(
      target_lf_fishery = target_fishery,
      fishery = target_fishery,
      length = as.numeric(length)
    )
}

make_selectivity_df <- function(data, parameters, par_list, target_fishery) {
  par_init <- parameters$par_sel[target_fishery, ]
  par_est <- par_list$par_sel[target_fishery, ]
  sel_init <- if (data$sel_type_f[target_fishery] == 1L) {
    sel_logistic(data$len_mid, par_init)
  } else {
    sel_double_normal(data$len_mid, par_init)
  }
  sel_est <- if (data$sel_type_f[target_fishery] == 1L) {
    sel_logistic(data$len_mid, par_est)
  } else {
    sel_double_normal(data$len_mid, par_est)
  }
  tibble(
    target_lf_fishery = target_fishery,
    fishery = target_fishery,
    length = rep(data$len_mid, 2),
    selectivity = c(as.numeric(sel_init), as.numeric(sel_est)),
    type = rep(c("initial", "estimated"), each = length(data$len_mid))
  )
}

make_spawning_biomass_df <- function(rep, target_fishery) {
  sb <- as.numeric(rep$spawning_biomass_y)
  tibble(
    target_lf_fishery = target_fishery,
    time_step = seq_along(sb),
    spawning_biomass = sb,
    relative_to_initial = sb / sb[1],
    relative_to_B0 = sb / as.numeric(rep$B0)
  )
}

make_metrics <- function(built, opt, fit_grad, rep, par_list, target_fishery) {
  cpue_obs <- built$data$cpue_data$value
  cpue_pred <- as.numeric(rep$cpue_pred)
  cpue_log_resid <- log(cpue_obs) - log(cpue_pred)
  early_n <- min(8L, length(cpue_obs))

  lf_df <- make_lf_df(built$data, rep, target_fishery)
  lf_wide <- lf_df |>
    select(target_lf_fishery, obs_id, year, length, type, proportion) |>
    pivot_wider(names_from = type, values_from = proportion)

  mean_len <- lf_wide |>
    group_by(obs_id, year) |>
    summarise(
      obs_mean_length = sum(length * observed),
      pred_mean_length = sum(length * predicted),
      .groups = "drop"
    )

  sel_active <- par_list$par_sel[target_fishery, built$active_sel_cols]
  sel_start <- built$parameters$par_sel[target_fishery, built$active_sel_cols]
  sel_start_nat <- selectivity_natural(
    built$data,
    built$parameters$par_sel[target_fishery, ],
    target_fishery
  )
  sel_est_nat <- selectivity_natural(
    built$data,
    par_list$par_sel[target_fishery, ],
    target_fishery
  )
  sb <- as.numeric(rep$spawning_biomass_y)

  tibble(
    target_lf_fishery = target_fishery,
    sel_type = if (built$data$sel_type_f[target_fishery] == 1L) "logistic" else "double_normal",
    active_sel_cols = paste(built$active_sel_cols, collapse = ","),
    n_par = length(built$obj$par),
    n_lf_obs = built$data$n_lf,
    mean_lf_n = mean(built$data$lf_n),
    convergence = opt$convergence,
    message = opt$message,
    nll = built$obj$fn(opt$par),
    max_gradient = fit_grad,
    B0 = as.numeric(rep$B0),
    lp_prior = as.numeric(rep$lp_prior),
    lp_penalty = as.numeric(rep$lp_penalty),
    lp_rec = as.numeric(rep$lp_rec),
    lp_cpue = sum(as.numeric(rep$lp_cpue)),
    lp_lf = sum(as.numeric(rep$lp_lf)),
    lp_wf = sum(as.numeric(rep$lp_wf)),
    log_cpue_rmse = sqrt(mean(cpue_log_resid^2)),
    early_log_cpue_rmse = sqrt(mean(cpue_log_resid[seq_len(early_n)]^2)),
    lf_prop_rmse = sqrt(mean((lf_wide$observed - lf_wide$predicted)^2)),
    lf_prop_mae = mean(abs(lf_wide$observed - lf_wide$predicted)),
    lf_mean_length_rmse = sqrt(mean((mean_len$obs_mean_length - mean_len$pred_mean_length)^2)),
    lf_mean_length_bias = mean(mean_len$obs_mean_length - mean_len$pred_mean_length),
    early_spawning_biomass_ratio = sb[min(9L, length(sb))] / sb[1],
    final_spawning_biomass_ratio = tail(sb, 1) / as.numeric(rep$B0),
    q_est = paste(round(as.numeric(par_list$log_cpue_q), 4), collapse = ","),
    sel_start = paste(round(as.numeric(sel_start), 4), collapse = ","),
    sel_est = paste(round(as.numeric(sel_active), 4), collapse = ","),
    sel_start_nat = sel_start_nat,
    sel_est_nat = sel_est_nat
  )
}

write_outputs <- function(results, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  metrics <- bind_rows(map(results, "metrics"))
  cpue <- bind_rows(map(results, "cpue"))
  lf <- bind_rows(map(results, "lf"))
  sel <- bind_rows(map(results, "selectivity"))
  sb <- bind_rows(map(results, "spawning_biomass"))

  write.csv(metrics, file.path(out_dir, "lf_metrics.csv"), row.names = FALSE)
  write.csv(cpue, file.path(out_dir, "lf_cpue_predictions.csv"), row.names = FALSE)
  write.csv(lf, file.path(out_dir, "lf_predictions.csv"), row.names = FALSE)
  write.csv(sel, file.path(out_dir, "lf_selectivity.csv"), row.names = FALSE)
  write.csv(sb, file.path(out_dir, "lf_spawning_biomass.csv"), row.names = FALSE)
  saveRDS(results, file.path(out_dir, "lf_results.rds"))

  plot_cpue(cpue, out_dir)
  plot_lf(lf, out_dir)
  plot_mean_length(lf, out_dir)
  plot_selectivity(sel, out_dir)
  plot_spawning_biomass(sb, out_dir)

  invisible(list(metrics = metrics, cpue = cpue, lf = lf, selectivity = sel, spawning_biomass = sb))
}

plot_cpue <- function(cpue, out_dir) {
  p <- cpue |>
    filter(obs_id <= 32) |>
    ggplot(aes(x = obs_id)) +
    geom_point(aes(y = value), size = 0.9, alpha = 0.8) +
    geom_line(aes(y = pred), colour = "red3", linewidth = 0.5) +
    facet_wrap(~ target_lf_fishery, scales = "free_y", ncol = 2) +
    labs(x = "CPUE observation", y = "CPUE", title = "Early CPUE fit by one-LF run")
  ggsave(file.path(out_dir, "lf_cpue_first32.png"), p, width = 10, height = 10, dpi = 150)
}

plot_lf <- function(lf, out_dir) {
  for (f in sort(unique(lf$target_lf_fishery))) {
    plot_year_data <- lf |>
      filter(target_lf_fishery == f) |>
      distinct(obs_id, year) |>
      arrange(year)
    plot_years <- tail(plot_year_data$year, min(24L, nrow(plot_year_data)))

    plot_data <- lf |>
      filter(target_lf_fishery == f, year %in% plot_years)
    obs_data <- plot_data |>
      filter(type == "observed")
    pred_data <- plot_data |>
      filter(type == "predicted")

    p <- plot_data |>
      ggplot(aes(x = length, y = proportion)) +
      geom_col(
        data = obs_data,
        fill = "grey75",
        width = 2
      ) +
      geom_line(
        data = pred_data,
        colour = "red3",
        linewidth = 0.6
      ) +
      facet_wrap(~ year, scales = "free_y", ncol = 6) +
      labs(x = "Length (cm)", y = "Proportion", title = paste("LF fishery", f)) +
      theme(strip.text = element_text(size = 7), axis.text = element_text(size = 6))
    ggsave(file.path(out_dir, paste0("lf_fit_f", f, ".png")), p, width = 12, height = 8, dpi = 150)
  }
}

plot_mean_length <- function(lf, out_dir) {
  mean_len <- lf |>
    group_by(target_lf_fishery, obs_id, year, type) |>
    summarise(mean_length = sum(length * proportion), .groups = "drop")

  p <- mean_len |>
    ggplot(aes(x = year, y = mean_length, colour = type)) +
    geom_line(linewidth = 0.5) +
    geom_point(size = 0.7) +
    facet_wrap(~ target_lf_fishery, scales = "free_y", ncol = 2) +
    scale_colour_manual(values = c(observed = "black", predicted = "red3")) +
    labs(x = "Model time step", y = "Mean length (cm)", colour = NULL)
  ggsave(file.path(out_dir, "lf_mean_length.png"), p, width = 10, height = 10, dpi = 150)
}

plot_selectivity <- function(sel, out_dir) {
  p <- sel |>
    ggplot(aes(x = length, y = selectivity, colour = type, linetype = type)) +
    geom_line(linewidth = 0.7) +
    facet_wrap(~ target_lf_fishery, ncol = 2) +
    scale_colour_manual(values = c(initial = "grey35", estimated = "red3")) +
    labs(x = "Length (cm)", y = "Selectivity", colour = NULL, linetype = NULL)
  ggsave(file.path(out_dir, "lf_selectivity.png"), p, width = 10, height = 10, dpi = 150)
}

plot_spawning_biomass <- function(sb, out_dir) {
  p <- sb |>
    ggplot(aes(x = time_step, y = relative_to_initial)) +
    geom_hline(yintercept = 1, colour = "grey60", linewidth = 0.3) +
    geom_line(colour = "steelblue4", linewidth = 0.6) +
    facet_wrap(~ target_lf_fishery, scales = "free_y", ncol = 2) +
    labs(x = "Model time step", y = "Spawning biomass / initial")
  ggsave(file.path(out_dir, "lf_spawning_biomass_relative_initial.png"), p, width = 10, height = 10, dpi = 150)
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out_dir <- if (length(args) >= 1L) args[[1]] else "inst/scripts/bet_lf_selectivity_diagnostics"
  eval_max <- as.integer(Sys.getenv("OPAL_BET_LF_EVAL_MAX", "1000"))
  iter_max <- as.integer(Sys.getenv("OPAL_BET_LF_ITER_MAX", as.character(eval_max)))
  max_restarts <- as.integer(Sys.getenv("OPAL_BET_LF_MAX_RESTARTS", "2"))
  fisheries <- Sys.getenv("OPAL_BET_LF_FISHERIES", "8,9,10,11,12,13,14")
  fisheries <- as.integer(trimws(strsplit(fisheries, ",", fixed = TRUE)[[1]]))

  source_opal_r()
  inputs <- load_bet_inputs()
  lf_wide <- make_lf_wide(inputs$lf)

  cat(
    "Running one-LF BET selectivity diagnostics for fisheries ",
    paste(fisheries, collapse = ", "),
    " (eval_max=", eval_max,
    ", iter_max=", iter_max,
    ", max_restarts=", max_restarts,
    ", log_B0_start=", num_env("OPAL_BET_LF_LOG_B0_START", 15),
    ", log_B0_lower=", num_env("OPAL_BET_LF_LOG_B0_LOWER", 14),
    ", width_denom_lower=", num_env("OPAL_BET_LF_WIDTH_DENOM_LOWER", 4),
    ")\n",
    sep = ""
  )
  flush(stdout())

  results <- vector("list", length(fisheries))
  for (i in seq_along(fisheries)) {
    f <- fisheries[i]
    cat("[", i, "/", length(fisheries), "] LF fishery ", f, "\n", sep = "")
    flush(stdout())
    results[[i]] <- fit_one_lf(
      base_data = inputs$data,
      wcpo_pars = inputs$parameters,
      lf_wide = lf_wide,
      target_fishery = f,
      eval_max = eval_max,
      iter_max = iter_max,
      max_restarts = max_restarts,
      trace = TRUE
    )
  }

  outputs <- write_outputs(results, out_dir)
  print(outputs$metrics |>
    select(
      target_lf_fishery, sel_type, active_sel_cols, n_par, n_lf_obs,
      convergence, max_gradient, nll, B0, lp_cpue, lp_lf,
      log_cpue_rmse, early_log_cpue_rmse,
      lf_prop_rmse, lf_mean_length_rmse, lf_mean_length_bias,
      early_spawning_biomass_ratio, sel_start, sel_est,
      sel_start_nat, sel_est_nat
    ))
}

is_direct_script <- function(script) {
  file_args <- commandArgs(trailingOnly = FALSE)
  file_args <- file_args[startsWith(file_args, "--file=")]
  if (!length(file_args)) return(FALSE)
  normalizePath(sub("^--file=", "", file_args[[1]]), mustWork = FALSE) ==
    normalizePath(script, mustWork = FALSE)
}

if (is_direct_script("inst/scripts/bet_lf_selectivity_diagnostics.R")) {
  main()
}
