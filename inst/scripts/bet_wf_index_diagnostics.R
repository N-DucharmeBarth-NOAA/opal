#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(RTMB)
  library(dplyr)
  library(tidyr)
  library(purrr)
})

combo <- new.env(parent = globalenv())
sys.source("inst/scripts/bet_lf_combo_diagnostics.R", combo)

make_wf_wide <- function(wf) {
  wf |>
    pivot_wider(
      id_cols = c(fishery, year, month, ts),
      names_from = bin,
      values_from = value,
      values_fill = 0
    ) |>
    arrange(fishery, ts)
}

parse_int_values <- function(value, default) {
  if (!nzchar(value)) return(default)
  as.integer(strsplit(value, ",", fixed = TRUE)[[1]])
}

default_sel_specs <- function(index_sel_cols = c(1L, 2L)) {
  list(
    `9` = c(1L, 3L, 4L),
    `10` = c(1L, 3L, 4L),
    `14` = 1L,
    `15` = as.integer(index_sel_cols)
  )
}

parse_sel_specs <- function(value, default) {
  if (!nzchar(value)) return(default)
  groups <- strsplit(value, ";", fixed = TRUE)[[1]]
  out <- list()
  for (group in groups) {
    parts <- strsplit(group, ":", fixed = TRUE)[[1]]
    stopifnot("Selectivity specs must use fishery:col,col syntax" = length(parts) == 2L)
    out[[parts[[1]]]] <- as.integer(strsplit(parts[[2]], ",", fixed = TRUE)[[1]])
  }
  out
}

parse_named_numeric <- function(value) {
  if (!nzchar(value)) return(numeric())
  groups <- strsplit(value, ";", fixed = TRUE)[[1]]
  out <- numeric()
  for (group in groups) {
    parts <- strsplit(group, ":", fixed = TRUE)[[1]]
    stopifnot("Named numeric specs must use name:value syntax" = length(parts) == 2L)
    out[[parts[[1]]]] <- as.numeric(parts[[2]])
  }
  out
}

format_sel_specs <- function(sel_specs) {
  paste(
    paste0(names(sel_specs), ":", map_chr(sel_specs, ~ paste(.x, collapse = ","))),
    collapse = ";"
  )
}

prepare_lf_wf_data <- function(data, lf_wide, wf_wide,
                               lf_fisheries = c(9L, 10L, 14L),
                               wf_fisheries = 15L,
                               lf_var_adjust_scalar = 80,
                               lf_var_adjust_override = numeric(),
                               wf_var_adjust_scalar = 2000) {
  data <- combo$helper$add_length_bins(data)
  data <- combo$helper$set_quarter_cpue(data)

  lf_var_adjust <- rep(lf_var_adjust_scalar, data$n_fishery)
  if (length(lf_var_adjust_override)) {
    override_f <- as.integer(names(lf_var_adjust_override))
    lf_var_adjust[override_f] <- as.numeric(lf_var_adjust_override)
  }
  data <- prep_lf_data(
    data = data,
    lf_wide = lf_wide,
    lf_keep_fisheries = lf_fisheries,
    lf_var_adjust = lf_var_adjust,
    lf_switch = 1L
  )

  data$wt_bin_start <- 1
  data$wt_bin_width <- 1
  data$n_wt <- 200L
  wf_var_adjust <- rep(wf_var_adjust_scalar, data$n_fishery)
  data <- prep_wf_data(
    data = data,
    wf_wide = wf_wide,
    wf_keep_fisheries = wf_fisheries,
    wf_switch = 1L,
    wf_var_adjust = wf_var_adjust
  )

  data$lf_switch <- 1L
  data$wf_switch <- 1L
  data$priors <- list()
  data
}

make_index_map <- function(data, parameters,
                           estimate_index_sel = TRUE,
                           index_sel_cols = c(1L, 2L),
                           sel_specs = NULL,
                           estimate_lf14_peak = TRUE) {
  map_sel <- matrix(NA_integer_, nrow(parameters$par_sel), ncol(parameters$par_sel))
  level <- 1L

  if (is.null(sel_specs)) {
    sel_specs <- default_sel_specs(index_sel_cols = index_sel_cols)
    if (!isTRUE(estimate_lf14_peak)) sel_specs[["14"]] <- NULL
    if (!isTRUE(estimate_index_sel)) sel_specs[["15"]] <- NULL
  }

  for (f_chr in names(sel_specs)) {
    f <- as.integer(f_chr)
    cols <- as.integer(sel_specs[[f_chr]])
    map_sel[f, cols] <- seq.int(level, length.out = length(cols))
    level <- level + length(cols)
  }

  map_lf_tau <- rep(NA_integer_, data$n_fishery)
  map_wf_tau <- rep(NA_integer_, data$n_fishery)

  map <- list(
    log_h = factor(NA),
    log_sigma_r = factor(NA),
    log_cpue_q = factor(seq_len(data$n_index)),
    cpue_creep = factor(rep(NA, data$n_index)),
    log_cpue_tau = factor(rep(NA, data$n_index)),
    log_cpue_omega = factor(rep(NA, data$n_index)),
    log_lf_tau = factor(map_lf_tau),
    log_wf_tau = factor(map_wf_tau),
    log_L1 = factor(NA),
    log_L2 = factor(NA),
    log_k = factor(NA),
    log_CV1 = factor(NA),
    log_CV2 = factor(NA),
    par_sel = factor(map_sel),
    rdev_y = factor(seq_along(parameters$rdev_y))
  )

  list(map = map, map_sel = map_sel)
}

build_fit <- function(base_data, wcpo_pars, lf_wide, wf_wide,
                      wf_var_adjust_scalar,
                      estimate_index_sel = TRUE,
                      index_sel_cols = c(1L, 2L),
                      lf_fisheries = c(9L, 10L, 14L),
                      wf_fisheries = 15L,
                      lf_var_adjust_override = numeric(),
                      sel_specs = NULL) {
  data <- prepare_lf_wf_data(
    data = base_data,
    lf_wide = lf_wide,
    wf_wide = wf_wide,
    lf_fisheries = lf_fisheries,
    wf_fisheries = wf_fisheries,
    lf_var_adjust_override = lf_var_adjust_override,
    wf_var_adjust_scalar = wf_var_adjust_scalar
  )
  parameters <- combo$helper$make_parameters(data, wcpo_pars)
  map_info <- make_index_map(
    data = data,
    parameters = parameters,
    estimate_index_sel = estimate_index_sel,
    index_sel_cols = index_sel_cols,
    sel_specs = sel_specs
  )
  obj <- MakeADFun(
    func = cmb(opal_model, data),
    parameters = parameters,
    map = map_info$map
  )
  obj$env$tracemgc <- FALSE
  bounds <- combo$helper$make_bounds(obj, data, map_info$map, map_info$map_sel)

  list(
    data = data,
    parameters = parameters,
    map = map_info$map,
    map_sel = map_info$map_sel,
    obj = obj,
    bounds = bounds
  )
}

fit_one <- function(base_data, wcpo_pars, lf_wide, wf_wide,
                    wf_var_adjust_scalar,
                    estimate_index_sel = TRUE,
                    index_sel_cols = c(1L, 2L),
                    lf_fisheries = c(9L, 10L, 14L),
                    wf_fisheries = 15L,
                    lf_var_adjust_override = numeric(),
                    sel_specs = NULL,
                    label = NULL,
                    eval_max = 1000L, iter_max = 1000L,
                    max_restarts = 4L) {
  if (is.null(sel_specs)) sel_specs <- default_sel_specs(index_sel_cols = index_sel_cols)
  if (is.null(label) || !nzchar(label)) {
    label <- paste0("wf", wf_var_adjust_scalar)
  }
  built <- build_fit(
    base_data = base_data,
    wcpo_pars = wcpo_pars,
    lf_wide = lf_wide,
    wf_wide = wf_wide,
    wf_var_adjust_scalar = wf_var_adjust_scalar,
    estimate_index_sel = estimate_index_sel,
    index_sel_cols = index_sel_cols,
    lf_fisheries = lf_fisheries,
    wf_fisheries = wf_fisheries,
    lf_var_adjust_override = lf_var_adjust_override,
    sel_specs = sel_specs
  )
  obj <- built$obj
  opt <- list(par = obj$par, convergence = NA_integer_, message = "not run")
  control <- list(eval.max = eval_max, iter.max = iter_max)

  cat(
    "wf_adjust=", wf_var_adjust_scalar,
    ", lf_fisheries=", paste(lf_fisheries, collapse = ","),
    ", lf_adjust_override=", paste(paste(names(lf_var_adjust_override), lf_var_adjust_override, sep = ":"), collapse = ";"),
    ", wf_fisheries=", paste(wf_fisheries, collapse = ","),
    ", sel_specs=", format_sel_specs(sel_specs),
    ": start nll=", signif(obj$fn(obj$par), 8), "\n",
    sep = ""
  )
  flush(stdout())

  for (i in seq_len(max_restarts)) {
    opt <- nlminb(
      start = opt$par,
      objective = obj$fn,
      gradient = obj$gr,
      lower = built$bounds$lower,
      upper = built$bounds$upper,
      control = control
    )
    grad <- obj$gr(opt$par)
    raw_grad <- max(abs(grad))
    kkt_grad <- combo$max_kkt_gradient(
      opt$par, grad, built$bounds$lower, built$bounds$upper
    )
    cat(
      "  restart ", i,
      ": convergence=", opt$convergence,
      ", nll=", signif(opt$objective, 8),
      ", raw_grad=", signif(raw_grad, 6),
      ", kkt_grad=", signif(kkt_grad, 6),
      "\n",
      sep = ""
    )
    flush(stdout())
    if (is.finite(kkt_grad) && kkt_grad < 1e-2) break
  }

  obj$env$last.par.best <- opt$par
  rep <- obj$report(opt$par)
  par_list <- obj$env$parList(opt$par)
  grad <- obj$gr(opt$par)
  list(
    label = label,
    wf_var_adjust = wf_var_adjust_scalar,
    lf_var_adjust_override = lf_var_adjust_override,
    estimate_index_sel = estimate_index_sel,
    index_sel_cols = index_sel_cols,
    lf_fisheries = lf_fisheries,
    wf_fisheries = wf_fisheries,
    sel_specs = sel_specs,
    fit = built,
    opt = opt,
    report = rep,
    par_list = par_list,
    raw_grad = max(abs(grad)),
    kkt_grad = combo$max_kkt_gradient(
      opt$par, grad, built$bounds$lower, built$bounds$upper
    )
  )
}

make_wf_df <- function(data, rep, label) {
  map_dfr(seq_along(data$wf_fishery_f), function(j) {
    f <- data$wf_fishery_f[j]
    wf_pred <- rep$wf_pred[[j]]
    bmin <- data$wf_minbin[f]
    bmax <- data$wf_maxbin[f]
    wt <- data$wt_mid[bmin:bmax]
    rows <- data$wf_row_fi[[as.character(f)]]
    obs_prop <- as.matrix(data$wf_obs_in[rows, bmin:bmax, drop = FALSE])
    obs_prop <- obs_prop / rowSums(obs_prop)
    years <- data$wf_year_fi[[as.character(f)]]

    pred_df <- as.data.frame(wf_pred)
    names(pred_df) <- wt
    obs_df <- as.data.frame(obs_prop)
    names(obs_df) <- wt

    bind_rows(
      obs_df |>
        mutate(obs_id = row_number(), year = years, type = "observed"),
      pred_df |>
        mutate(obs_id = row_number(), year = years, type = "predicted")
    ) |>
      pivot_longer(
        cols = -c(obs_id, year, type),
        names_to = "weight",
        values_to = "proportion"
      ) |>
      mutate(
        label = label,
        fishery = f,
        weight = as.numeric(weight)
      )
  })
}

selectivity_status <- function(result, fishery) {
  combo$selectivity_bound_status(
    list(
      fit = result$fit,
      opt = result$opt,
      par_list = result$par_list
    ),
    fishery
  )
}

make_metrics <- function(result) {
  data <- result$fit$data
  rep <- result$report
  cpue_obs <- data$cpue_data$value
  cpue_pred <- as.numeric(rep$cpue_pred)
  cpue_log_resid <- log(cpue_obs) - log(cpue_pred)
  early_n <- min(8L, length(cpue_obs))
  sb <- as.numeric(rep$spawning_biomass_y)

  sel_fisheries <- as.integer(names(result$sel_specs))
  sel_est <- map_chr(sel_fisheries, function(f) {
    paste0(
      "F", f, ":",
      combo$helper$selectivity_natural(data, result$par_list$par_sel[f, ], f)
    )
  })
  sel_status <- map_chr(sel_fisheries, function(f) {
    status_df <- selectivity_status(result, f)
    if (!nrow(status_df)) return(paste0("F", f, ":fixed"))
    paste0(
      "F", f, ":",
      status_df |>
        mutate(status = case_when(
          at_lower ~ paste0(param, "=lower"),
          at_upper ~ paste0(param, "=upper"),
          TRUE ~ paste0(param, "=free")
        )) |>
        pull(status) |>
        paste(collapse = ";")
    )
  })

  tibble(
    label = result$label,
    wf_var_adjust = result$wf_var_adjust,
    estimate_index_sel = result$estimate_index_sel,
    index_sel_cols = paste(result$index_sel_cols, collapse = ","),
    lf_fisheries = paste(result$lf_fisheries, collapse = ","),
    wf_fisheries = paste(result$wf_fisheries, collapse = ","),
    sel_specs = format_sel_specs(result$sel_specs),
    n_par = length(result$fit$obj$par),
    n_lf_obs = data$n_lf,
    n_wf_obs = data$n_wf,
    wf_eff_n_total = sum(data$wf_n),
    convergence = result$opt$convergence,
    message = result$opt$message,
    nll = result$fit$obj$fn(result$opt$par),
    raw_max_gradient = result$raw_grad,
    kkt_max_gradient = result$kkt_grad,
    B0 = as.numeric(rep$B0),
    lp_prior = as.numeric(rep$lp_prior),
    lp_penalty = as.numeric(rep$lp_penalty),
    lp_rec = as.numeric(rep$lp_rec),
    lp_cpue = sum(as.numeric(rep$lp_cpue)),
    lp_lf = sum(as.numeric(rep$lp_lf)),
    lp_wf = sum(as.numeric(rep$lp_wf)),
    log_cpue_rmse = sqrt(mean(cpue_log_resid^2)),
    early_log_cpue_rmse = sqrt(mean(cpue_log_resid[seq_len(early_n)]^2)),
    final_spawning_biomass_ratio = tail(sb, 1) / as.numeric(rep$B0),
    sel_est = paste(sel_est, collapse = " | "),
    sel_status = paste(sel_status, collapse = " | ")
  )
}

make_wf_metrics <- function(result) {
  wf_df <- make_wf_df(result$fit$data, result$report, result$label)
  wf_wide <- wf_df |>
    select(label, fishery, obs_id, year, weight, type, proportion) |>
    pivot_wider(names_from = type, values_from = proportion)

  mean_wt <- wf_wide |>
    group_by(label, fishery, obs_id, year) |>
    summarise(
      obs_mean_weight = sum(weight * observed),
      pred_mean_weight = sum(weight * predicted),
      .groups = "drop"
    )

  left_join(
    wf_wide |>
      group_by(label, fishery) |>
      summarise(
        wf_prop_rmse = sqrt(mean((observed - predicted)^2)),
        wf_prop_mae = mean(abs(observed - predicted)),
        .groups = "drop"
      ),
    mean_wt |>
      group_by(label, fishery) |>
      summarise(
        wf_mean_weight_rmse = sqrt(mean((obs_mean_weight - pred_mean_weight)^2)),
        wf_mean_weight_bias = mean(obs_mean_weight - pred_mean_weight),
        .groups = "drop"
      ),
    by = c("label", "fishery")
  )
}

make_lf_metrics <- function(result) {
  result2 <- result
  result2$combo <- result$label
  combo$make_lf_metrics(result2)
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out_dir <- if (length(args) >= 1L) args[[1]] else "inst/scripts/bet_wf_index_diagnostics"
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  combo$helper$source_opal_r()
  inputs <- combo$helper$load_bet_inputs()
  load("data/wcpo_bet_wf.rda")
  wf_wide <- make_wf_wide(wcpo_bet_wf)
  lf_wide <- combo$helper$make_lf_wide(inputs$lf)

  wf_adjust_values <- as.numeric(strsplit(
    Sys.getenv("OPAL_BET_WF_ADJUST_VALUES", unset = "500,2000,5000,10000"),
    ",",
    fixed = TRUE
  )[[1]])
  eval_max <- as.integer(combo$helper$num_env("OPAL_BET_WF_EVAL_MAX", 1000))
  iter_max <- as.integer(combo$helper$num_env("OPAL_BET_WF_ITER_MAX", 1000))
  max_restarts <- as.integer(combo$helper$num_env("OPAL_BET_WF_MAX_RESTARTS", 4))
  index_sel_cols <- parse_int_values(Sys.getenv("OPAL_BET_INDEX_SEL_COLS", unset = "1,2"), c(1L, 2L))
  lf_fisheries <- parse_int_values(Sys.getenv("OPAL_BET_LF_FISHERIES", unset = "9,10,14"), c(9L, 10L, 14L))
  wf_fisheries <- parse_int_values(Sys.getenv("OPAL_BET_WF_FISHERIES", unset = "15"), 15L)
  sel_specs <- parse_sel_specs(
    Sys.getenv("OPAL_BET_SEL_SPECS", unset = ""),
    default_sel_specs(index_sel_cols = index_sel_cols)
  )
  lf_var_adjust_override <- parse_named_numeric(
    Sys.getenv("OPAL_BET_LF_VAR_ADJUST_OVERRIDE", unset = "")
  )
  scenario_label <- Sys.getenv("OPAL_BET_SCENARIO_LABEL", unset = "")

  results <- map(wf_adjust_values, function(w) {
    fit_one(
      base_data = inputs$data,
      wcpo_pars = inputs$parameters,
      lf_wide = lf_wide,
      wf_wide = wf_wide,
      wf_var_adjust_scalar = w,
      index_sel_cols = index_sel_cols,
      lf_fisheries = lf_fisheries,
      wf_fisheries = wf_fisheries,
      lf_var_adjust_override = lf_var_adjust_override,
      sel_specs = sel_specs,
      label = if (nzchar(scenario_label)) paste0(scenario_label, "_wf", w) else NULL,
      eval_max = eval_max,
      iter_max = iter_max,
      max_restarts = max_restarts
    )
  })

  metrics <- bind_rows(map(results, make_metrics))
  lf_metrics <- bind_rows(map(results, make_lf_metrics))
  wf_metrics <- bind_rows(map(results, make_wf_metrics))
  saveRDS(results, file.path(out_dir, "wf_index_results.rds"))
  write.csv(metrics, file.path(out_dir, "wf_index_metrics.csv"), row.names = FALSE)
  write.csv(lf_metrics, file.path(out_dir, "wf_index_lf_metrics.csv"), row.names = FALSE)
  write.csv(wf_metrics, file.path(out_dir, "wf_index_wf_metrics.csv"), row.names = FALSE)

  print(metrics |>
    select(
      label, wf_var_adjust, wf_eff_n_total, convergence, raw_max_gradient,
      kkt_max_gradient, B0, lp_cpue, lp_lf, lp_wf, early_log_cpue_rmse,
      final_spawning_biomass_ratio, sel_est, sel_status
    ))
  print(lf_metrics |>
    select(combo, fishery, lf_prop_rmse, lf_mean_length_rmse,
           lf_mean_length_bias, bound_status, sel_est_nat))
  print(wf_metrics)
}

is_direct_script <- function(script) {
  file_args <- commandArgs(trailingOnly = FALSE)
  file_args <- file_args[startsWith(file_args, "--file=")]
  if (!length(file_args)) return(FALSE)
  normalizePath(sub("^--file=", "", file_args[[1]]), mustWork = FALSE) ==
    normalizePath(script, mustWork = FALSE)
}

if (is_direct_script("inst/scripts/bet_wf_index_diagnostics.R")) {
  main()
}
