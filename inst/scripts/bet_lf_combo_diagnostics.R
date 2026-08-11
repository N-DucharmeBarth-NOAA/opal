#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(RTMB)
  library(dplyr)
  library(tidyr)
  library(purrr)
})

helper <- new.env(parent = globalenv())
sys.source("inst/scripts/bet_lf_selectivity_diagnostics.R", helper)

prepare_combo_lf_data <- function(data, lf_wide, fisheries,
                                  lf_var_adjust_scalar = 80) {
  lf_var_adjust <- rep(lf_var_adjust_scalar, data$n_fishery)
  data <- prep_lf_data(
    data = data,
    lf_wide = lf_wide,
    lf_keep_fisheries = fisheries,
    lf_var_adjust = lf_var_adjust,
    lf_switch = 1L
  )
  data$lf_switch <- 1L
  data$wf_switch <- 0L
  data$n_wf <- 0L
  data$priors <- list()
  data
}

make_combo_map <- function(data, parameters, fisheries) {
  map_sel <- matrix(NA_integer_, nrow(parameters$par_sel), ncol(parameters$par_sel))
  level <- 1L
  for (f in fisheries) {
    cols <- helper$selectivity_columns(data, f)
    map_sel[f, cols] <- seq.int(level, length.out = length(cols))
    level <- level + length(cols)
  }

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

  list(map = map, map_sel = map_sel)
}

build_combo_fit <- function(base_data, wcpo_pars, lf_wide, fisheries) {
  data <- base_data |>
    helper$add_length_bins() |>
    helper$set_quarter_cpue() |>
    prepare_combo_lf_data(lf_wide = lf_wide, fisheries = fisheries)

  parameters <- helper$make_parameters(data, wcpo_pars)
  map_info <- make_combo_map(data, parameters, fisheries)
  obj <- MakeADFun(
    func = cmb(opal_model, data),
    parameters = parameters,
    map = map_info$map
  )
  obj$env$tracemgc <- FALSE
  bounds <- helper$make_bounds(obj, data, map_info$map, map_info$map_sel)

  list(
    data = data,
    parameters = parameters,
    map = map_info$map,
    map_sel = map_info$map_sel,
    obj = obj,
    bounds = bounds
  )
}

max_kkt_gradient <- function(par, grad, lower, upper, tol = 1e-7) {
  out <- abs(grad)
  at_lower <- is.finite(lower) & par <= lower + tol
  at_upper <- is.finite(upper) & par >= upper - tol
  out[at_lower] <- pmax(0, -grad[at_lower])
  out[at_upper] <- pmax(0, grad[at_upper])
  max(out)
}

fit_combo <- function(base_data, wcpo_pars, lf_wide, fisheries,
                      eval_max = 1000L, iter_max = 1000L,
                      max_restarts = 4L, use_hessian = FALSE) {
  built <- build_combo_fit(base_data, wcpo_pars, lf_wide, fisheries)
  obj <- built$obj
  opt <- list(par = obj$par, convergence = NA_integer_, message = "not run")
  control <- list(eval.max = eval_max, iter.max = iter_max)

  helper$print_progress(
    "combo ", paste(fisheries, collapse = ","),
    ": start nll=", signif(tryCatch(obj$fn(obj$par), error = function(e) NaN), 8),
    "\n"
  )

  for (i in seq_len(max_restarts)) {
    args <- list(
      start = opt$par,
      objective = obj$fn,
      gradient = obj$gr,
      lower = built$bounds$lower,
      upper = built$bounds$upper,
      control = control
    )
    if (isTRUE(use_hessian)) args$hessian <- obj$he
    opt <- do.call(nlminb, args)
    grad <- obj$gr(opt$par)
    raw_grad <- max(abs(grad))
    kkt_grad <- max_kkt_gradient(opt$par, grad, built$bounds$lower, built$bounds$upper)
    helper$print_progress(
      "  restart ", i,
      ": convergence=", opt$convergence,
      ", nll=", signif(opt$objective, 8),
      ", raw_grad=", signif(raw_grad, 6),
      ", kkt_grad=", signif(kkt_grad, 6),
      "\n"
    )
    if (is.finite(kkt_grad) && kkt_grad < 1e-2) break
  }

  obj$env$last.par.best <- opt$par
  rep <- obj$report(opt$par)
  par_list <- obj$env$parList(opt$par)
  grad <- obj$gr(opt$par)

  list(
    combo = paste(fisheries, collapse = "_"),
    fisheries = fisheries,
    fit = built,
    opt = opt,
    report = rep,
    par_list = par_list,
    raw_grad = max(abs(grad)),
    kkt_grad = max_kkt_gradient(opt$par, grad, built$bounds$lower, built$bounds$upper)
  )
}

make_combo_lf_df <- function(data, rep, combo) {
  map_dfr(seq_along(data$lf_fishery_f), function(j) {
    f <- data$lf_fishery_f[j]
    lf_pred <- rep$lf_pred[[j]]
    obs_counts <- data$lf_obs_data$obs[[j]]
    obs_prop <- obs_counts / rowSums(obs_counts)
    bmin <- data$lf_minbin[f]
    bmax <- data$lf_maxbin[f]
    len <- data$len_mid[bmin:bmax]
    years <- data$lf_year_fi[[as.character(f)]]

    pred_df <- as.data.frame(lf_pred)
    names(pred_df) <- len
    obs_df <- as.data.frame(obs_prop)
    names(obs_df) <- len

    bind_rows(
      obs_df |>
        mutate(obs_id = row_number(), year = years, type = "observed"),
      pred_df |>
        mutate(obs_id = row_number(), year = years, type = "predicted")
    ) |>
      pivot_longer(
        cols = -c(obs_id, year, type),
        names_to = "length",
        values_to = "proportion"
      ) |>
      mutate(
        combo = combo,
        fishery = f,
        length = as.numeric(length)
      )
  })
}

selectivity_bound_status <- function(result, fishery) {
  built <- result$fit
  obj <- built$obj
  map_sel <- built$map_sel
  cols <- which(!is.na(map_sel[fishery, ]))
  if (!length(cols)) {
    return(tibble(
      param = character(),
      estimate = numeric(),
      lower = numeric(),
      upper = numeric(),
      at_lower = logical(),
      at_upper = logical()
    ))
  }
  par_sel_pos <- which(names(obj$par) == "par_sel")
  par_sel_level_pos <- setNames(par_sel_pos, levels(built$map$par_sel))
  labels <- c("peak", "plateau", "asc_width", "desc_width", "init", "final")

  map_dfr(cols, function(p) {
    level <- as.character(map_sel[fishery, p])
    idx <- par_sel_level_pos[[level]]
    est <- unname(result$opt$par[idx])
    lo <- built$bounds$lower[idx]
    up <- built$bounds$upper[idx]
    tibble(
      param = labels[p],
      estimate = est,
      lower = lo,
      upper = up,
      at_lower = is.finite(lo) && est <= lo + 1e-6,
      at_upper = is.finite(up) && est >= up - 1e-6
    )
  })
}

make_metrics <- function(result) {
  data <- result$fit$data
  rep <- result$report
  par_list <- result$par_list
  cpue_obs <- data$cpue_data$value
  cpue_pred <- as.numeric(rep$cpue_pred)
  cpue_log_resid <- log(cpue_obs) - log(cpue_pred)
  early_n <- min(8L, length(cpue_obs))
  sb <- as.numeric(rep$spawning_biomass_y)

  tibble(
    combo = result$combo,
    fisheries = paste(result$fisheries, collapse = ","),
    n_lf_fisheries = length(result$fisheries),
    n_par = length(result$fit$obj$par),
    n_lf_obs = data$n_lf,
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
    log_cpue_rmse = sqrt(mean(cpue_log_resid^2)),
    early_log_cpue_rmse = sqrt(mean(cpue_log_resid[seq_len(early_n)]^2)),
    early_spawning_biomass_ratio = sb[min(9L, length(sb))] / sb[1],
    final_spawning_biomass_ratio = tail(sb, 1) / as.numeric(rep$B0),
    q_est = paste(round(as.numeric(par_list$log_cpue_q), 4), collapse = ",")
  )
}

make_lf_metrics <- function(result) {
  lf_df <- make_combo_lf_df(result$fit$data, result$report, result$combo)
  lf_wide <- lf_df |>
    select(combo, fishery, obs_id, year, length, type, proportion) |>
    pivot_wider(names_from = type, values_from = proportion)

  mean_len <- lf_wide |>
    group_by(combo, fishery, obs_id, year) |>
    summarise(
      obs_mean_length = sum(length * observed),
      pred_mean_length = sum(length * predicted),
      .groups = "drop"
    )

  lf_stats <- lf_wide |>
    group_by(combo, fishery) |>
    summarise(
      lf_prop_rmse = sqrt(mean((observed - predicted)^2)),
      lf_prop_mae = mean(abs(observed - predicted)),
      .groups = "drop"
    )

  mean_stats <- mean_len |>
    group_by(combo, fishery) |>
    summarise(
      lf_mean_length_rmse = sqrt(mean((obs_mean_length - pred_mean_length)^2)),
      lf_mean_length_bias = mean(obs_mean_length - pred_mean_length),
      .groups = "drop"
    )

  metrics <- left_join(lf_stats, mean_stats, by = c("combo", "fishery"))
  metrics |>
    mutate(
      sel_est_nat = map_chr(fishery, function(f) {
        helper$selectivity_natural(
          result$fit$data,
          result$par_list$par_sel[f, ],
          f
        )
      }),
      bound_status = map_chr(fishery, function(f) {
        status_df <- selectivity_bound_status(result, f)
        if (!nrow(status_df)) return("fixed")
        status_df |>
          mutate(status = case_when(
            at_lower ~ paste0(param, "=lower"),
            at_upper ~ paste0(param, "=upper"),
            TRUE ~ paste0(param, "=free")
          )) |>
          pull(status) |>
          paste(collapse = ";")
      })
    )
}

benchmark_hessian <- function(base_data, wcpo_pars, lf_wide) {
  built <- build_combo_fit(base_data, wcpo_pars, lf_wide, c(9L, 10L))
  no_hessian <- system.time(
    opt_no_hessian <- nlminb(
      start = built$obj$par,
      objective = built$obj$fn,
      gradient = built$obj$gr,
      lower = built$bounds$lower,
      upper = built$bounds$upper,
      control = list(eval.max = 5, iter.max = 5)
    )
  )
  with_hessian <- system.time(
    opt_with_hessian <- nlminb(
      start = built$obj$par,
      objective = built$obj$fn,
      gradient = built$obj$gr,
      hessian = built$obj$he,
      lower = built$bounds$lower,
      upper = built$bounds$upper,
      control = list(eval.max = 5, iter.max = 5)
    )
  )

  tibble(
    use_hessian = c(FALSE, TRUE),
    elapsed_sec = c(no_hessian[["elapsed"]], with_hessian[["elapsed"]]),
    objective = c(opt_no_hessian$objective, opt_with_hessian$objective),
    convergence = c(opt_no_hessian$convergence, opt_with_hessian$convergence)
  )
}

parse_combo_list <- function(value) {
  if (!nzchar(value)) return(NULL)
  strsplit(value, ";", fixed = TRUE)[[1]] |>
    map(function(group) {
      as.integer(strsplit(group, ",", fixed = TRUE)[[1]])
    })
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out_dir <- if (length(args) >= 1L) args[[1]] else "inst/scripts/bet_lf_combo_diagnostics"
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  helper$source_opal_r()
  inputs <- helper$load_bet_inputs()
  lf_wide <- helper$make_lf_wide(inputs$lf)
  eval_max <- as.integer(helper$num_env("OPAL_BET_COMBO_EVAL_MAX", 1000))
  iter_max <- as.integer(helper$num_env("OPAL_BET_COMBO_ITER_MAX", 1000))
  max_restarts <- as.integer(helper$num_env("OPAL_BET_COMBO_MAX_RESTARTS", 4))

  combos <- parse_combo_list(Sys.getenv("OPAL_BET_COMBO_LIST", unset = ""))
  if (is.null(combos)) {
    combos <- list(
      c(9L, 10L),
      c(8L, 9L, 10L),
      c(9L, 10L, 11L),
      c(9L, 10L, 12L),
      c(9L, 10L, 13L),
      c(9L, 10L, 14L),
      c(9L, 10L, 11L, 12L, 13L, 14L)
    )
  }

  if (identical(Sys.getenv("OPAL_BET_COMBO_BENCH_HESSIAN"), "true")) {
    hessian_benchmark <- benchmark_hessian(inputs$data, inputs$parameters, lf_wide)
    write.csv(hessian_benchmark, file.path(out_dir, "hessian_benchmark.csv"), row.names = FALSE)
    print(hessian_benchmark)
  }
  cat(
    "Running ", length(combos), " LF combos",
    " (eval.max=", eval_max,
    ", iter.max=", iter_max,
    ", max_restarts=", max_restarts,
    ")\n",
    sep = ""
  )
  flush(stdout())

  results <- map(combos, function(fisheries) {
    fit_combo(
      base_data = inputs$data,
      wcpo_pars = inputs$parameters,
      lf_wide = lf_wide,
      fisheries = fisheries,
      eval_max = eval_max,
      iter_max = iter_max,
      max_restarts = max_restarts,
      use_hessian = FALSE
    )
  })

  combo_metrics <- bind_rows(map(results, make_metrics))
  lf_metrics <- bind_rows(map(results, make_lf_metrics))
  saveRDS(results, file.path(out_dir, "combo_results.rds"))
  write.csv(combo_metrics, file.path(out_dir, "combo_metrics.csv"), row.names = FALSE)
  write.csv(lf_metrics, file.path(out_dir, "combo_lf_metrics.csv"), row.names = FALSE)

  print(combo_metrics |>
    select(
      combo, convergence, raw_max_gradient, kkt_max_gradient, B0,
      lp_cpue, lp_lf, early_log_cpue_rmse, final_spawning_biomass_ratio
    ))
  print(lf_metrics |>
    select(combo, fishery, lf_prop_rmse, lf_mean_length_rmse,
           lf_mean_length_bias, bound_status, sel_est_nat))
}

is_direct_script <- function(script) {
  file_args <- commandArgs(trailingOnly = FALSE)
  file_args <- file_args[startsWith(file_args, "--file=")]
  if (!length(file_args)) return(FALSE)
  normalizePath(sub("^--file=", "", file_args[[1]]), mustWork = FALSE) ==
    normalizePath(script, mustWork = FALSE)
}

if (is_direct_script("inst/scripts/bet_lf_combo_diagnostics.R")) {
  main()
}
