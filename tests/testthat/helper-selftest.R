apply_sim_obs <- function(data, sim) {
  stopifnot(all(c("cpue_log_obs", "lf_obs_flat") %in% names(sim)))

  data$cpue_data$value <- exp(as.numeric(sim$cpue_log_obs))
  data$lf_obs_flat <- as.numeric(sim$lf_obs_flat)
  data$lf_obs_ints <- as.integer(round(data$lf_obs_flat))

  offset <- 0L
  proportions <- numeric(length(data$lf_obs_flat))
  for (fishery_index in seq_along(data$lf_fishery_f)) {
    fishery <- data$lf_fishery_f[fishery_index]
    n_bin <- data$lf_maxbin[fishery] - data$lf_minbin[fishery] + 1L
    for (observation_index in seq_len(data$lf_n_f[fishery_index])) {
      indices <- offset + seq_len(n_bin)
      observation <- data$lf_obs_flat[indices]
      proportion <- observation / sum(observation)
      proportions[indices] <- (proportion + data$lf_addtocomp) /
        sum(proportion + data$lf_addtocomp)
      offset <- offset + n_bin
    }
  }
  data$lf_obs_prop <- proportions
  data
}

run_selftest <- function(n_sim, start = c("truth", "default"),
                         seed0 = 20260916L) {
  start <- match.arg(start)
  inputs <- opaka_inputs()
  om <- opaka_obj(inputs)
  om_opt <- fit_quickstart(om, inputs$parameters)
  truth_par <- om$env$last.par.best
  truth <- om$report(truth_par)
  terminal_year <- inputs$data$n_year + 1L

  rows <- lapply(seq_len(n_sim), function(simulation) {
    set.seed(seed0 + simulation)
    data <- apply_sim_obs(inputs$data, om$simulate(par = truth_par))
    em <- opaka_obj(inputs, data = data)
    initial <- if (start == "truth") om_opt$par else em$par
    opt <- tryCatch(
      fit_quickstart(em, inputs$parameters, start = initial),
      error = identity
    )
    if (inherits(opt, "error")) {
      return(data.frame(
        sim = simulation, converged = FALSE, max_gr = NA_real_,
        re_B0 = NA_real_, re_sb_term = NA_real_, re_dep_term = NA_real_,
        max_abs_re_sb = NA_real_
      ))
    }

    report <- em$report(em$env$last.par.best)
    max_gradient <- max(abs(em$gr(opt$par)))
    data.frame(
      sim = simulation,
      converged = opt$convergence == 0 && max_gradient <= 0.01,
      max_gr = max_gradient,
      re_B0 = report$B0 / truth$B0 - 1,
      re_sb_term = report$spawning_biomass_y[terminal_year] /
        truth$spawning_biomass_y[terminal_year] - 1,
      re_dep_term = report$static_depletion_y[terminal_year] /
        truth$static_depletion_y[terminal_year] - 1,
      max_abs_re_sb = max(abs(
        report$spawning_biomass_y / truth$spawning_biomass_y - 1
      ))
    )
  })
  do.call(rbind, rows)
}

summarize_selftest <- function(results, n_boot = 2000L, seed = 20260916L) {
  stopifnot(is.data.frame(results), n_boot > 0)
  required <- c(
    "converged", "re_B0", "re_sb_term", "re_dep_term", "max_abs_re_sb"
  )
  stopifnot(all(required %in% names(results)))

  converged <- results[results$converged, , drop = FALSE]
  metrics <- c("re_B0", "re_sb_term", "re_dep_term", "max_abs_re_sb")
  if (!nrow(converged)) {
    return(data.frame(
      metric = metrics,
      n_sim = nrow(results),
      n_converged = 0L,
      convergence_rate = 0,
      median = NA_real_,
      mc_se = NA_real_
    ))
  }

  set.seed(seed)
  do.call(rbind, lapply(metrics, function(metric) {
    values <- converged[[metric]]
    boot_medians <- replicate(n_boot, median(sample(values, replace = TRUE)))
    data.frame(
      metric = metric,
      n_sim = nrow(results),
      n_converged = nrow(converged),
      convergence_rate = nrow(converged) / nrow(results),
      median = median(values),
      mc_se = stats::sd(boot_medians)
    )
  }))
}