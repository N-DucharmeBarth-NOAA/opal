opaka_quickstart_inputs <- function() {
  data(opaka_data, package = "opal", envir = environment())
  data(opaka_lf, package = "opal", envir = environment())
  data(opaka_parameters, package = "opal", envir = environment())

  model_data <- opaka_data
  model_data$n_index <- 2L
  model_data$priors <- NULL
  model_data$A1 <- 0
  model_data$A2 <- model_data$A1 + model_data$n_age - 1L
  model_data$min_age <- 0
  model_data$age_a <- seq(model_data$min_age, by = 1, length.out = model_data$n_age)
  model_data$sex_ratio <- rep(1, model_data$n_age)

  lf_wide <- tidyr::pivot_wider(
    opaka_lf,
    id_cols = c(fishery, year, month, ts),
    names_from = bin,
    values_from = value,
    values_fill = 0
  )
  lf_wide <- dplyr::arrange(lf_wide, fishery, ts)
  model_data <- prep_lf_data(
    model_data,
    lf_wide,
    lf_keep_fisheries = c(1, 3),
    lf_var_adjust = rep(1, model_data$n_fishery),
    lf_addtocomp = 1e-10
  )
  model_data$lf_switch <- 1L

  bias_years <- c(1931, 1974, 2021.88, 2022.43)
  model_data$init_bias_adj_a <- get_bias_adj_vector(
    years = 1949 - (model_data$age_a - 1),
    do_rec_bias_ramp = 1,
    bias_years = bias_years,
    max_bias_adj = 0.869
  )
  model_data$bias_adj_y <- get_bias_adj_vector(
    years = (1949 - 1) + model_data$years,
    do_rec_bias_ramp = 1,
    bias_years = bias_years,
    max_bias_adj = 0.869
  )

  parameters <- list(
    log_B0 = 9,
    log_h = as.numeric(opaka_parameters$log_h),
    log_sigma_r = as.numeric(opaka_parameters$log_sigma_r),
    log_cpue_q = as.numeric(opaka_parameters$log_cpue_q),
    cpue_creep = as.numeric(rep(opaka_parameters$cpue_creep, model_data$n_index)),
    log_cpue_tau = rep(log(0.1), model_data$n_index),
    log_cpue_omega = as.numeric(rep(opaka_parameters$log_cpue_omega, model_data$n_index)),
    log_lf_tau = rep(log(0.1), model_data$n_fishery),
    log_L1 = as.numeric(opaka_parameters$log_L1),
    log_L2 = as.numeric(opaka_parameters$log_L2),
    log_k = as.numeric(opaka_parameters$log_k),
    log_CV1 = as.numeric(opaka_parameters$log_CV1),
    log_CV2 = as.numeric(opaka_parameters$log_CV2),
    par_sel = convert_ss3_selex_to_rtmb(
      as.matrix(opaka_parameters$par_sel), model_data$sel_type_f, model_data$len_mid
    ),
    log_init_F_f = c(log(0.0119122), rep(log(1e-8), model_data$n_fishery - 1L)),
    rdev_y = c(opaka_parameters$rdev_y[-1], 0),
    init_rdev_a = c(opaka_parameters$rdev_y[1], rep(0, model_data$n_age - 1L))
  )

  selectivity_map <- matrix(NA_integer_, nrow(parameters$par_sel), ncol(parameters$par_sel))
  selectivity_map[1, 1:2] <- 1:2
  selectivity_map[3, c(1, 4, 6)] <- 3:5
  map <- list(
    log_h = factor(NA),
    log_sigma_r = factor(NA),
    cpue_creep = factor(rep(NA, model_data$n_index)),
    log_cpue_tau = factor(rep(NA, model_data$n_index)),
    log_cpue_omega = factor(rep(NA, model_data$n_index)),
    log_lf_tau = factor(rep(NA, length(parameters$log_lf_tau))),
    log_L1 = factor(NA), log_L2 = factor(NA), log_k = factor(NA),
    log_CV1 = factor(NA), log_CV2 = factor(NA),
    par_sel = factor(as.vector(selectivity_map)),
    log_init_F_f = factor(rep(NA, model_data$n_fishery)),
    init_rdev_a = factor(seq_len(model_data$n_age)),
    rdev_y = factor(seq_len(length(parameters$rdev_y)))
  )
  model_data$priors <- get_priors(parameters, model_data)

  list(data = model_data, parameters = parameters, map = map)
}