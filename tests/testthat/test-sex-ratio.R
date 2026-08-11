library(RTMB)

# Use smaller benchmark dataset for faster tests.
data("opaka_data", envir = environment())
data_base <- opaka_data
data_base$n_index <- 2L

# set ages to anchor growth curve
data_base$A1 <- 0
data_base$A2 <- data_base$A1 + (data_base$n_age-1)
data_base$min_age <- 0
data_base$age_a <- seq(data_base$min_age, by = 1, length.out = data_base$n_age)
data_base$sex_ratio <- rep(1, data_base$n_age)

# LF data
data("opaka_lf", envir = environment())
lf_wide <- opaka_lf %>%
  tidyr::pivot_wider(id_cols = c(fishery, year, month, ts), names_from = bin, values_from = value, values_fill = 0) %>%
  arrange(fishery, ts)
data_base <- prep_lf_data(data_base, lf_wide, lf_keep_fisheries = c(1, 3), lf_var_adjust = rep(1, data_base$n_fishery))
data_base$lf_switch <- 0L

# Bias ramps
data_base$init_bias_adj_a <- get_bias_adj_vector(1949 - (data_base$age_a - 1), 1, c(1931, 1974, 2021.88, 2022.43), 0.869)
data_base$bias_adj_y <- get_bias_adj_vector((1949-1) + data_base$years, 1, c(1931, 1974, 2021.88, 2022.43), 0.869)

# Params base
data("opaka_parameters", envir = environment())
init_rdev_a <- rep(0.0, data_base$n_age)
init_rdev_a[1] <- opaka_parameters$rdev_y[1]
rdev_y_shifted <- c(opaka_parameters$rdev_y[2:length(opaka_parameters$rdev_y)], 0.0)

params_base <- list(
  log_B0 = 9,
  log_h = as.numeric(opaka_parameters$log_h),
  log_sigma_r = as.numeric(opaka_parameters$log_sigma_r),
  log_cpue_q = as.numeric(opaka_parameters$log_cpue_q),
  cpue_creep = as.numeric(rep(opaka_parameters$cpue_creep, data_base$n_index)),
  log_cpue_tau = rep(log(0.1), data_base$n_index),
  log_cpue_omega = as.numeric(rep(opaka_parameters$log_cpue_omega, data_base$n_index)),
  log_lf_tau = as.numeric(log(rep(0.1, data_base$n_fishery))),
  log_L1 = as.numeric(opaka_parameters$log_L1),
  log_L2 = as.numeric(opaka_parameters$log_L2),
  log_k = as.numeric(opaka_parameters$log_k),
  log_CV1 = as.numeric(opaka_parameters$log_CV1),
  log_CV2 = as.numeric(opaka_parameters$log_CV2),
  par_sel = convert_ss3_selex_to_rtmb(as.matrix(opaka_parameters$par_sel), data_base$sel_type_f, data_base$len_mid),
  log_init_F_f = c(log(0.0119122), rep(log(1e-8), data_base$n_fishery - 1)),
  rdev_y = rdev_y_shifted,
  init_rdev_a = init_rdev_a
)

# Construct map from opakapaka vignette example.
map_sel <- matrix(NA, nrow(params_base$par_sel), ncol(params_base$par_sel))
map_sel[1,1:2] <- c(1,2)
map_sel[3,c(1,4,6)] <- c(3,4,5)
map_sel <- as.vector(map_sel)
map_lf_tau <- rep(NA, length(params_base$log_lf_tau))
map_base <- list(
  log_h = factor(NA),
  log_sigma_r = factor(NA),
  cpue_creep = as.factor(rep(NA, data_base$n_index)),
  log_cpue_tau = as.factor(rep(NA, data_base$n_index)),
  log_cpue_omega = as.factor(rep(NA, data_base$n_index)),
  log_lf_tau = factor(map_lf_tau),
  log_L1 = factor(NA),
  log_L2 = factor(NA),
  log_k = factor(NA),
  log_CV1 = factor(NA),
  log_CV2 = factor(NA),
  par_sel = factor(map_sel),
  log_init_F_f = rep(factor(NA), data_base$n_fishery),
  init_rdev_a = as.factor(seq_len(data_base$n_age))
)

obj_base <- MakeADFun(func = cmb(opal_model, data_base), parameters = params_base, map = map_base, silent = TRUE)
rep_base <- obj_base$report()
SPR0_base <- rep_base$B0 / rep_base$R0


# Test 1: age-based sex-ratio scaling.
test_that("sex_ratio directly scales SPR0 and R0 appropriately", {
  data_mod <- data_base
  data_mod$sex_ratio <- rep(0.5, data_base$n_age)

  obj_mod <- MakeADFun(func = cmb(opal_model, data_mod), parameters = params_base, map = map_base, silent = TRUE)
  rep_mod <- obj_mod$report()
  SPR0_mod <- rep_mod$B0 / rep_mod$R0

  # When B0 is fixed, halving sex_ratio should exactly halve SPR0 and double R0
  expect_equal(rep_mod$B0, rep_base$B0)
  expect_equal(SPR0_mod, SPR0_base * 0.5)
  expect_equal(rep_mod$R0, rep_base$R0 * 2)
})


test_that("sex_ratio length-based inputs resolve to age-based equivalent via PLA", {
  # Uniform ratio, but supplied as length-based vector (length = n_len)
  data_mod <- data_base
  data_mod$sex_ratio <- rep(0.5, data_base$n_len)

  obj_mod <- MakeADFun(func = cmb(opal_model, data_mod), parameters = params_base, map = map_base, silent = TRUE)
  rep_mod <- obj_mod$report()

  SPR0_mod <- rep_mod$B0 / rep_mod$R0

  expect_equal(SPR0_mod, SPR0_base * 0.5)
  expect_equal(rep_mod$R0, rep_base$R0 * 2)
})
