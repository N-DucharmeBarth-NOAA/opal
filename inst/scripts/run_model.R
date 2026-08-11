library(tidyverse)
library(opal)
library(RTMB)

theme_set(theme_bw())

data(wcpo_bet_data)
data <- wcpo_bet_data

ages <- data$age_a
real_age <- ages / 4

data(wcpo_bet_lf)

# Pivot to wide format: one row per fishery x timestep, bins as columns
# Note: LF data exist for fisheries 8-14, but we initially use only fisheries 8 & 9
lf_wide <- wcpo_bet_lf %>%
  tidyr::pivot_wider(
    id_cols = c(fishery, year, month, ts),
    names_from = bin,
    values_from = value,
    values_fill = 0
  ) %>%
  arrange(fishery, ts)

# define variance adjustment scalars for each fishery if using lf_switch = 1 (multinomial)
var_adjust_scalars <- rep(1, data$n_fishery)
# var_adjust_scalars <- 1 / rep(20000, data$n_fishery)
# var_adjust_scalars[c(1, 4, 5, 6, 15)] <- 1 / 40000

# data <- prep_lf_data(data, lf_wide, lf_keep_fisheries = c(8, 9),
#                      lf_var_adjust = var_adjust_scalars)
data <- prep_lf_data(data, lf_wide, lf_keep_fisheries = NULL,
                     lf_var_adjust = var_adjust_scalars)

data(wcpo_bet_wf)

# Weight bin scalars (1 kg bins, 1–200 kg)
data$wt_bin_start <- 1
data$wt_bin_width <- 1
data$n_wt         <- 200L

# Pivot to wide format: one row per fishery x timestep, bins as columns
wf_wide <- wcpo_bet_wf %>%
  tidyr::pivot_wider(
    id_cols     = c(fishery, year, month, ts),
    names_from  = bin,
    values_from = value,
    values_fill = 0
  ) %>%
  dplyr::arrange(fishery, ts)

data <- prep_wf_data(data, wf_wide, wf_keep_fisheries = c(2), wf_switch = 1L)
  
data(wcpo_bet_parameters)

parameters <- list(
  log_B0 = 20,
  log_h = as.numeric(wcpo_bet_parameters$log_h),
  log_sigma_r = as.numeric(wcpo_bet_parameters$log_sigma_r),
  log_cpue_q = as.numeric(wcpo_bet_parameters$log_cpue_q),
  cpue_creep = as.numeric(wcpo_bet_parameters$cpue_creep),
  log_cpue_tau = as.numeric(wcpo_bet_parameters$log_cpue_tau),
  log_cpue_omega = as.numeric(wcpo_bet_parameters$log_cpue_omega),
  log_lf_tau = as.numeric(log(rep(0.1, data$n_fishery))),
  log_wf_tau = rep(0, data$n_fishery),
  log_L1 = as.numeric(wcpo_bet_parameters$log_L1),
  log_L2 = as.numeric(wcpo_bet_parameters$log_L2),
  log_k = as.numeric(wcpo_bet_parameters$log_k),
  log_CV1 = as.numeric(wcpo_bet_parameters$log_CV1),
  log_CV2 = as.numeric(wcpo_bet_parameters$log_CV2),
  par_sel = as.matrix(wcpo_bet_parameters$par_sel),
  rdev_y = as.numeric(wcpo_bet_parameters$rdev_y)
)

data$priors <- get_priors(parameters = parameters, data = data)
evaluate_priors(parameters = parameters, priors = data$priors)

map_sel <- matrix(NA, nrow(parameters$par_sel), ncol(parameters$par_sel))
map_sel[8:14, 1] <- 1:7 # just comment these lines out if you dont want to estimate sels
map_sel[8:14, 3] <- 8:14
map_sel[8:14, 4] <- 15:21
map_rdev <- rep(NA, length(parameters$rdev_y))

map_lf_tau <- rep(NA, length(parameters$log_lf_tau))
# map_lf_tau[c(8, 9)] <- NA

map <- list(
  # log_B0 = factor(NA),
  log_h = factor(NA),
  log_sigma_r = factor(NA),
  # log_cpue_q = factor(NA),
  cpue_creep = factor(NA),
  log_cpue_tau = factor(NA),
  log_cpue_omega = factor(NA),
  log_lf_tau = factor(map_lf_tau),
  log_wf_tau = factor(rep(NA, data$n_fishery)),  # fixed initially
  log_L1  = factor(NA),
  log_L2  = factor(NA),
  log_k   = factor(NA),
  log_CV1 = factor(NA),
  log_CV2 = factor(NA),
  par_sel = factor(map_sel)
  # rdev_y = as.factor(map_rdev),
)

# data$lf_switch <- 0 # skip length comps (removal only)
data$lf_switch <- 1 # multinomial likelihood on flat counts (default)
# data$lf_switch <- 2 # fails while fitting, need to sort out log_lf_tau for this, use simulate to tune and find
# data$lf_switch <- 3 # fails - wants integers - think this should be an issue to RTMBdist guys

# Note: wf_switch was already set by prep_wf_data(); it is re-stated here for
# clarity alongside the equivalent lf_switch assignment.
data$wf_switch <- 1 # skip weight comps (removal only)
# data$wf_switch <- 1 # multinomial likelihood on flat counts (default)
# data$wf_switch <- 2 # Dirichlet
# data$wf_switch <- 3 # Dirichlet-multinomial

obj <- MakeADFun(func = cmb(opal_model, data), parameters = parameters, map = map)
unique(names(obj$par))
obj$fn()
obj$gr()

Lwr <- rep(-Inf, length(obj$par))
Upr <- rep(Inf, length(obj$par))
Lwr[grep("log_B0", names(obj$par))] <- log(1)
Upr[grep("log_B0", names(obj$par))] <- log(exp(22))
Lwr[grep("log_cpue_q", names(obj$par))] <- log(0.1)
Upr[grep("log_cpue_q", names(obj$par))] <- log(10)
Lwr[grep("log_lf_tau", names(obj$par))] <- rep(-9, length(grep("log_lf_tau", names(obj$par))))
Upr[grep("log_lf_tau", names(obj$par))] <- rep(9, length(grep("log_lf_tau", names(obj$par))))
Lwr[grep("rdev_y", names(obj$par))] <- rep(-5, length(grep("rdev_y", names(obj$par))))
Upr[grep("rdev_y", names(obj$par))] <- rep(5, length(grep("rdev_y", names(obj$par))))
# Lwr[grep("par_sel", names(obj$par))] <- rep(-7, length(grep("par_sel", names(obj$par))))
# Upr[grep("par_sel", names(obj$par))] <- rep(7, length(grep("par_sel", names(obj$par))))
bounds <- data.frame(par = names(obj$par), lower = Lwr, upper = Upr)

control <- list(eval.max = 10000, iter.max = 10000)

opt <- nlminb(start = obj$par, objective = obj$fn, gradient = obj$gr,
              hessian = obj$he, lower = Lwr, upper = Upr, control = control)
opt <- nlminb(start = opt$par, objective = obj$fn, gradient = obj$gr,
              hessian = obj$he, lower = Lwr, upper = Upr, control = control)
