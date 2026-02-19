#' CPUE index likelihood
#'
#' Computes the likelihood for a standardized CPUE index using a log-linear model.
#'
#' @param data Integer switch to activate the likelihood.
#' @param parameters a \code{vector} of year indices for CPUE observations.
#' @param number_ysa a 3D \code{array} year, season, age of numbers-at-age.
#' @param sel_fya a 3D \code{array} of selectivity by fishery, year, and age.
#' @return a \code{list} with predicted CPUE, residuals, and likelihood vector.
#' @importFrom RTMB ADoverload dnorm
#' @export
#' 
get_cpue_like <- function(data, parameters, number_ysa, sel_fya, creep_init = 1) {
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  getAll(data, parameters, warn = FALSE)
  cpue_tau <- exp(log_cpue_tau)
  cpue_omega <- exp(log_cpue_omega)
  n_cpue <- nrow(cpue_data)
  cpue_adjust <- cpue_log_pred <- lp <- numeric(n_cpue)
  cpue_adjust[1] <- creep_init
  for (i in 2:n_cpue) cpue_adjust[i] <- cpue_adjust[i - 1] + cpue_creep
  cpue_sigma <- sqrt(cpue_data$se^2 + cpue_tau^2)
  for (i in seq_len(n_cpue)) {
    y <- cpue_data$ts[i]
    # s <- cpue_data$season[i]
    f <- cpue_data$fishery[i]
    cpue_n <- number_ysa[y, 1,] * sel_fya[f, y,]
    if (cpue_data$units[i] == 1) cpue_n <- cpue_n * weight_fya[f, y,] # 1=weight, 2=numbers
    sum_n <- sum(cpue_n) + 1e-6
    cpue_log_pred[i] <- log(cpue_adjust[i]) + cpue_omega * log(sum_n)
  }
  cpue_log_pred <- cpue_log_pred - log(mean(exp(cpue_log_pred))) + log_cpue_q
  cpue_log_obs <- log(cpue_data$value)
  cpue_log_obs <- OBS(cpue_log_obs)
  if (cpue_switch > 0) {
    lp[] <- -dnorm(x = cpue_log_obs, mean = cpue_log_pred, sd = cpue_sigma, log = TRUE)
  }
  cpue_pred <- exp(cpue_log_pred)
  REPORT(cpue_pred)
  REPORT(cpue_sigma)
  return(lp = lp)
}

#' Length Composition Likelihood
#'
#' Computes likelihood for observed length compositions using ALKs.
#'
#' @param lf_switch Vectors of indices for observations.
#' @param removal_switch_f Vectors of indices for observations.
#' @param lf_year,lf_season,lf_fishery Vectors of indices for observations.
#' @param lf_minbin Minimum size bin to be aggregated.
#' @param lf_obs Matrix of observed length proportions.
#' @param lf_n Vector of effective sample sizes.
#' @param par_log_lf_alpha Vector of effective sample sizes.
#' @param catch_pred_fya 3D array of predicted catch.
#' @param alk_ysal 4D array year, season, age, length_bin of ALKs.
#' @return List with predicted compositions and likelihood contributions.
#' @importFrom RTMB dmultinom
#' @importFrom RTMBdist ddirichlet ddirmult
#' @export
#' 
get_length_like <- function(data, parameters, log_lf_alpha, catch_pred_fya, alk_ysal) {
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  getAll(data, parameters, warn = FALSE)
  n_lf <- nrow(lf_obs)
  n_bins <- 25
  n_age <- dim(catch_pred_fya)[3]
  lp <- numeric(n_lf)
  lf_pred <- matrix(0, n_lf, n_bins)
  for (i in seq_len(n_lf)) {
    y <- lf_year[i]
    s <- lf_season[i]
    f <- lf_fishery[i]
    catch_a <- catch_pred_fya[f, y,]
    pred <- catch_a %*% alk_ysal[y, s,, 1:n_bins]
    pred <- pred[1,] / sum(pred)
    obs <- lf_obs[i,]
    mbin <- lf_minbin[f]
    if (mbin > 1) {
      pred[mbin] <- sum(pred[1:mbin])
      obs[mbin] <- sum(obs[1:mbin])
      obs <- obs[mbin:n_bins]
      pred <- pred[mbin:n_bins]
    }
    lf_pred[i, mbin:n_bins] <- pred
    if (removal_switch_f[f] == 0 & lf_n[i] > 0) {
      if (lf_switch == 1) {
        # obs <- obs + 1e-6
        obs <- obs * lf_n[i]
        pred <- pred + 1e-6
        pred <- pred / sum(pred)
        lp[i] <- -dmultinom(x = obs, prob = pred, log = TRUE)
      }
      if (lf_switch == 2) {
        obs <- obs + 1e-6
        obs <- obs / sum(obs)
        pred <- pred + 1e-6
        pred <- pred / sum(pred) * lf_n[i] * exp(par_log_lf_alpha[f])
        lp[i] <- -ddirichlet(x = obs, alpha = pred, log = TRUE)
      }
      if (lf_switch == 3) {
        obs <- obs * lf_n[i]
        pred <- pred + 1e-6
        pred <- pred / sum(pred) * exp(par_log_lf_alpha[f])
        lp[i] <- -ddirmult(x = obs, size = lf_n[i], alpha = pred, log = TRUE)
      }
      if (lf_switch == 9) {
        pred <- pred + 1e-6 # I added this is as I was having Inf issues with pred
        lp[i] <- -lf_n[i] * sum(obs * log(pred))
        obs <- obs + 1e-6
        lp[i] <- lp[i] + lf_n[i] * sum(obs * log(obs))
      }
    }
  }
  REPORT(lf_pred)
  return(lp)
}
