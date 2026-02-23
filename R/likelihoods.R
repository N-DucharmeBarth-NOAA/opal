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
#' Computes likelihood for observed length compositions using a probability-of-length-at-age
#' (PLA) matrix. Gradients propagate through to growth parameters when they are estimated.
#'
#' @param data a \code{list} of data inputs.
#'   \describe{
#'     \item{lf_maxbin}{Integer vector \[n_fishery\]. Per-fishery maximum bin index.
#'       Bins from lf_maxbin\[f\] to n_len are aggregated into bin lf_maxbin\[f\].
#'       Default n_len (no upper aggregation).}
#'   }
#' @param parameters a \code{list} of parameter values.
#' @param catch_pred_fya 3D \code{array} \code{[n_fishery, n_year, n_age]} of predicted
#'   catch-at-age from \code{do_dynamics()}.
#' @param pla Matrix \code{[n_len, n_age]} probability-of-length-at-age from
#'   \code{get_pla()}. On the AD tape when growth parameters are estimated.
#' @return a \code{vector} of negative log-likelihood contributions, one per observation.
#' @importFrom RTMB ADoverload dmultinom OBS REPORT
#' @importFrom RTMBdist ddirichlet ddirmult
#' @export
#' 
#         # if (lf_switch == 9) { # Multinomial offset form
#         #   #   -n_eff * sum(obs * log(pred)) + n_eff * sum(obs * log(obs))
#         #   # Equivalent to n_eff * KL(obs || pred); minimum NLL = 0 when pred == obs.
#         #   # Small constant added to both terms for numerical safety and exact
#         #   # cancellation at perfect fit.
#         #   n_eff <- lf_n[i] * exp(log_lf_tau[f])
#         #   lp[i] <- -n_eff * sum(obs * log(pred + 1e-8))
#         #   lp[i] <- lp[i] + n_eff * sum(obs * log(obs + 1e-8))
get_length_like <- function(lf_obs_flat, lf_obs_ints, lf_obs_prop,
                            catch_pred_fya, pla,
                            lf_n_f, lf_fishery_f, lf_year_fi, lf_n_fi,
                            lf_minbin, lf_maxbin, removal_switch_f,
                            lf_switch, n_len, n_lf, log_lf_tau) {
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  
  n_f <- length(lf_n_f)
  lf_pred <- vector("list", n_f)
  lp <- numeric(n_lf)
  
  # OBS-mark only the vector used by the active likelihood
  if (lf_switch == 1) lf_obs_flat <- OBS(lf_obs_flat) # unrounded counts
  if (lf_switch == 2) lf_obs_prop <- OBS(lf_obs_prop) # proportions
  if (lf_switch == 3) lf_obs_ints <- OBS(lf_obs_ints) # integer counts
  
  idx <- 0L
  obs_offset <- 0L
  for (j in seq_len(n_f)) {
    f <- lf_fishery_f[j]
    bmin <- lf_minbin[f]
    bmax <- lf_maxbin[f]
    nbins <- bmax - bmin + 1L
    lf_pred[[j]] <- matrix(0, lf_n_f[j], nbins)
    for (i in seq_len(lf_n_f[j])) {
      idx <- idx + 1L
      y <- lf_year_fi[[j]][i]
      catch_a <- catch_pred_fya[f, y, ]
      pred <- c(pla %*% catch_a)
      if (bmin > 1) pred[bmin] <- sum(pred[1:bmin])
      if (bmax < n_len) pred[bmax] <- sum(pred[bmax:n_len])
      pred <- pred[bmin:bmax]
      pred <- pred + 1e-8
      pred <- pred / sum(pred)
      lf_pred[[j]][i, ] <- pred
      n_i <- lf_n_fi[[j]][i]
      if (removal_switch_f[f] == 0 & n_i > 0) {
        if (lf_switch == 1) { # Multinomial
          obs_i <- lf_obs_flat[(obs_offset + 1):(obs_offset + nbins)]
          lp[idx] <- -RTMB::dmultinom(x = obs_i, prob = pred, log = TRUE)
        }
        if (lf_switch == 2) { # Dirichlet
          obs_i <- lf_obs_prop[(obs_offset + 1):(obs_offset + nbins)]
          alpha_i <- pred * n_i * exp(log_lf_tau[f])
          lp[idx] <- -RTMBdist::ddirichlet(x = obs_i, alpha = alpha_i, log = TRUE)
        }
        if (lf_switch == 3) { # Dirichlet-multinomial
          obs_i <- lf_obs_ints[(obs_offset + 1):(obs_offset + nbins)]
          alpha_i <- pred * exp(log_lf_tau[f])
          lp[idx] <- -RTMBdist::ddirmult(x = obs_i, size = n_i, alpha = alpha_i, log = TRUE)
        }
      }
      obs_offset <- obs_offset + nbins
    }
  }
  
  REPORT(lf_pred)
  return(lp)
}
