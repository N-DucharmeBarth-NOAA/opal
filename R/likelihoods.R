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
get_length_like <- function(data, parameters, catch_pred_fya, pla) {
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  getAll(data, parameters, warn = FALSE)
  
  n_f <- length(lf_n_f)
  
  lf_year_fi <- split(lf_year, lf_fishery)
  lf_n_fi <- split(lf_n, lf_fishery)
  
  n_lf <- nrow(lf_obs_in)
  n_len_local <- ncol(lf_obs_in)
  # lp <- numeric(n_lf)
  # lf_pred <- matrix(0, n_lf, n_len_local)
  lp <- lf_obs <- lf_pred <- vector("list", n_f)
  
  for (j in seq_len(n_f)) {
    f <- lf_fishery_f[j]
    # Handle bin aggregation
    bmin <- lf_minbin[f]
    bmax <- lf_maxbin[f]
    # lf_obs[[f]] <- matrix(0, lf_n_f[j], bmax - bmin + 1)
    lf_obs[[j]] <- lf_pred[[j]] <- matrix(0, lf_n_f[j], bmax - bmin + 1)
    # f <- lf_fishery[i]
    for (i in seq_len(lf_n_f[j])) {
      y <- lf_year_fi[[j]][i]
      # Predicted catch-at-age for this fishery and timestep
      catch_a <- catch_pred_fya[f, y, ]
      # Convert to predicted numbers-at-length via PLA
      # pla is [n_len, n_age], catch_a is [n_age] -> result is [n_len]
      pred <- c(pla %*% catch_a)
      # Observed proportions for this row
      obs <- lf_obs_in[i, ]
      # Lower-tail aggregation: collapse bins 1:bmin into bin bmin
      if (bmin > 1) {
        pred[bmin] <- sum(pred[1:bmin])
        obs[bmin] <- sum(obs[1:bmin])
      }
      # Upper-tail aggregation: collapse bins bmax:n_len into bin bmax
      if (bmax < n_len_local) {
        pred[bmax] <- sum(pred[bmax:n_len_local])
        obs[bmax] <- sum(obs[bmax:n_len_local])
      }
      # Subset to active bins
      pred <- pred[bmin:bmax]
      obs <- obs[bmin:bmax]
      # Normalize predicted to proportions (epsilon for AD safety)
      pred <- pred + 1e-8
      pred <- pred / sum(pred)
      lf_obs[[j]][i, ] <- obs * lf_n_fi[[j]][i]# * exp(log_lf_tau[f])
      lf_pred[[j]][i, ] <- pred
      # if (lf_switch == 1) lf_obs[i, bmin:bmax] <- obs * lf_n[i] * exp(log_lf_tau[f])
      # lf_pred[i, bmin:bmax] <- pred # Store predicted proportions for diagnostics
    }
  }
  
  # lf_obs <- OBS(lf_obs)
  
  # Evaluate likelihood
  
  # for (i in seq_len(n_lf)) {
  for (j in seq_len(n_f)) {
    f <- lf_fishery_f[j]
    bmin <- lf_minbin[f]
    bmax <- lf_maxbin[f]
    lp[[j]] <- numeric(lf_n_f[j])
    for (i in seq_len(lf_n_f[j])) {
      if (removal_switch_f[f] == 0 & lf_n_fi[[j]][i] > 0) {
        if (lf_switch == 1) { # Multinomial
          # obs <- obs * lf_n[i] * exp(log_lf_tau[f])
          # lp[i] <- -RTMB::dmultinom(x = obs, prob = pred, log = TRUE) # REMOVE LATER, if you source this files then reverts to stats::dmultinom which return zero as log like
          lp[[j]][i] <- -RTMB::dmultinom(x = lf_obs[[j]][i, ], prob = lf_pred[[j]][i, ], log = TRUE) # REMOVE LATER, if you source this files then reverts to stats::dmultinom which return zero as log like
        }
        # if (lf_switch == 2) { # Dirichlet
        #   obs <- obs + 1e-8
        #   obs <- obs / sum(obs)
        #   pred <- pred + 1e-8
        #   pred <- pred / sum(pred)
        #   lf_alpha <- pred * lf_n[i] * exp(log_lf_tau[f])
        #   lp[i] <- -ddirichlet(x = obs, alpha = lf_alpha, log = TRUE)
        # }
        # if (lf_switch == 3) { # Dirichlet-multinomial
        #   obs <- obs * lf_n[i]
        #   lf_alpha <- pred * exp(log_lf_tau[f])
        #   lp[i] <- -ddirmult(x = obs, size = lf_n[i], alpha = lf_alpha, log = TRUE)
        # }
        # if (lf_switch == 9) { # Multinomial offset form
        #   #   -n_eff * sum(obs * log(pred)) + n_eff * sum(obs * log(obs))
        #   # Equivalent to n_eff * KL(obs || pred); minimum NLL = 0 when pred == obs.
        #   # Small constant added to both terms for numerical safety and exact
        #   # cancellation at perfect fit.
        #   n_eff <- lf_n[i] * exp(log_lf_tau[f])
        #   lp[i] <- -n_eff * sum(obs * log(pred + 1e-8))
        #   lp[i] <- lp[i] + n_eff * sum(obs * log(obs + 1e-8))
      }
    }
  }
  
  REPORT(lf_pred)
  return(unlist(lp))
}
