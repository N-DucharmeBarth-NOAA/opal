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
#' @param parameters a \code{list} of parameter values.
#' @param catch_pred_fya 3D \code{array} \code{[n_fishery, n_year, n_age]} of predicted
#'   catch-at-age from \code{do_dynamics()}.
#' @param pla Matrix \code{[n_len, n_age]} probability-of-length-at-age from
#'   \code{get_pla()}. On the AD tape when growth parameters are estimated.
#' @return a \code{vector} of negative log-likelihood contributions, one per observation.
#' @importFrom RTMB ADoverload REPORT dmultinom
#' @export
#' 
get_length_like <- function(data, parameters, catch_pred_fya, pla) {
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  getAll(data, parameters, warn = FALSE)

  n_lf <- nrow(lf_obs)
  n_len_local <- ncol(lf_obs)
  lp <- numeric(n_lf)
  lf_pred <- matrix(0, n_lf, n_len_local)

  for (i in seq_len(n_lf)) {
    f <- lf_fishery[i]
    y <- lf_year[i]

    # Predicted catch-at-age for this fishery and timestep
    catch_a <- catch_pred_fya[f, y, ]

    # Convert to predicted numbers-at-length via PLA
    # pla is [n_len, n_age], catch_a is [n_age] -> result is [n_len]
    pred <- c(pla %*% catch_a)

    # Observed proportions for this row
    obs <- lf_obs[i, ]

    # Handle minimum bin aggregation
    mbin <- lf_minbin[f]
    if (mbin > 1) {
      pred[mbin] <- sum(pred[1:mbin])
      obs[mbin] <- sum(obs[1:mbin])
      pred <- pred[mbin:n_len_local]
      obs <- obs[mbin:n_len_local]
    }

    # Normalize predicted to proportions (epsilon for AD safety)
    pred <- pred + 1e-8
    pred <- pred / sum(pred)

    # Store predicted proportions for diagnostics
    lf_pred[i, mbin:n_len_local] <- pred

    # Evaluate likelihood
    if (removal_switch_f[f] == 0 & lf_n[i] > 0) {

      # Effective sample size = raw sample size * variance adjustment factor
      n_eff <- lf_n[i] * lf_var_adj[f]

      if (lf_switch == 1) {
        # Multinomial via dmultinom
        obs_counts <- obs * n_eff
        lp[i] <- -dmultinom(x = obs_counts, prob = pred, log = TRUE)
      }

      if (lf_switch == 9) {
        # Multinomial offset form:
        #   -n_eff * sum(obs * log(pred)) + n_eff * sum(obs * log(obs))
        # Equivalent to n_eff * KL(obs || pred); minimum NLL = 0 when pred == obs.
        # Small constant added to both terms for numerical safety and exact
        # cancellation at perfect fit.
        lp[i] <- -n_eff * sum(obs * log(pred + 1e-8))
        lp[i] <- lp[i] + n_eff * sum(obs * log(obs + 1e-8))
      }

      # --- Alternative likelihoods (inactive, retained for future use) ---
      # if (lf_switch == 2) {
      #   # Dirichlet
      #   obs <- obs + 1e-8
      #   obs <- obs / sum(obs)
      #   alpha <- pred * n_eff * exp(par_log_lf_alpha[f])
      #   lp[i] <- -ddirichlet(x = obs, alpha = alpha, log = TRUE)
      # }
      # if (lf_switch == 3) {
      #   # Dirichlet-multinomial
      #   obs_counts <- obs * n_eff
      #   alpha <- pred * exp(par_log_lf_alpha[f])
      #   lp[i] <- -ddirmult(x = obs_counts, size = n_eff, alpha = alpha, log = TRUE)
      # }
    }
  }

  REPORT(lf_pred)
  return(lp)
}
