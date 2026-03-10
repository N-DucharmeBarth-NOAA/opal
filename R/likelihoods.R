#' CPUE index likelihood (multi-index)
#'
#' Computes the likelihood for one or more standardised CPUE/survey indices
#' using a log-linear model. Each index has its own catchability (q),
#' extra variance (tau), power parameter (omega), and effort creep.
#' Mean-centering of predicted CPUE is performed within each index.
#'
#' @param data a \code{list} of data inputs. Must contain:
#'   \describe{
#'     \item{cpue_data}{data.frame with columns \code{ts}, \code{fishery},
#'       \code{value}, \code{se}, \code{units}, and \code{index}.}
#'     \item{cpue_switch}{integer switch (0 = skip likelihood).}
#'     \item{n_index}{integer number of distinct indices.}
#'   }
#' @param parameters a \code{list} of parameter values. Must contain:
#'   \describe{
#'     \item{log_cpue_q}{numeric vector \code{[n_index]}.}
#'     \item{log_cpue_tau}{numeric vector \code{[n_index]}.}
#'     \item{log_cpue_omega}{numeric vector \code{[n_index]}.}
#'     \item{cpue_creep}{numeric vector \code{[n_index]}.}
#'   }
#' @param number_ysa a 3D \code{array} `[n_year, n_season, n_age]` of numbers-at-age.
#' @param sel_fya a 3D \code{array} `[n_fishery, n_year, n_age]` of selectivity by fishery, year, and age.
#' @param weight_fya a 3D \code{array} `[n_fishery, n_year, n_age]` of weight-at-age by fishery and year.
#' @param creep_init scalar initialization value for creeping adjustment (default 1).
#' @return numeric vector of length \code{nrow(cpue_data)} with per-observation
#'   negative log-likelihood contributions.
#' @importFrom RTMB ADoverload dnorm
#' @export
#' 
get_cpue_like <- function(data, parameters, number_ysa, sel_fya, weight_fya, creep_init = 1) {
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  getAll(data, parameters, warn = FALSE)

  n_cpue <- nrow(cpue_data)
  cpue_log_pred <- lp <- numeric(n_cpue)
  cpue_sigma <- numeric(n_cpue)
  if (!("index" %in% names(cpue_data))) cpue_data$index <- rep(1L, n_cpue)
  if (!exists("n_index", inherits = FALSE)) n_index <- max(cpue_data$index)

  for (i in seq_len(n_cpue)) {
    y <- cpue_data$ts[i]
    f <- cpue_data$fishery[i]
    cpue_n <- number_ysa[y, 1, ] * sel_fya[f, y, ]
    if (cpue_data$units[i] == 1) cpue_n <- cpue_n * weight_fya[f, y,] # 1=weight, 2=numbers
    sum_n <- sum(cpue_n) + 1e-6
    idx <- cpue_data$index[i]
    cpue_log_pred[i] <- exp(log_cpue_omega[idx]) * log(sum_n)
  }

  for (idx in seq_len(n_index)) {
    rows <- which(cpue_data$index == idx)
    n_idx <- length(rows)
    if (n_idx == 0) next

    cpue_adjust_idx <- numeric(n_idx)
    cpue_adjust_idx[1] <- creep_init
    if (n_idx > 1) {
      for (j in 2:n_idx) {
        cpue_adjust_idx[j] <- cpue_adjust_idx[j - 1] + cpue_creep[idx]
      }
    }
    for (j in seq_len(n_idx)) {
      cpue_log_pred[rows[j]] <- cpue_log_pred[rows[j]] + log(cpue_adjust_idx[j])
    }

    pred_idx <- cpue_log_pred[rows]
    center <- log(mean(exp(pred_idx)))
    for (j in seq_len(n_idx)) {
      cpue_log_pred[rows[j]] <- pred_idx[j] - center + log_cpue_q[idx]
    }

    tau_idx <- exp(log_cpue_tau[idx])
    for (j in seq_len(n_idx)) {
      cpue_sigma[rows[j]] <- sqrt(cpue_data$se[rows[j]]^2 + tau_idx^2)
    }
  }

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
#' @param lf_obs_flat numeric vector of unrounded length composition counts (used when lf_switch=1).
#' @param lf_obs_ints integer vector of integer length composition counts (used when lf_switch=3).
#' @param lf_obs_prop numeric vector of length composition proportions (used when lf_switch=2).
#' @param catch_pred_fya 3D \code{array} `[n_fishery, n_year, n_age]` of predicted
#'   catch-at-age from \code{do_dynamics()}.
#' @param pla matrix `[n_len, n_age]` probability-of-length-at-age from
#'   \code{get_pla()}. On the AD tape when growth parameters are estimated.
#' @param lf_n_f integer vector `[n_fishery]` of number of observations per fishery.
#' @param lf_fishery_f integer vector of fishery indices for each observation group.
#' @param lf_year_fi list of integer vectors of year indices for each observation.
#' @param lf_n_fi list of integer vectors of sample sizes for each observation.
#' @param lf_minbin integer vector `[n_fishery]` of minimum bin index per fishery.
#' @param lf_maxbin integer vector `[n_fishery]` of maximum bin index per fishery.
#' @param removal_switch_f integer vector `[n_fishery]` indicating if fishery is removed (0=included, 1=removed).
#' @param lf_switch integer likelihood type selector (1=multinomial, 2=dirichlet, 3=dirichlet-multinomial).
#' @param n_len integer number of length bins.
#' @param n_lf integer total number of length composition observations.
#' @param log_lf_tau numeric vector `[n_fishery]` of log-scale variance adjustment parameters.
#' @return a \code{numeric} vector of negative log-likelihood contributions, one per observation.
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
    f    <- lf_fishery_f[j]
    ys   <- lf_year_fi[[j]]            # all observed years for fishery f
  
    # ONE matmul: [n_obs_f × n_age] %*% [n_age × n_len] → [n_obs_f × n_len]
    catch_ya <- catch_pred_fya[f, ys, ]    # [n_obs_f × n_age]
    pred_yl  <- catch_ya %*% t(pla)        # [n_obs_f × n_len]  ← replaces the inner loop matvecs

    bmin <- lf_minbin[f]
    bmax <- lf_maxbin[f]
    nbins <- bmax - bmin + 1L
    lf_pred[[j]] <- matrix(0, lf_n_f[j], nbins)
    for (i in seq_len(lf_n_f[j])) {
      idx <- idx + 1L
      y <- lf_year_fi[[j]][i]
      pred <- pred_yl[i, ]                 # simple row slice, no matvec
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

#' Weight Composition Likelihood
#'
#' Computes likelihood for observed weight compositions using the PLA
#' matrix and a precomputed rebinning matrix. The prediction pipeline is:
#' catch_at_age -> PLA -> pred_at_length -> rebin_matrix -> pred_at_weight.
#' Gradients propagate through PLA to growth parameters when estimated.
#'
#' @param wf_obs_flat numeric vector of unrounded weight comp counts (wf_switch=1).
#' @param wf_obs_ints integer vector of integer weight comp counts (wf_switch=3).
#' @param wf_obs_prop numeric vector of weight comp proportions (wf_switch=2).
#' @param catch_pred_fya 3D array `[n_fishery, n_year, n_age]` of predicted
#'   catch-at-age from do_dynamics().
#' @param pla matrix `[n_len, n_age]` probability-of-length-at-age from
#'   get_pla(). On the AD tape when growth parameters are estimated.
#' @param wf_rebin_matrix matrix `[n_wt, n_len]` precomputed rebinning weights
#'   from prep_wf_data().
#' @param wf_n_f integer vector `[n_fishery]` of observation counts per fishery.
#' @param wf_fishery_f integer vector of fishery indices with WF data.
#' @param wf_year_fi list of integer vectors of year indices per fishery.
#' @param wf_n_fi list of integer vectors of sample sizes per fishery.
#' @param wf_minbin integer vector `[n_fishery]` minimum weight bin index.
#' @param wf_maxbin integer vector `[n_fishery]` maximum weight bin index.
#' @param removal_switch_f integer vector `[n_fishery]` removal flags.
#' @param wf_switch integer likelihood type (1=multinomial, 2=Dirichlet, 3=DM).
#' @param n_wt integer number of weight bins.
#' @param n_wf integer total number of WF observations.
#' @param log_wf_tau numeric vector `[n_fishery]` log-scale variance adjustment.
#' @return numeric vector of negative log-likelihood contributions, one per observation.
#' @importFrom RTMB ADoverload dmultinom OBS REPORT
#' @importFrom RTMBdist ddirichlet ddirmult
#' @export
#' 
get_weight_like <- function(wf_obs_flat, wf_obs_ints, wf_obs_prop,
                            catch_pred_fya, pla, wf_rebin_matrix,
                            wf_n_f, wf_fishery_f, wf_year_fi, wf_n_fi,
                            wf_minbin, wf_maxbin, removal_switch_f,
                            wf_switch, n_wt, n_wf, log_wf_tau) {
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")

  n_f <- length(wf_n_f)
  wf_pred <- vector("list", n_f)
  lp <- numeric(n_wf)

  # OBS-mark only the vector used by the active likelihood
  if (wf_switch == 1) wf_obs_flat <- OBS(wf_obs_flat) # unrounded counts
  if (wf_switch == 2) wf_obs_prop <- OBS(wf_obs_prop) # proportions
  if (wf_switch == 3) wf_obs_ints <- OBS(wf_obs_ints) # integer counts

  idx <- 0L
  obs_offset <- 0L
  for (j in seq_len(n_f)) {
    f <- wf_fishery_f[j]
    bmin <- wf_minbin[f]
    bmax <- wf_maxbin[f]
    nbins <- bmax - bmin + 1L
    wf_pred[[j]] <- matrix(0, wf_n_f[j], nbins)

    # Precompute combined age->weight projection once per fishery
    # rebin_pla: [n_wt x n_age] = wf_rebin_matrix [n_wt x n_len] %*% pla [n_len x n_age]
    rebin_pla <- wf_rebin_matrix %*% pla

    # ONE matmul per fishery: [n_obs_f x n_age] %*% [n_age x n_wt] -> [n_obs_f x n_wt]
    ys        <- wf_year_fi[[j]]
    catch_ya  <- catch_pred_fya[f, ys, ]        # [n_obs_f x n_age]
    pred_yw   <- catch_ya %*% t(rebin_pla)      # [n_obs_f x n_wt]

    for (i in seq_len(wf_n_f[j])) {
      idx <- idx + 1L
      pred_at_weight <- pred_yw[i, ]            # simple row slice, no matvec
      if (bmin > 1) pred_at_weight[bmin] <- sum(pred_at_weight[1:bmin])
      if (bmax < n_wt) pred_at_weight[bmax] <- sum(pred_at_weight[bmax:n_wt])
      pred <- pred_at_weight[bmin:bmax]
      pred <- pred + 1e-8
      pred <- pred / sum(pred)
      wf_pred[[j]][i, ] <- pred
      n_i <- wf_n_fi[[j]][i]
      if (removal_switch_f[f] == 0 & n_i > 0) {
        if (wf_switch == 1) { # Multinomial
          obs_i <- wf_obs_flat[(obs_offset + 1):(obs_offset + nbins)]
          lp[idx] <- -RTMB::dmultinom(x = obs_i, prob = pred, log = TRUE)
        }
        if (wf_switch == 2) { # Dirichlet
          obs_i <- wf_obs_prop[(obs_offset + 1):(obs_offset + nbins)]
          alpha_i <- pred * n_i * exp(log_wf_tau[f])
          lp[idx] <- -RTMBdist::ddirichlet(x = obs_i, alpha = alpha_i, log = TRUE)
        }
        if (wf_switch == 3) { # Dirichlet-multinomial
          obs_i <- wf_obs_ints[(obs_offset + 1):(obs_offset + nbins)]
          alpha_i <- pred * exp(log_wf_tau[f])
          lp[idx] <- -RTMBdist::ddirmult(x = obs_i, size = n_i, alpha = alpha_i, log = TRUE)
        }
      }
      obs_offset <- obs_offset + nbins
    }
  }

  REPORT(wf_pred)
  return(lp)
}
