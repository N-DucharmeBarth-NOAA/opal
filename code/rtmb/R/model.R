utils::globalVariables(c(
  "par_log_B0", "par_log_h", "par_log_sigma_r", 
  "par_log_cpue_q", "par_cpue_creep", "par_log_cpue_sigma", "par_log_cpue_omega", 
  "par_rdev_y", 
  "n_age", "min_age", "max_age", 
  "first_yr", "last_yr", "first_yr_catch", "n_year", "n_season", "n_fishery",
  "length_m50", "length_m95", "length_mu_ysa", "length_sd_a",
  "removal_switch_f", "weight_fya", "alk_ysal", "dl_yal", "catch_obs_ysf", "af_sliced_ysfa",
  "cpue_switch", "cpue_years", "cpue_lfs", "cpue_n", "cpue_a1", "cpue_a2", "cpue_obs", "cpue_sd",
  "lf_switch", "lf_year", "lf_season", "lf_fishery", "lf_minbin", "lf_obs", "lf_n",
  "priors", "M_a"
))

#' The globals
#' 
#' @return a \code{list} of functions to be passed to \code{sample_sparse_tmb} when doing MCMC.
#' @export
#' 
bet_globals <- function() {
  list(
    posfun = posfun, 
    get_M = get_M, 
    get_rho = get_rho, 
    get_selectivity = get_selectivity, 
    get_selectivity_prior = get_selectivity_prior,
    get_initial_numbers = get_initial_numbers, 
    get_recruitment = get_recruitment, 
    get_harvest_rate = get_harvest_rate, 
    get_length_like = get_length_like, 
    get_cpue_length_like = get_cpue_length_like, 
    get_cpue_like = get_cpue_like, 
    get_recruitment_prior = get_recruitment_prior, 
    evaluate_priors = evaluate_priors)
}

#' The BET model
#' 
#' Obtain the negative log-likelihood (NLL) value from the sbt model.
#' 
#' @param parameters a \code{list} of parameter values.
#' @param data a \code{list} of data inputs.
#' @return the negative log-likelihood (NLL) value.
#' @importFrom RTMB ADoverload getAll REPORT ADREPORT
#' @export
#' 
bet_model <- function(parameters, data) {
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  getAll(data, parameters, warn = FALSE)
  
  # Natural mortality ----
  # M_a is now read from data (extracted from MFCL)
  # M_a <- rep(exp(log_M), n_age)  # OLD: constant M
  # M_a is already available from data via getAll()

  # Selectivity ----
  
  # NOT CODED
  #sel_fya <- get_selectivity(n_age, max_age, first_yr, first_yr_catch, sel_min_age_f, sel_max_age_f, sel_end_f, sel_change_year_fy, par_log_sel_fya)
  sel_fya <- array(1, dim = c(n_fishery, n_year, n_age))
  for (f in seq_len(n_fishery)) {
    for (y in seq_len(n_year)) {
      sel_fya[f, y,] <- maturity_a
    }
  }
  
  # Recruitment ----
  
  # tau_ac2 <- get_rho(first_yr, last_yr, rdev_y)
  tau_ac2 <- 0.6
  lp_rec <- get_recruitment_prior(rdev_y, log_sigma_r, tau_ac2)
  
  # Main population loop ----
  
  B0 <- exp(log_B0)
  h <- exp(log_h)
  init <- get_initial_numbers(B0 = B0, h = h, M_a = M_a, maturity_a = maturity_a)
  R0 <- init$R0
  alpha <- init$alpha
  beta <- init$beta
  
  dyn <- do_dynamics(data, parameters,
                     B0 = B0, R0 = R0, alpha = alpha, beta = beta, h = h,
                     M_a = M_a, init_number_a = init$Ninit, sel_fya = sel_fya)
  
  number_ysa <- dyn$number_ysa
  spawning_biomass_y <- dyn$spawning_biomass_y
  lp_penalty <- dyn$lp_penalty
  
  # plot(dyn$spawning_biomass_y)
  # plot(rowSums(dyn$number_ysa[,1,]))
  # plot(catch_obs_ysf - catch_pred_ysf[1:n_year,,])
  # points(catch_pred_ysf, pch = 2, col = 2)

  # Likelihoods and priors ----

  x <- get_cpue_like(data, parameters, number_ysa, sel_fya)
  lp_cpue <- x$lp
  # cpue_pred <- x$pred

  # x <- get_length_like(lf_switch, removal_switch_f, lf_year, lf_season, lf_fishery, lf_minbin, lf_obs, lf_n, par_log_lf_alpha, catch_pred_fya, alk_ysal)
  # lp_lf <- x$lp
  # lf_pred <- x$pred
  # 
  # x <- get_cpue_length_like(lf_switch, cpue_years, cpue_lfs, cpue_n, par_log_lf_alpha, number_ysa, sel_fya, alk_ysal)
  # lp_cpue_lf <- x$lp
  # cpue_lf_pred <- x$pred
  # 
  lp_prior <- evaluate_priors(parameters, priors)

  # nll <- lp_prior + lp_penalty + lp_rec + sum(lp_cpue_lf) + sum(lp_cpue) + sum(lp_lf)
  nll <- lp_prior + lp_penalty + lp_rec + sum(lp_cpue)
  
  # Reporting ----
  
  REPORT(B0)
  REPORT(R0)
  REPORT(alpha)
  REPORT(beta)
  # ADREPORT(sigma_r)
  # REPORT(tau_ac2)
  # REPORT(par_rdev_y)
  # REPORT(rec_dev_y)
  # REPORT(rdev_y)
  # REPORT(cpue_lf_pred)
  # REPORT(lf_pred)
  REPORT(spawning_biomass_y)
  REPORT(number_ysa)
  REPORT(sel_fya)
  # REPORT(lp_prior)
  # REPORT(lp_penalty)
  # REPORT(lp_rec)
  # REPORT(lp_cpue)
  # REPORT(lp_cpue_lf)
  # REPORT(lp_lf)

  return(nll)
}
