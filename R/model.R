utils::globalVariables(c(
  "log_B0", "log_h", "sigma_r", 
  "log_cpue_q", "cpue_creep", "log_cpue_sigma", "log_cpue_omega", 
  "rdev_y", 
  "par_sel",
  "log_L1", "log_L2", "log_k", "log_CV1", "log_CV2",
  "n_age", "min_age", "max_age", 
  "first_yr", "last_yr", "first_yr_catch", "n_year", "n_season", "n_fishery",
  "M",
  "A1", "A2", "lw_a", "lw_b", "maturity", "len_bin_start", "len_bin_width",
  "length_m50", "length_m95", "length_mu_ysa", "length_sd_a",
  "removal_switch_f", "alk_ysal", "dl_yal", "catch_obs_ysf", "af_sliced_ysfa",
  "cpue_switch", "cpue_years", "cpue_n", "cpue_obs", "cpue_sd",
  "lf_switch", "lf_year", "lf_season", "lf_fishery", "lf_minbin", "lf_obs", "lf_n",
  "priors"
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
    get_growth = get_growth,
    get_sd_at_age = get_sd_at_age,
    get_weight_at_length = get_weight_at_length,
    get_maturity_at_age = get_maturity_at_age,
    resolve_bio_vector = resolve_bio_vector,
    get_selectivity = get_selectivity,
    sel_logistic = sel_logistic,
    sel_double_normal = sel_double_normal,
    get_pla = get_pla,
    get_initial_numbers = get_initial_numbers, 
    get_recruitment = get_recruitment, 
    get_harvest_rate = get_harvest_rate, 
    get_length_like = get_length_like, 
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

  # Derive length bins from scalars (single source of truth) ----

  len_lower <- seq(from = len_bin_start, by = len_bin_width, length.out = n_len)
  len_upper <- len_lower + len_bin_width
  len_mid   <- len_lower + len_bin_width / 2

  # Growth module ----

  # Back-transform growth/variability parameters
  L1  <- exp(log_L1)
  L2  <- exp(log_L2)
  k   <- exp(log_k)
  CV1 <- exp(log_CV1)
  CV2 <- exp(log_CV2)

  # Module 1: Mean length-at-age (Schnute VB)
  mu_a <- get_growth(n_age, A1, A2, L1, L2, k)

  # Module 2: SD of length-at-age (linear CV interpolation)
  sd_a <- get_sd_at_age(mu_a, L1, L2, CV1, CV2)

  # Shared PLA — computed once and reused for weight, maturity, selectivity
  pla <- get_pla(len_lower, len_upper, mu_a, sd_a)

  # Module 3: Weight-at-length → weight-at-age
  wt_at_len <- get_weight_at_length(len_mid, lw_a, lw_b)
  weight_a  <- as.vector(t(pla) %*% wt_at_len)
  # Replicate weight across fisheries and years (AD-safe: use loop + [<- overload)
  weight_fya_mod <- array(0, dim = c(n_fishery, n_year, n_age))
  for (f in seq_len(n_fishery)) {
    for (y in seq_len(n_year)) {
      weight_fya_mod[f, y, ] <- weight_a
    }
  }

  # Module 4: Resolve biology vectors to age-basis via PLA ----
  # Accepts either age-basis (length n_age) or length-basis (length n_len) vectors.
  # Length-basis vectors are converted using: vec_a = t(pla) %*% vec_l
  maturity_a <- resolve_bio_vector(maturity, n_age, n_len, pla, "maturity")
  M_a        <- resolve_bio_vector(M, n_age, n_len, pla, "M")

  # Selectivity ----

  # mu_a and sd_a from growth module so AD gradients propagate if growth is estimated
  sel_fya <- get_selectivity(data, par_sel, pla, len_mid)

  # Main population loop ----

  B0 <- exp(log_B0)
  h <- exp(log_h)
  init <- get_initial_numbers(B0 = B0, h = h, M_a = M_a, maturity_a = maturity_a)
  R0 <- init$R0
  alpha <- init$alpha
  beta <- init$beta
  sigma_r <- exp(log_sigma_r)
  REPORT(B0)
  REPORT(R0)
  REPORT(alpha)
  REPORT(beta)
  REPORT(sigma_r)

  dyn <- do_dynamics(data, parameters,
                     B0 = B0, R0 = R0, alpha = alpha, beta = beta, h = h, sigma_r = sigma_r,
                     M_a = M_a, maturity_a = maturity_a, weight_fya = weight_fya_mod,
                     init_number_a = init$Ninit, sel_fya = sel_fya)

  number_ysa <- dyn$number_ysa
  lp_penalty <- dyn$lp_penalty

  # plot(spawning_biomass_y)
  # plot(rowSums(dyn$number_ysa[,1,]))
  # plot(catch_obs_ysf - catch_pred_ysf)
  # points(catch_pred_ysf, pch = 2, col = 2)

  # Priors ----

  # lp_prior <- evaluate_priors(parameters, priors)
  lp_rec <- get_recruitment_prior(rdev_y, sigma_r)
  lp_prior <- 0

  # Likelihoods ----

  lp_cpue <- get_cpue_like(data, parameters, number_ysa, sel_fya)
  # lp_lf <- get_length_like(lf_switch, removal_switch_f, lf_year, lf_season, lf_fishery, lf_minbin, lf_obs, lf_n, par_log_lf_alpha, catch_pred_fya, alk_ysal)

  nll <- lp_prior + lp_penalty + lp_rec + sum(lp_cpue)# + sum(lp_lf)

  # Reporting ----

  REPORT(number_ysa)
  REPORT(sel_fya)

  REPORT(lp_prior)
  REPORT(lp_penalty)
  REPORT(lp_rec)
  REPORT(lp_cpue)
  # REPORT(lp_lf)

  return(nll)
}
