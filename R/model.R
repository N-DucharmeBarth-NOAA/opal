utils::globalVariables(c(
  "log_B0", "log_h", "log_sigma_r", "sigma_r", 
  "log_cpue_q", "log_cpue_omega",  
  "cpue_creep",
  "rdev_y", 
  "par_sel",
  "log_L1", "log_L2", "log_k", "log_CV1", "log_CV2",
  "n_age", "min_age", "max_age", 
  "first_yr", "last_yr", "first_yr_catch", "n_year", "n_season", "n_fishery",
  "M",
  "A1", "A2", "lw_a", "lw_b", "maturity", "fecundity", "len_bin_start", "len_bin_width",
  "length_m50", "length_m95", "length_mu_ysa", "length_sd_a",
  "removal_switch_f", "alk_ysal", "dl_yal", "catch_obs_ysf", "af_sliced_ysfa",
  "cpue_switch", "cpue_data", "n_index",
  "lf_switch", "lf_year", "lf_season", "lf_fishery", "lf_minbin", "lf_maxbin", "lf_obs", "lf_n",
  "lf_var_adj", "lf_addtocomp",
  "wf_switch", "wf_obs_flat", "wf_obs_ints", "wf_obs_prop",
  "wf_n_f", "wf_fishery_f", "wf_fishery", "wf_year", "wf_n",
  "wf_minbin", "wf_maxbin", "wf_rebin_matrix", "n_wf", "n_wt",
  "wt_bin_start", "wt_bin_width",
  "log_wf_tau", "wf_addtocomp",
  "priors", "sel_fa_external",
  "length_dynamics", "plus_group_growth", "growth_sd_mode", "sd_floor",
  "len_lower", "len_upper", "len_mid", "weight", "spawning_potential",
  "sex_ratio", "fecundity", "init_rdev_a", "init_bias_adj_a",
  "log_init_F_f", "bias_adj_y", "log_lf_tau", "log_wf_tau",
  "lf_year", "lf_fishery", "lf_n", "lf_obs_flat", "lf_obs_ints", "lf_obs_prop",
  "lf_n_f", "lf_fishery_f", "removal_switch_f", "n_lf", "n_len",
  "wf_year_fi", "wf_n_fi", "wf_obs_flat", "wf_obs_ints", "wf_obs_prop",
  "wf_n_f", "wf_fishery_f"
))

#' The opal globals
#' 
#' @return a \code{list} of functions to be passed to \code{sample_sparse_tmb} when doing MCMC.
#' @export
#' 
opal_globals <- function() {
  list(
    posfun = posfun, 
    get_M = get_M, 
    get_rho = get_rho, 
    get_growth = get_growth,
    get_sd_at_age = get_sd_at_age,
    get_weight_at_length = get_weight_at_length,
    get_maturity_at_age = get_maturity_at_age,
    get_recruit_length_dist = get_recruit_length_dist,
    get_growth_matrix = get_growth_matrix,
    resolve_bio_vector = resolve_bio_vector,
    get_selectivity = get_selectivity,
    get_selectivity_length = get_selectivity_length,
    sel_logistic = sel_logistic,
    sel_double_normal = sel_double_normal,
    get_pla = get_pla,
    get_initial_numbers = get_initial_numbers,
    get_initial_numbers_length = get_initial_numbers_length,
    get_recruitment = get_recruitment,
    get_harvest_rate = get_harvest_rate,
    do_dynamics_length = do_dynamics_length,
    get_length_like = get_length_like,
    get_weight_like = get_weight_like,
    rebin_counts = rebin_counts,
    rebin_matrix = rebin_matrix,
    get_cpue_like = get_cpue_like,
    get_cpue_like_length = get_cpue_like_length,
    get_recruitment_prior = get_recruitment_prior, 
    evaluate_priors = evaluate_priors)
}

#' The opal model
#' 
#' Obtain the negative log-likelihood (NLL) value from the opal model.
#' 
#' @param parameters a \code{list} of parameter values.
#' @param data a \code{list} of data inputs.
#' @return the negative log-likelihood (NLL) value.
#' @importFrom RTMB ADoverload getAll REPORT ADREPORT
#' @export
#' 
opal_model <- function(parameters, data) {
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  getAll(data, parameters, warn = FALSE)
  has_init_rdev_a <- exists("init_rdev_a", inherits = FALSE) && !is.null(init_rdev_a)
  has_init_bias_adj_a <- exists("init_bias_adj_a", inherits = FALSE) && !is.null(init_bias_adj_a)
  if (!exists("cpue_switch", inherits = FALSE)) cpue_switch <- 0L
  if (!exists("lf_switch", inherits = FALSE)) lf_switch <- 0L
  if (!exists("n_lf", inherits = FALSE)) n_lf <- 0L
  if (!exists("lf_addtocomp", inherits = FALSE)) lf_addtocomp <- 1e-08
  if (!exists("wf_switch", inherits = FALSE)) wf_switch <- 0L
  if (!exists("n_wf", inherits = FALSE)) n_wf <- 0L
  if (!exists("wf_addtocomp", inherits = FALSE)) wf_addtocomp <- 1e-08
  if (!exists("log_init_F_f", inherits = FALSE)) log_init_F_f <- rep(log(1e-8), n_fishery)
  if (!has_init_rdev_a) init_rdev_a <- NULL
  if (!exists("bias_adj_y", inherits = FALSE)) bias_adj_y <- rep(1.0, n_year)
  if (!has_init_bias_adj_a) init_bias_adj_a <- NULL
  if (!exists("sex_ratio", inherits = FALSE)) sex_ratio <- rep(1.0, n_age)
  # Length-dynamics switch and companions (missing => age-only path, unchanged)
  if (!exists("length_dynamics", inherits = FALSE)) length_dynamics <- 0L
  if (!exists("plus_group_growth", inherits = FALSE)) plus_group_growth <- 1L
  if (!exists("growth_sd_mode", inherits = FALSE)) growth_sd_mode <- "conditional"
  if (!exists("sd_floor", inherits = FALSE)) sd_floor <- 1e-3

  # Growth module ----

  # Back-transform growth/variability parameters
  L1 <- exp(log_L1)
  L2 <- exp(log_L2)

  # Module 1: Mean length-at-age (Schnute VB)
  mu_a <- get_growth(n_age, A1, A2, L1, L2, log_k)

  # Module 2: SD of length-at-age (linear CV interpolation)
  sd_a <- get_sd_at_age(mu_a, L1, L2, log_CV1, log_CV2)

  # Length-bin geometry: derive from bin scalars if not already supplied
  # (prep_lf_data() adds these; derive here so standalone runs work too).
  if (!exists("len_lower", inherits = FALSE)) {
    len_lower <- seq(from = len_bin_start, by = len_bin_width, length.out = n_len)
  }
  if (!exists("len_upper", inherits = FALSE)) len_upper <- len_lower + len_bin_width
  if (!exists("len_mid", inherits = FALSE)) len_mid <- len_lower + len_bin_width / 2

  # Shared PLA — computed once and reused for weight, maturity, selectivity
  pla <- get_pla(len_lower, len_upper, mu_a, sd_a)

  # ============================================================================
  # Length-dynamics engine (opt-in) ----
  # Joint age-length state advanced by a growth transition matrix. Biology and
  # selectivity stay on the length grid (no PLA collapse); LF/WF likelihoods
  # consume the length marginal directly (identity PLA); CPUE uses NL_yl.
  # ============================================================================
  if (length_dynamics == 1L) {
    B0 <- exp(log_B0)
    h <- exp(log_h)
    sigma_r <- exp(log_sigma_r)
    init_F_f <- exp(log_init_F_f)
    if (!exists("log_lf_tau", inherits = FALSE)) log_lf_tau <- rep(0.0, n_fishery)
    if (!exists("log_wf_tau", inherits = FALSE)) log_wf_tau <- rep(0.0, n_fishery)

    # Biology at length ----
    if (exists("weight", inherits = FALSE) && length(weight) == n_len) {
      weight_l <- weight
    } else {
      weight_l <- get_weight_at_length(len_mid, lw_a, lw_b)
    }
    if (length(maturity) != n_len)
      stop("length_dynamics requires 'maturity' on a length basis (length n_len).")
    if (length(fecundity) != n_len)
      stop("length_dynamics requires 'fecundity' on a length basis (length n_len).")
    maturity_l  <- maturity
    fecundity_l <- fecundity
    sex_ratio_l <- if (length(sex_ratio) == n_len) sex_ratio else rep(1.0, n_len)
    if (exists("spawning_potential", inherits = FALSE) && length(spawning_potential) == n_len) {
      spawning_potential_l <- spawning_potential
    } else {
      spawning_potential_l <- sex_ratio_l * maturity_l * fecundity_l
    }
    # Natural mortality basis (age vector n_age, or length vector n_len)
    if (length(M) == n_len && n_len != n_age) {
      M_basis <- "length"; M_vec <- M
    } else {
      M_basis <- "age"; M_vec <- M
    }

    # Selectivity at length ----
    if (exists("sel_fa_external", inherits = FALSE) && !is.null(sel_fa_external)) {
      sel_dims <- dim(sel_fa_external)
      if (is.null(sel_dims)) {
        stop("'sel_fa_external' must be a matrix or 3-D array.")
      } else if (length(sel_dims) == 2L) {
        if (!all(sel_dims == c(n_fishery, n_len)))
          stop("length-mode 'sel_fa_external' must have dimensions [n_fishery, n_len].")
        sel_fyl <- array(0, dim = c(n_fishery, n_year, n_len))
        for (f in seq_len(n_fishery)) for (y in seq_len(n_year)) sel_fyl[f, y, ] <- sel_fa_external[f, ]
      } else if (length(sel_dims) == 3L) {
        if (!all(sel_dims == c(n_fishery, n_year, n_len)))
          stop("length-mode 'sel_fa_external' must have dimensions [n_fishery, n_year, n_len].")
        sel_fyl <- sel_fa_external
      } else {
        stop("'sel_fa_external' must be a matrix or 3-D array.")
      }
    } else {
      sel_fyl <- get_selectivity_length(data, par_sel, len_mid)
    }

    # Growth transition matrix and recruit length distribution ----
    G <- get_growth_matrix(len_lower, len_upper, len_mid, sd_a, L1, L2, log_k, A1, A2,
                           plus_group_growth = plus_group_growth,
                           sd_mode = growth_sd_mode, sd_floor = sd_floor)
    recruit_dist_l <- get_recruit_length_dist(len_lower, len_upper, mu_a[1], sd_a[1])

    # Weight-at-length by fishery/year
    weight_fyl <- array(0, dim = c(n_fishery, n_year, n_len))
    for (f in seq_len(n_fishery)) for (y in seq_len(n_year)) weight_fyl[f, y, ] <- weight_l

    # Initial equilibrium and dynamics ----
    sel_fl0 <- if (n_fishery == 1L) as.vector(sel_fyl[1, 1, ]) else sel_fyl[, 1, ]
    init <- get_initial_numbers_length(B0 = B0, h = h, M_a = M_vec,
                                       spawning_potential_l = spawning_potential_l,
                                       pla = pla, G = G, recruit_dist_l = recruit_dist_l,
                                       init_F_f = init_F_f, sel_fl = sel_fl0,
                                       init_rdev_a = init_rdev_a, sigma_r = sigma_r,
                                       init_bias_adj_a = init_bias_adj_a,
                                       plus_group_growth = plus_group_growth, M_basis = M_basis)
    R0 <- init$R0; alpha <- init$alpha; beta <- init$beta

    dyn <- do_dynamics_length(data, parameters,
                              B0 = B0, R0 = R0, alpha = alpha, beta = beta, h = h,
                              sigma_r = sigma_r, M_a = M_vec,
                              spawning_potential_l = spawning_potential_l,
                              weight_fyl = weight_fyl, recruit_dist_l = recruit_dist_l, G = G,
                              init_number_al = init$Ninit_al, init_number0_al = init$Ninit0_al,
                              sel_fyl = sel_fyl, bias_adj_y = bias_adj_y,
                              plus_group_growth = plus_group_growth, M_basis = M_basis)

    number_ysal <- dyn$number_ysal
    NL_yl <- dyn$NL_yl
    lp_penalty <- dyn$lp_penalty
    catch_pred_fyl <- dyn$catch_pred_fyl

    # Composition source: catch-at-length, or vulnerable numbers for survey fleets.
    comp_pred_fyl <- catch_pred_fyl
    for (f in seq_len(n_fishery)) {
      if (sum(catch_obs_ysf[, , f]) <= 0) {
        for (y in seq_len(n_year)) comp_pred_fyl[f, y, ] <- NL_yl[y, ] * sel_fyl[f, y, ]
      }
    }

    # Priors ----
    lp_rec <- get_recruitment_prior(rdev_y, sigma_r)
    if (has_init_rdev_a) {
      lp_init_rec <- get_recruitment_prior(init_rdev_a, sigma_r)
    } else {
      lp_init_rec <- 0
    }
    if (exists("priors", inherits = FALSE) && !is.null(priors) && length(priors) > 0) {
      lp_prior <- evaluate_priors(parameters, priors)
    } else {
      lp_prior <- 0
    }

    # Likelihoods (identity PLA => length marginal consumed directly) ----
    pla_identity <- diag(n_len)
    if (cpue_switch > 0) {
      lp_cpue <- get_cpue_like_length(cpue_data, parameters, NL_yl, sel_fyl, weight_l, cpue_switch)
    } else {
      lp_cpue <- 0
    }
    if (lf_switch > 0 && n_lf > 0) {
      if (!exists("lf_year_fi", inherits = FALSE)) lf_year_fi <- split(lf_year, lf_fishery)
      if (!exists("lf_n_fi", inherits = FALSE)) lf_n_fi <- split(lf_n, lf_fishery)
      lp_lf <- get_length_like(
        lf_obs_flat = lf_obs_flat, lf_obs_ints = lf_obs_ints, lf_obs_prop = lf_obs_prop,
        catch_pred_fya = comp_pred_fyl, pla = pla_identity,
        lf_n_f = lf_n_f, lf_fishery_f = lf_fishery_f, lf_year_fi = lf_year_fi,
        lf_n_fi = lf_n_fi, lf_minbin = lf_minbin, lf_maxbin = lf_maxbin,
        removal_switch_f = removal_switch_f, lf_switch = lf_switch,
        n_len = n_len, n_lf = n_lf, log_lf_tau = log_lf_tau, lf_addtocomp = lf_addtocomp)
    } else {
      lp_lf <- 0
    }
    if (wf_switch > 0 && n_wf > 0) {
      if (!exists("wf_year_fi", inherits = FALSE)) wf_year_fi <- split(wf_year, wf_fishery)
      if (!exists("wf_n_fi", inherits = FALSE)) wf_n_fi <- split(wf_n, wf_fishery)
      lp_wf <- get_weight_like(
        wf_obs_flat = wf_obs_flat, wf_obs_ints = wf_obs_ints, wf_obs_prop = wf_obs_prop,
        catch_pred_fya = comp_pred_fyl, pla = pla_identity, wf_rebin_matrix = wf_rebin_matrix,
        wf_n_f = wf_n_f, wf_fishery_f = wf_fishery_f, wf_year_fi = wf_year_fi,
        wf_n_fi = wf_n_fi, wf_minbin = wf_minbin, wf_maxbin = wf_maxbin,
        removal_switch_f = removal_switch_f, wf_switch = wf_switch,
        n_wt = n_wt, n_wf = n_wf, log_wf_tau = log_wf_tau, wf_addtocomp = wf_addtocomp)
    } else {
      lp_wf <- 0
    }

    nll <- lp_prior + lp_penalty + lp_rec + lp_init_rec + sum(lp_cpue) + sum(lp_lf) + sum(lp_wf)

    REPORT(sel_fyl)
    REPORT(weight_l)
    REPORT(weight_fyl)
    REPORT(maturity_l)
    REPORT(fecundity_l)
    REPORT(spawning_potential_l)
    REPORT(recruit_dist_l)
    REPORT(G)
    REPORT(comp_pred_fyl)
    REPORT(pla)
    REPORT(lp_prior)
    REPORT(lp_penalty)
    REPORT(lp_rec)
    REPORT(lp_init_rec)
    REPORT(lp_cpue)
    REPORT(lp_lf)
    REPORT(lp_wf)
    REPORT(B0)
    REPORT(R0)
    REPORT(alpha)
    REPORT(beta)
    REPORT(sigma_r)
    REPORT(init_F_f)
    return(nll)
  }

  # Module 3: Weight-at-age ----
  # If a pre-computed weight vector is supplied in the data , use it directly via resolve_bio_vector.
  # This allows age-basis vectors (length n_age) to pass through unchanged,
  # or length-basis vectors (length n_len) to be converted via the PLA.
  # Otherwise, derive weight-at-age internally from the length-weight
  # relationship and the PLA.
  if (exists("weight", inherits = FALSE)) {
    weight_a <- resolve_bio_vector(weight, n_age, n_len, pla, "weight")
  } else {
    wt_at_len <- get_weight_at_length(len_mid, lw_a, lw_b)
    weight_a <- c(t(pla) %*% wt_at_len)
  }

  # Replicate weight across fisheries and years (AD-safe: use loop + [<- overload)
  weight_fya_mod <- array(0, dim = c(n_fishery, n_year, n_age))
  for (f in seq_len(n_fishery)) {
    for (y in seq_len(n_year)) {
      weight_fya_mod[f, y, ] <- weight_a
    }
  }

  # Module 4: Resolve biology vectors to age-basis via PLA ----
  # Accepts either age-basis (length n_age) or length-basis (length n_len)
  # vectors. Length-basis vectors are converted using: vec_a = t(pla) %*% vec_l
  maturity_a <- resolve_bio_vector(maturity, n_age, n_len, pla, "maturity")
  M_a <- resolve_bio_vector(M, n_age, n_len, pla, "M")
  fecundity_a <- resolve_bio_vector(fecundity, n_age, n_len, pla, "fecundity")
  sex_ratio_a <- resolve_bio_vector(sex_ratio, n_age, n_len, pla, "sex_ratio")

  # Spawning potential-at-age ----
  # If a pre-computed spawning_potential vector is supplied in the data (e.g.,
  # platoon-weighted and sex-ratio-adjusted average from SS3), use it directly
  # via resolve_bio_vector.
  #
  # Otherwise, compute spawning potential from its components:
  #   spawning_potential_a = sex_ratio_a * maturity_a * fecundity_a
  # where sex_ratio is the fraction female at age (typically 0.5 for all ages).
  # This requires a sex_ratio vector in the data object.
  if (exists("spawning_potential", inherits = FALSE)) {
    spawning_potential_a <- resolve_bio_vector(spawning_potential, n_age, n_len, pla, "spawning_potential")
  } else {
    spawning_potential_a <- sex_ratio_a * maturity_a * fecundity_a
  }
  
  # Selectivity ----
  if (exists("sel_fa_external", inherits = FALSE) && !is.null(sel_fa_external)) {
    sel_dims <- dim(sel_fa_external)
    if (is.null(sel_dims)) {
      stop("'sel_fa_external' must be a matrix or 3-D array.")
    } else if (length(sel_dims) == 2L) {
      if (!all(dim(sel_fa_external) == c(n_fishery, n_age))) {
        stop("'sel_fa_external' must have dimensions [n_fishery, n_age].")
      }
      sel_fya <- array(0, dim = c(n_fishery, n_year, n_age))
      for (f in seq_len(n_fishery)) {
        for (y in seq_len(n_year)) {
          sel_fya[f, y, ] <- sel_fa_external[f, ]
        }
      }
    } else if (length(sel_dims) == 3L) {
      if (!all(sel_dims == c(n_fishery, n_year, n_age))) {
        stop("'sel_fa_external' must have dimensions [n_fishery, n_year, n_age].")
      }
      sel_fya <- sel_fa_external
    } else {
      stop("'sel_fa_external' must be a matrix or 3-D array.")
    }
  } else {
    # mu_a and sd_a from growth module so AD gradients propagate if growth is estimated
    sel_fya <- get_selectivity(data, par_sel, pla, len_mid)
  }
  
  # Main population loop ----
  
  B0 <- exp(log_B0)
  h <- exp(log_h)
  sigma_r <- exp(log_sigma_r)
  init_F_f <- exp(log_init_F_f)
  init <- get_initial_numbers(B0 = B0, h = h, M_a = M_a, spawning_potential_a = spawning_potential_a,
                              init_F_f = init_F_f, sel_fa = sel_fya[, 1, ],
                              init_rdev_a = init_rdev_a, sigma_r = sigma_r,
                              init_bias_adj_a = init_bias_adj_a)
  R0 <- init$R0
  alpha <- init$alpha
  beta <- init$beta
  
  dyn <- do_dynamics(data, parameters,
                     B0 = B0, R0 = R0, alpha = alpha, beta = beta, h = h, sigma_r = sigma_r,
                     M_a = M_a, spawning_potential_a = spawning_potential_a, weight_fya = weight_fya_mod,
                     init_number_a = init$Ninit, init_number0_a = init$Ninit0,
                     sel_fya = sel_fya, bias_adj_y = bias_adj_y)
  
  number_ysa <- dyn$number_ysa
  lp_penalty <- dyn$lp_penalty
  catch_pred_fya <- dyn$catch_pred_fya
  comp_pred_fya <- catch_pred_fya
  for (f in seq_len(n_fishery)) {
    if (sum(catch_obs_ysf[, , f]) <= 0) {
      for (y in seq_len(n_year)) {
        comp_pred_fya[f, y, ] <- number_ysa[y, 1, ] * sel_fya[f, y, ]
      }
    }
  }
  
  # plot(spawning_biomass_y)
  # plot(rowSums(dyn$number_ysa[,1,]))
  # plot(catch_obs_ysf - catch_pred_ysf)
  # points(catch_pred_ysf, pch = 2, col = 2)
  
  # Priors ----
  
  lp_rec <- get_recruitment_prior(rdev_y, sigma_r)
  if (has_init_rdev_a) {
    lp_init_rec <- get_recruitment_prior(init_rdev_a, sigma_r)
  } else {
    lp_init_rec <- 0
  }
  if (exists("priors", inherits = FALSE) && !is.null(priors) && length(priors) > 0) {
    lp_prior <- evaluate_priors(parameters, priors)
  } else {
    lp_prior <- 0
  }
  
  # Likelihoods ----
  
  # CPUE likelihood ----
  if (cpue_switch > 0) {
    lp_cpue <- get_cpue_like(cpue_data, parameters, number_ysa, sel_fya, weight_fya_mod,cpue_switch)
  } else {
    lp_cpue <- 0
  }
  
  # Length composition likelihood ----
  if (lf_switch > 0 && n_lf > 0) {
    if (!exists("lf_year_fi", inherits = FALSE)) lf_year_fi <- split(lf_year, lf_fishery)
    if (!exists("lf_n_fi", inherits = FALSE)) lf_n_fi <- split(lf_n, lf_fishery)
    lp_lf <- get_length_like(
      lf_obs_flat = lf_obs_flat,
      lf_obs_ints = lf_obs_ints,
      lf_obs_prop = lf_obs_prop,
      catch_pred_fya = comp_pred_fya,
      pla = pla,
      lf_n_f = lf_n_f,
      lf_fishery_f = lf_fishery_f,
      lf_year_fi = lf_year_fi,
      lf_n_fi = lf_n_fi,
      lf_minbin = lf_minbin,
      lf_maxbin = lf_maxbin,
      removal_switch_f = removal_switch_f,
      lf_switch = lf_switch,
      n_len = n_len,
      n_lf = n_lf, log_lf_tau = log_lf_tau,
      lf_addtocomp = lf_addtocomp
    )
  } else {
    lp_lf <- 0
  }
  # Weight composition likelihood ----
  if (wf_switch > 0 && n_wf > 0) {
    if (!exists("wf_year_fi", inherits = FALSE)) wf_year_fi <- split(wf_year, wf_fishery)
    if (!exists("wf_n_fi", inherits = FALSE)) wf_n_fi <- split(wf_n, wf_fishery)
    lp_wf <- get_weight_like(
      wf_obs_flat = wf_obs_flat,
      wf_obs_ints = wf_obs_ints,
      wf_obs_prop = wf_obs_prop,
      catch_pred_fya = comp_pred_fya,
      pla = pla,
      wf_rebin_matrix = wf_rebin_matrix,
      wf_n_f = wf_n_f,
      wf_fishery_f = wf_fishery_f,
      wf_year_fi = wf_year_fi,
      wf_n_fi = wf_n_fi,
      wf_minbin = wf_minbin,
      wf_maxbin = wf_maxbin,
      removal_switch_f = removal_switch_f,
      wf_switch = wf_switch,
      n_wt = n_wt,
      n_wf = n_wf,
      log_wf_tau = log_wf_tau,
      wf_addtocomp = wf_addtocomp
    )
  } else {
    lp_wf <- 0
  }
  nll <- lp_prior + lp_penalty + lp_rec + lp_init_rec + sum(lp_cpue) + sum(lp_lf) + sum(lp_wf)
  
  # Reporting ----
  
  REPORT(number_ysa)
  REPORT(sel_fya)
  
  REPORT(lp_prior)
  REPORT(lp_penalty)
  REPORT(lp_rec)
  REPORT(lp_init_rec)
  REPORT(lp_cpue)
  REPORT(lp_lf)
  REPORT(lp_wf)
  
  REPORT(B0)
  REPORT(R0)
  REPORT(alpha)
  REPORT(beta)
  REPORT(sigma_r)
  REPORT(maturity_a)
  REPORT(fecundity_a)
  REPORT(spawning_potential_a)
  REPORT(M_a)
  REPORT(weight_fya_mod)
  REPORT(init_F_f)
  if (!has_init_rdev_a) init_rdev_a <- rep(0.0, n_age)
  REPORT(init_rdev_a)
  
  return(nll)
}
