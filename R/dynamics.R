#' Compute unfished equilibrium quantities from natural mortality only
#'
#' Internal helper used by \code{\link{get_initial_numbers}} to calculate
#' unfished survivorship-per-recruit, \eqn{R_0}, and Beverton-Holt parameters.
#'
#' @param B0 Unfished spawning biomass.
#' @param h Beverton-Holt steepness parameter.
#' @param M_a a \code{vector} of natural mortality at age.
#' @param spawning_potential_a a \code{vector} of spawning potential at age
#'   (maturity x fecundity).
#' @return A list with \code{rel_N}, \code{R0}, \code{alpha}, and \code{beta}.
#' @importFrom RTMB ADoverload
#' @keywords internal
#' 
get_unfished_init <- function(B0, h, M_a, spawning_potential_a) {
  "[<-" <- ADoverload("[<-")
  n_age <- length(M_a)
  rel_N <- numeric(n_age)
  rel_N[1] <- 1
  if (n_age > 1) {
    for (a in 2:n_age) rel_N[a] <- rel_N[a - 1] * exp(-M_a[a - 1])
  }
  rel_N[n_age] <- rel_N[n_age] / (1 - exp(-M_a[n_age]))
  R0    <- B0 / sum(spawning_potential_a * rel_N)
  alpha <- (4 * h * R0) / (5 * h - 1)
  beta  <- (B0 * (1 - h)) / (5 * h - 1)
  return(list(rel_N = rel_N, R0 = R0, alpha = alpha, beta = beta))
}

#' Initial numbers and Beverton-Holt parameters
#'
#' Computes the initial equilibrium numbers-at-age, unfished recruitment (R0),
#' and Beverton-Holt stock-recruitment parameters.
#'
#' @param B0 Unfished spawning biomass.
#' @param h Beverton-Holt steepness parameter.
#' @param M_a a \code{vector} of natural mortality at age.
#' @param spawning_potential_a a \code{vector} of spawning potential at age
#'   (maturity x fecundity).
#' @param init_F_f an optional \code{vector} of initial fishing mortality by
#'   fishery.
#' @param sel_fa an optional matrix of selectivity-at-age with dimensions
#'   \code{[n_fishery, n_age]}.
#' @param init_rdev_a an optional \code{vector} of initial age deviations.
#' @param sigma_r recruitment standard deviation used in lognormal correction.
#' @param init_bias_adj_a an optional \code{vector} of bias adjustment scalars
#'   for initial age deviations. Defaults to zero so fixed zero initial
#'   deviations do not alter the equilibrium initial age structure.
#' @return A list containing:
#' \describe{
#'   \item{Ninit}{Initial numbers-at-age (vector).}
#'   \item{Ninit0}{Initial unfished numbers-at-age (vector).}
#'   \item{R0}{Unfished recruitment (scalar).}
#'   \item{alpha}{BH alpha parameter.}
#'   \item{beta}{BH beta parameter.}
#' }
#' @importFrom RTMB ADoverload
#' @export
#'
get_initial_numbers <- function(B0, h, M_a, spawning_potential_a,
                                init_F_f = NULL, sel_fa = NULL,
                                init_rdev_a = NULL, sigma_r = 0.6,
                                init_bias_adj_a = NULL) {
  "[<-" <- ADoverload("[<-")
  n_age <- length(M_a)

  # Unfished survivorship for R0, alpha, beta
  rel_N0 <- numeric(n_age) + B0 * 0
  rel_N0[1] <- 1
  if (n_age > 1) {
    for (a in 2:n_age) rel_N0[a] <- rel_N0[a - 1] * exp(-M_a[a - 1])
  }
  rel_N0[n_age] <- rel_N0[n_age] / (1 - exp(-M_a[n_age]))

  SPR0  <- sum(spawning_potential_a * rel_N0)
  R0    <- B0 / SPR0
  alpha <- (4 * h * R0) / (5 * h - 1)
  beta  <- (B0 * (1 - h)) / (5 * h - 1)

  # Fished survivorship for Ninit
  Z_a <- M_a + B0 * 0
  if (!is.null(init_F_f) && !is.null(sel_fa)) {
      if (is.null(dim(sel_fa))) {
        ## single fishery sel_fa may be passed as vector, not matrix
        Z_a <- Z_a + init_F_f[1L] * sel_fa
      } else {
        for (f in seq_along(init_F_f)) {
          Z_a <- Z_a + init_F_f[f] * sel_fa[f, ]
        }
      }
  }

  rel_N <- numeric(n_age) + B0 * 0
  rel_N[1] <- 1
  if (n_age > 1) {
    for (a in 2:n_age) rel_N[a] <- rel_N[a - 1] * exp(-Z_a[a - 1])
  }
  rel_N[n_age] <- rel_N[n_age] / (1 - exp(-Z_a[n_age]))

  # Fished equilibrium recruitment
  SPR_eq <- sum(spawning_potential_a * rel_N)
  R_eq   <- alpha - (beta / SPR_eq)

  Ninit <- R_eq * rel_N
  Ninit0 <- R0 * rel_N0

  if (!is.null(init_rdev_a)) {
    if (is.null(init_bias_adj_a)) init_bias_adj_a <- rep(0.0, n_age)
    for (a in seq_len(n_age)) {
      Ninit[a] <- Ninit[a] * exp(init_rdev_a[a] - init_bias_adj_a[a] * 0.5 * sigma_r^2)
      Ninit0[a] <- Ninit0[a] * exp(init_rdev_a[a] - init_bias_adj_a[a] * 0.5 * sigma_r^2)
    }
  }

  return(list(Ninit = Ninit, Ninit0 = Ninit0, R0 = R0, alpha = alpha, beta = beta))
}

#' Initial numbers-at-age-and-length (length engine)
#'
#' Length-structured analogue of \code{\link{get_initial_numbers}}, returning the
#' equilibrium numbers-at-age-and-length \code{[n_age, n_len]}.
#'
#' With \code{M_basis = "age"} (default) the initial state is seeded from the
#' static age-length key: the age totals, \code{R0} and Beverton-Holt parameters
#' come straight from \code{\link{get_initial_numbers}} (spawning potential and
#' selectivity collapsed to age via \code{pla}), then each age's numbers are
#' distributed across length by \code{pla[, a]}. This makes the length engine's
#' equilibrium reduce \emph{exactly} to the age-only equilibrium and avoids the
#' discretisation drift of iterating a near-degenerate distribution through the
#' growth matrix. (At a fished initial equilibrium with strong size-selective
#' \code{init_F}, distributing by \code{pla} does not distort the within-age
#' length composition for the removed sizes; this is a small first-timestep
#' approximation the dynamics immediately begin correcting. It is exact when
#' \code{init_F = 0}, as for the opakapaka prototype.)
#'
#' With \code{M_basis = "length"} natural mortality is length-specific, so
#' survivorship is path-dependent and the equilibrium is built by propagating
#' \code{recruit_dist_l} forward through the growth transition matrix \code{G}
#' (subject to bin-discretisation error, worst for coarse bins / old ages).
#'
#' @param B0 Unfished spawning biomass.
#' @param h Beverton-Holt steepness.
#' @param M_a Natural mortality. Length \code{n_age} (age basis, default) or
#'   \code{n_len} (length basis, set \code{M_basis = "length"}).
#' @param spawning_potential_l Numeric vector (length \code{n_len}) of spawning
#'   potential at length (maturity x fecundity x sex ratio).
#' @param pla Numeric matrix \code{[n_len, n_age]} age-length key (columns sum to
#'   1), from \code{\link{get_pla}}.
#' @param G Numeric array \code{[n_age, n_len, n_len]} growth transition matrix
#'   (from \code{\link{get_growth_matrix}}). Used only when \code{M_basis = "length"}.
#' @param recruit_dist_l Numeric vector (length \code{n_len}) recruit length
#'   distribution (from \code{\link{get_recruit_length_dist}}), sums to 1. Used
#'   only when \code{M_basis = "length"}.
#' @param init_F_f Optional numeric vector of initial F by fishery.
#' @param sel_fl Optional selectivity-at-length, matrix \code{[n_fishery, n_len]}
#'   or vector \code{[n_len]} for a single fishery.
#' @param init_rdev_a Optional numeric vector (length \code{n_age}) of initial
#'   age deviations.
#' @param sigma_r Recruitment SD used in the lognormal bias correction.
#' @param init_bias_adj_a Optional numeric vector (length \code{n_age}) of bias
#'   adjustment scalars. Defaults to zero.
#' @param plus_group_growth Logical (or 0/1). Plus-group behaviour for the
#'   \code{M_basis = "length"} path: \code{TRUE} closes the equilibrium with an
#'   \code{(I - S G_plus)^{-1}} solve; \code{FALSE} closes it elementwise.
#' @param M_basis Either \code{"age"} (default) or \code{"length"}.
#' @return A list with \code{Ninit_al}, \code{Ninit0_al} (both \code{[n_age, n_len]}),
#'   \code{R0}, \code{alpha}, \code{beta}.
#' @importFrom RTMB ADoverload colSums
#' @export
#'
get_initial_numbers_length <- function(B0, h, M_a, spawning_potential_l, pla, G,
                                       recruit_dist_l = NULL,
                                       init_F_f = NULL, sel_fl = NULL,
                                       init_rdev_a = NULL, sigma_r = 0.6,
                                       init_bias_adj_a = NULL,
                                       plus_group_growth = TRUE,
                                       M_basis = c("age", "length")) {
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  M_basis <- match.arg(M_basis)
  n_len <- nrow(pla)
  n_age <- ncol(pla)

  if (M_basis == "age") {
    # PLA-seed: exact age-only equilibrium, distributed across length by the ALK.
    sp_a <- as.vector(t(pla) %*% spawning_potential_l)
    sel_fa <- NULL
    if (!is.null(sel_fl)) {
      if (is.null(dim(sel_fl))) {
        sel_fa <- as.vector(t(pla) %*% sel_fl)
      } else {
        sel_fa <- array(0, dim = c(nrow(sel_fl), n_age))
        for (f in seq_len(nrow(sel_fl))) sel_fa[f, ] <- as.vector(t(pla) %*% sel_fl[f, ])
      }
    }
    age <- get_initial_numbers(B0 = B0, h = h, M_a = M_a, spawning_potential_a = sp_a,
                               init_F_f = init_F_f, sel_fa = sel_fa,
                               init_rdev_a = init_rdev_a, sigma_r = sigma_r,
                               init_bias_adj_a = init_bias_adj_a)
    Ninit_al  <- array(0, dim = c(n_age, n_len))
    Ninit0_al <- array(0, dim = c(n_age, n_len))
    for (a in seq_len(n_age)) {
      Ninit_al[a, ]  <- age$Ninit[a]  * pla[, a]
      Ninit0_al[a, ] <- age$Ninit0[a] * pla[, a]
    }
    return(list(Ninit_al = Ninit_al, Ninit0_al = Ninit0_al,
                R0 = age$R0, alpha = age$alpha, beta = age$beta))
  }

  # M_basis == "length": propagate the recruit distribution through G.
  pg <- as.logical(plus_group_growth)

  M_al <- array(0, dim = c(n_age, n_len))
  for (a in seq_len(n_age)) M_al[a, ] <- M_a + B0 * 0        # M_a holds M_l (length n_len)

  F_l <- numeric(n_len) + B0 * 0
  if (!is.null(init_F_f) && !is.null(sel_fl)) {
    if (is.null(dim(sel_fl))) {
      F_l <- F_l + init_F_f[1L] * sel_fl
    } else {
      for (f in seq_along(init_F_f)) F_l <- F_l + init_F_f[f] * sel_fl[f, ]
    }
  }

  build_nphi <- function(Z_al) {
    nphi <- array(0, dim = c(n_age, n_len))
    nphi[1, ] <- recruit_dist_l
    if (n_age > 2) {
      for (a in 2:(n_age - 1)) {
        surv <- exp(-Z_al[a - 1, ]) * nphi[a - 1, ]
        nphi[a, ] <- as.vector(G[a - 1, , ] %*% surv)
      }
    }
    inflow <- as.vector(G[n_age - 1, , ] %*% (exp(-Z_al[n_age - 1, ]) * nphi[n_age - 1, ]))
    s_plus <- exp(-Z_al[n_age, ])
    if (pg) {
      # Plus-group residents survive then grow via G_plus each year:
      #   x = inflow + G_plus %*% (s_plus * x). Solve by fixed-point iteration
      #   (AD-safe: matmul + elementwise only; converges as s_plus < 1).
      Gp <- G[n_age, , ]
      x <- inflow
      for (iter in seq_len(300)) x <- inflow + as.vector(Gp %*% (s_plus * x))
      nphi[n_age, ] <- x
    } else {
      nphi[n_age, ] <- inflow / (1 - s_plus)
    }
    nphi
  }

  ZF <- array(0, dim = c(n_age, n_len))
  for (a in seq_len(n_age)) ZF[a, ] <- M_al[a, ] + F_l

  nphi0 <- build_nphi(M_al)
  nphiF <- build_nphi(ZF)

  SPR0   <- sum(colSums(nphi0) * spawning_potential_l)
  R0     <- B0 / SPR0
  alpha  <- (4 * h * R0) / (5 * h - 1)
  beta   <- (B0 * (1 - h)) / (5 * h - 1)
  SPR_eq <- sum(colSums(nphiF) * spawning_potential_l)
  R_eq   <- alpha - (beta / SPR_eq)

  Ninit_al  <- R_eq * nphiF
  Ninit0_al <- R0 * nphi0

  if (!is.null(init_rdev_a)) {
    if (is.null(init_bias_adj_a)) init_bias_adj_a <- rep(0.0, n_age)
    for (a in seq_len(n_age)) {
      fac <- exp(init_rdev_a[a] - init_bias_adj_a[a] * 0.5 * sigma_r^2)
      Ninit_al[a, ]  <- Ninit_al[a, ] * fac
      Ninit0_al[a, ] <- Ninit0_al[a, ] * fac
    }
  }

  return(list(Ninit_al = Ninit_al, Ninit0_al = Ninit0_al,
              R0 = R0, alpha = alpha, beta = beta))
}

#' Population dynamics
#'
#' Runs the core age- and season-structured population dynamics loop. Starts
#' from initial equilibrium numbers (derived from B0 and h),
#' applies seasonal harvest, natural mortality, spawning, and recruitment
#' (Beverton-Holt with log-normal deviates), and computes predicted catches and
#' harvest rates.
#'
#' All derived biology arrays (\code{M_a}, \code{spawning_potential_a},
#' \code{weight_fya}) are passed as explicit arguments rather than read from
#' \code{data}. This ensures that AD gradients propagate correctly if any of
#' these quantities carry estimated parameters in the future (e.g., growth
#' parameters estimated via the PLA, or natural mortality via the Lorenzen
#' equation).
#'
#' @param data A \code{list} of model data.  Must contain at minimum:
#'   \code{first_yr}, \code{first_yr_catch}, \code{n_year}, \code{n_season},
#'   \code{n_fishery}, \code{n_age}, \code{catch_obs_ysf},
#'   \code{catch_units_f}.
#' @param parameters A \code{list} of model parameters.  Must contain at
#'   minimum: \code{rdev_y}.
#' @param B0 Numeric. Unfished equilibrium spawning biomass.
#' @param R0 Numeric. Unfished equilibrium recruitment.
#' @param alpha Numeric. Beverton-Holt stock-recruitment alpha parameter.
#' @param beta Numeric. Beverton-Holt stock-recruitment beta parameter.
#' @param h Numeric (0.2–1). Steepness of the Beverton-Holt
#'   stock-recruitment relationship.
#' @param sigma_r Numeric > 0. Standard deviation of log recruitment
#'   deviations.
#' @param M_a Numeric vector of length \code{n_age}. Natural mortality at age.
#'   Passed explicitly so AD gradients propagate if M is ever estimated.
#' @param spawning_potential_a Numeric vector of length \code{n_age}. Spawning
#'   potential at age (maturity x fecundity). Passed explicitly so AD gradients
#'   propagate if growth is ever estimated.
#' @param weight_fya Numeric array \code{[n_fishery, n_year, n_age]}. Mean
#'   weight at age by fishery and year.  Passed explicitly so AD gradients
#'   propagate if growth is ever estimated.
#' @param init_number_a Numeric vector of length \code{n_age}. Initial
#'   equilibrium numbers-at-age (from \code{\link{get_initial_numbers}}).
#' @param init_number0_a Numeric vector of length \code{n_age}. Initial
#'   unfished equilibrium numbers-at-age (from \code{\link{get_initial_numbers}}).
#' @param sel_fya Numeric array \code{[n_fishery, n_year, n_age]}.
#'   Fishery-specific selectivity at age by year (from
#'   \code{\link{get_selectivity}}).
#' @param bias_adj_y Numeric vector of length \code{n_year}. Recruitment bias
#'   adjustment scalar by year.
#' @return A named list with:
#' \describe{
#'   \item{number_ysa}{Numbers-at-age array \code{[n_year+1, n_season, n_age]}.}
#'   \item{number0_ysa}{Unfished numbers-at-age array \code{[n_year+1, n_season, n_age]}.}
#'   \item{lp_penalty}{Total penalty from \code{\link{posfun}} (harvest rate constraints).}
#'   \item{catch_pred_fya}{Predicted catch-at-age array \code{[n_fishery, n_year, n_age]}.}
#'   \item{spawning_biomass_y}{Spawning biomass trajectory under fishing.}
#'   \item{spawning_biomass0_y}{Spawning biomass trajectory in the dynamic unfished state.}
#'   \item{static_depletion_y}{Static depletion trajectory \code{spawning_biomass_y / B0}.}
#'   \item{dynamic_depletion_y}{Dynamic depletion trajectory \code{spawning_biomass_y / spawning_biomass0_y}.}
#' }
#' @importFrom RTMB ADoverload
#' @export
#'
do_dynamics <- function(data, parameters,
                        B0, R0, alpha, beta, h = 0.95, sigma_r = 0.6,
                        M_a, spawning_potential_a, weight_fya,
                        init_number_a, init_number0_a, sel_fya, bias_adj_y = NULL) {
  
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  getAll(data, parameters, warn = FALSE)
  if (is.null(bias_adj_y)) bias_adj_y <- rep(1.0, n_year)
  fy <- first_yr_catch - first_yr + 1
  n_age1 <- n_age - 1
  S_a <- exp(-M_a / n_season)
  number_ysa <- array(0, dim = c(n_year + 1, n_season, n_age))
  number_ysa[1, 1,] <- init_number_a
  spawning_biomass_y <- numeric(n_year + 1)
  spawning_biomass_y[1] <- sum(number_ysa[1, 1,] * spawning_potential_a)
  number0_ysa <- array(0, dim = c(n_year + 1, n_season, n_age))
  number0_ysa[1, 1,] <- init_number0_a
  spawning_biomass0_y <- numeric(n_year + 1)
  spawning_biomass0_y[1] <- sum(number0_ysa[1, 1,] * spawning_potential_a)
  hrate_ysa  <- array(0, dim = c(n_year + 1, n_season, n_age))
  hrate_ysfa  <- array(0, dim = c(n_year + 1, n_season, n_fishery, n_age))
  catch_pred_fya <- array(0, dim = c(n_fishery, n_year, n_age))
  catch_pred_ysf <- array(0, dim = c(n_year, n_season, n_fishery))
  lp_penalty <- 0
  eps_denom <- 1e-6
  F_f <- numeric(n_fishery)
  h_rate_fa <- array(0, dim = c(n_fishery, n_age))
  f_weight <- which(catch_units_f == 1)
  f_numbers <- which(catch_units_f != 1)
  
  for (y in seq_len(n_year)) {
    for (s in seq_len(n_season)) {
      if (y >= fy) {
        F_f[] <- 0
        h_rate_fa[] <- 0
        N_ys <- number_ysa[y, s,]

        for (f in f_weight) {
          if (catch_obs_ysf[y, s, f] > 0) {
            sel_N <- N_ys * sel_fya[f, y,]
            Nsum <- sum(sel_N * weight_fya[f, y,]) + eps_denom
            F_f[f] <- catch_obs_ysf[y, s, f] / Nsum
            h_rate_fa[f,] <- F_f[f] * sel_fya[f, y,]
          }
        }
        for (f in f_numbers) {
          if (catch_obs_ysf[y, s, f] > 0) {
            sel_N <- N_ys * sel_fya[f, y,]
            Nsum <- sum(sel_N) + eps_denom
            F_f[f] <- catch_obs_ysf[y, s, f] / Nsum
            h_rate_fa[f,] <- F_f[f] * sel_fya[f, y,]
          }
        }

        sum_F <- sum(F_f)
        tmp <- posfun(x = 1 - sum_F, eps = 0.001)
        lp_penalty <- lp_penalty + tmp$penalty
        hrate_ysfa[y, s,,] <- h_rate_fa
        hrate_ysa[y, s,] <- colSums(h_rate_fa)

        for (f in f_weight) {
          if (catch_obs_ysf[y, s, f] > 0) {
            catch_at_age_fs <- h_rate_fa[f,] * N_ys
            catch_pred_fya[f, y,] <- catch_pred_fya[f, y,] + catch_at_age_fs
            catch_pred_ysf[y, s, f] <- sum(catch_at_age_fs * weight_fya[f, y,])
          }
        }
        for (f in f_numbers) {
          if (catch_obs_ysf[y, s, f] > 0) {
            catch_at_age_fs <- h_rate_fa[f,] * N_ys
            catch_pred_fya[f, y,] <- catch_pred_fya[f, y,] + catch_at_age_fs
            catch_pred_ysf[y, s, f] <- sum(catch_at_age_fs)
          }
        }
      }
      if (s < n_season) {
        number_ysa[y, s + 1,] <- number_ysa[y, s,] * (1 - hrate_ysa[y, s,]) * S_a
        number0_ysa[y, s + 1,] <- number0_ysa[y, s,] * S_a
      }
    }

    number_ysa[y + 1, 1, 2:n_age] <- number_ysa[y, n_season, 1:n_age1] * (1 - hrate_ysa[y, n_season, 1:n_age1]) * S_a[1:n_age1]
    number_ysa[y + 1, 1, n_age] <- number_ysa[y + 1, 1, n_age] + (number_ysa[y, n_season, n_age] * (1 - hrate_ysa[y, n_season, n_age]) * S_a[n_age])
    spawning_biomass_y[y + 1] <- sum(number_ysa[y + 1, 1,] * spawning_potential_a)
    number0_ysa[y + 1, 1, 2:n_age] <- number0_ysa[y, n_season, 1:n_age1] * S_a[1:n_age1]
    number0_ysa[y + 1, 1, n_age] <- number0_ysa[y + 1, 1, n_age] + (number0_ysa[y, n_season, n_age] * S_a[n_age])
    spawning_biomass0_y[y + 1] <- sum(number0_ysa[y + 1, 1,] * spawning_potential_a)

    number_ysa[y + 1, 1, 1] <- get_recruitment(sbio = spawning_biomass_y[y + 1], rdev = rdev_y[y], B0 = B0, alpha = alpha, beta = beta, sigma_r = sigma_r, bias_adj = bias_adj_y[y])
    number0_ysa[y + 1, 1, 1] <- get_recruitment(sbio = spawning_biomass0_y[y + 1], rdev = rdev_y[y], B0 = B0, alpha = alpha, beta = beta, sigma_r = sigma_r, bias_adj = bias_adj_y[y])
  }
  static_depletion_y <- spawning_biomass_y / B0
  dynamic_depletion_y <- spawning_biomass_y / spawning_biomass0_y
  
  REPORT(catch_pred_ysf)
  REPORT(catch_pred_fya)
  REPORT(hrate_ysa)
  REPORT(hrate_ysfa)
  REPORT(number0_ysa)
  REPORT(spawning_biomass_y)
  REPORT(spawning_biomass0_y)
  REPORT(static_depletion_y)
  REPORT(dynamic_depletion_y)
  ADREPORT(spawning_biomass_y)
  ADREPORT(spawning_biomass0_y)
  ADREPORT(static_depletion_y)
  ADREPORT(dynamic_depletion_y)
  
  return(list(number_ysa = number_ysa, number0_ysa = number0_ysa, lp_penalty = lp_penalty,
              catch_pred_fya = catch_pred_fya,
              spawning_biomass_y = spawning_biomass_y,
              spawning_biomass0_y = spawning_biomass0_y,
              static_depletion_y = static_depletion_y,
              dynamic_depletion_y = dynamic_depletion_y))
}

#' Population dynamics (joint age-length engine)
#'
#' Length-structured analogue of \code{\link{do_dynamics}}. Carries a joint
#' numbers-at-age-and-length state \code{[n_year+1, n_season, n_age, n_len]} and
#' advances it each year by applying length-specific harvest and (age- or
#' length-based) natural mortality in the current length bin, then ageing and
#' redistributing across length via the growth transition matrix \code{G} (from
#' \code{\link{get_growth_matrix}}). Recruits enter age 1 spread over length by
#' \code{recruit_dist_l}. Because mortality is size-selective, the length-at-age
#' distribution evolves through time (unlike the static-PLA age-only engine).
#'
#' Removals use opal's Pope-style harvest-rate formulation (as in
#' \code{do_dynamics}), applied at length: for each fishery,
#' \eqn{F_f = C_f / \sum_l v_l w_l} (weight) or \eqn{C_f / \sum_l v_l} (numbers)
#' with vulnerable numbers \eqn{v_l = NL_l \, s_{f,l}}, and harvest rate
#' \eqn{h_{f,l} = F_f s_{f,l}}. Age and length marginals are produced for
#' downstream/reporting compatibility.
#'
#' @param data,parameters Model data and parameters (as for \code{do_dynamics}).
#' @param B0,R0,alpha,beta,h,sigma_r Population/stock-recruit scalars.
#' @param M_a Natural mortality: length \code{n_age} (age basis) or \code{n_len}
#'   (length basis; set \code{M_basis = "length"}).
#' @param spawning_potential_l Numeric vector (length \code{n_len}) spawning
#'   potential at length.
#' @param weight_fyl Numeric array \code{[n_fishery, n_year, n_len]} weight at
#'   length by fishery and year.
#' @param recruit_dist_l Numeric vector (length \code{n_len}) recruit length
#'   distribution (sums to 1).
#' @param G Numeric array \code{[n_age, n_len, n_len]} growth transition matrix.
#' @param init_number_al,init_number0_al Numeric matrices \code{[n_age, n_len]}
#'   initial (fished / unfished) equilibrium numbers, from
#'   \code{\link{get_initial_numbers_length}}.
#' @param sel_fyl Numeric array \code{[n_fishery, n_year, n_len]} selectivity at
#'   length by fishery and year.
#' @param bias_adj_y Numeric vector (length \code{n_year}) recruitment bias
#'   adjustment. Defaults to ones.
#' @param plus_group_growth Logical (or 0/1). If \code{TRUE} the plus-age-group
#'   keeps growing via \code{G[n_age,,]}; if \code{FALSE} its length is frozen.
#' @param M_basis Either \code{"age"} (default) or \code{"length"}.
#' @return A named list with \code{number_ysal}, \code{number0_ysal},
#'   \code{number_ysa}, \code{number0_ysa} (age marginals), \code{NL_yl} (length
#'   marginal, start-of-year), \code{lp_penalty}, \code{catch_pred_fyl},
#'   \code{catch_pred_fya}, \code{catch_pred_ysf}, \code{spawning_biomass_y},
#'   \code{spawning_biomass0_y}, \code{dynamic_depletion_y}.
#' @importFrom RTMB ADoverload getAll REPORT ADREPORT colSums rowSums
#' @export
#'
do_dynamics_length <- function(data, parameters,
                               B0, R0, alpha, beta, h = 0.95, sigma_r = 0.6,
                               M_a, spawning_potential_l, weight_fyl,
                               recruit_dist_l, G,
                               init_number_al, init_number0_al, sel_fyl,
                               bias_adj_y = NULL, plus_group_growth = TRUE,
                               M_basis = c("age", "length")) {
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  getAll(data, parameters, warn = FALSE)
  M_basis <- match.arg(M_basis)
  if (is.null(bias_adj_y)) bias_adj_y <- rep(1.0, n_year)
  fy <- first_yr_catch - first_yr + 1
  n_age1 <- n_age - 1
  n_len <- length(recruit_dist_l)

  # Seasonal survival multiplier per (age, length).
  S_al <- array(0, dim = c(n_age, n_len))
  if (M_basis == "age") {
    for (a in seq_len(n_age)) S_al[a, ] <- exp(-M_a[a] / n_season)
  } else {
    for (a in seq_len(n_age)) S_al[a, ] <- exp(-M_a / n_season)   # M_a holds M_l
  }

  number_ysal  <- array(0, dim = c(n_year + 1, n_season, n_age, n_len))
  number0_ysal <- array(0, dim = c(n_year + 1, n_season, n_age, n_len))
  number_ysal[1, 1, , ]  <- init_number_al
  number0_ysal[1, 1, , ] <- init_number0_al

  NL_yl <- array(0, dim = c(n_year + 1, n_len))
  NL_yl[1, ] <- colSums(init_number_al)

  spawning_biomass_y  <- numeric(n_year + 1)
  spawning_biomass0_y <- numeric(n_year + 1)
  spawning_biomass_y[1]  <- sum(colSums(init_number_al) * spawning_potential_l)
  spawning_biomass0_y[1] <- sum(colSums(init_number0_al) * spawning_potential_l)

  catch_pred_fyl <- array(0, dim = c(n_fishery, n_year, n_len))
  catch_pred_fya <- array(0, dim = c(n_fishery, n_year, n_age))
  catch_pred_ysf <- array(0, dim = c(n_year, n_season, n_fishery))
  hrate_ysl <- array(0, dim = c(n_year, n_season, n_len))

  lp_penalty <- 0
  eps_denom <- 1e-6
  f_weight  <- which(catch_units_f == 1)
  f_numbers <- which(catch_units_f != 1)

  for (y in seq_len(n_year)) {
    for (s in seq_len(n_season)) {
      F_f <- numeric(n_fishery)
      h_rate_fl <- array(0, dim = c(n_fishery, n_len))
      if (y >= fy) {
        N_al <- number_ysal[y, s, , ]
        NL   <- colSums(N_al)
        for (f in f_weight) {
          if (catch_obs_ysf[y, s, f] > 0) {
            vul_l <- NL * sel_fyl[f, y, ]
            Nsum  <- sum(vul_l * weight_fyl[f, y, ]) + eps_denom
            F_f[f] <- catch_obs_ysf[y, s, f] / Nsum
            h_rate_fl[f, ] <- F_f[f] * sel_fyl[f, y, ]
          }
        }
        for (f in f_numbers) {
          if (catch_obs_ysf[y, s, f] > 0) {
            vul_l <- NL * sel_fyl[f, y, ]
            Nsum  <- sum(vul_l) + eps_denom
            F_f[f] <- catch_obs_ysf[y, s, f] / Nsum
            h_rate_fl[f, ] <- F_f[f] * sel_fyl[f, y, ]
          }
        }
        sum_F <- sum(F_f)
        tmp <- posfun(x = 1 - sum_F, eps = 0.001)
        lp_penalty <- lp_penalty + tmp$penalty
        for (f in f_weight) {
          if (catch_obs_ysf[y, s, f] > 0) {
            catch_l <- h_rate_fl[f, ] * NL
            catch_pred_fyl[f, y, ] <- catch_pred_fyl[f, y, ] + catch_l
            catch_pred_fya[f, y, ] <- catch_pred_fya[f, y, ] + as.vector(N_al %*% h_rate_fl[f, ])
            catch_pred_ysf[y, s, f] <- sum(catch_l * weight_fyl[f, y, ])
          }
        }
        for (f in f_numbers) {
          if (catch_obs_ysf[y, s, f] > 0) {
            catch_l <- h_rate_fl[f, ] * NL
            catch_pred_fyl[f, y, ] <- catch_pred_fyl[f, y, ] + catch_l
            catch_pred_fya[f, y, ] <- catch_pred_fya[f, y, ] + as.vector(N_al %*% h_rate_fl[f, ])
            catch_pred_ysf[y, s, f] <- sum(catch_l)
          }
        }
      }
      hrate_l <- colSums(h_rate_fl)
      hrate_ysl[y, s, ] <- hrate_l
      if (s < n_season) {
        Nnext  <- array(0, dim = c(n_age, n_len))
        Nnext0 <- array(0, dim = c(n_age, n_len))
        for (a in seq_len(n_age)) {
          Nnext[a, ]  <- number_ysal[y, s, a, ] * (1 - hrate_l) * S_al[a, ]
          Nnext0[a, ] <- number0_ysal[y, s, a, ] * S_al[a, ]
        }
        number_ysal[y, s + 1, , ]  <- Nnext
        number0_ysal[y, s + 1, , ] <- Nnext0
      }
    }

    # Year boundary: last-season harvest + M, then age up and grow via G.
    hrate_last <- hrate_ysl[y, n_season, ]
    surv_al  <- array(0, dim = c(n_age, n_len))
    surv0_al <- array(0, dim = c(n_age, n_len))
    for (a in seq_len(n_age)) {
      surv_al[a, ]  <- number_ysal[y, n_season, a, ] * (1 - hrate_last) * S_al[a, ]
      surv0_al[a, ] <- number0_ysal[y, n_season, a, ] * S_al[a, ]
    }
    if (n_age1 >= 2) {
      for (a in 2:n_age1) {
        number_ysal[y + 1, 1, a, ]  <- as.vector(G[a - 1, , ] %*% surv_al[a - 1, ])
        number0_ysal[y + 1, 1, a, ] <- as.vector(G[a - 1, , ] %*% surv0_al[a - 1, ])
      }
    }
    number_ysal[y + 1, 1, n_age, ]  <- as.vector(G[n_age1, , ] %*% surv_al[n_age1, ]) +
      as.vector(G[n_age, , ] %*% surv_al[n_age, ])
    number0_ysal[y + 1, 1, n_age, ] <- as.vector(G[n_age1, , ] %*% surv0_al[n_age1, ]) +
      as.vector(G[n_age, , ] %*% surv0_al[n_age, ])

    NLnext  <- colSums(number_ysal[y + 1, 1, , ])
    NL0next <- colSums(number0_ysal[y + 1, 1, , ])
    spawning_biomass_y[y + 1]  <- sum(NLnext  * spawning_potential_l)
    spawning_biomass0_y[y + 1] <- sum(NL0next * spawning_potential_l)

    rec  <- get_recruitment(sbio = spawning_biomass_y[y + 1], rdev = rdev_y[y], B0 = B0, alpha = alpha, beta = beta, sigma_r = sigma_r, bias_adj = bias_adj_y[y])
    rec0 <- get_recruitment(sbio = spawning_biomass0_y[y + 1], rdev = rdev_y[y], B0 = B0, alpha = alpha, beta = beta, sigma_r = sigma_r, bias_adj = bias_adj_y[y])
    number_ysal[y + 1, 1, 1, ]  <- rec  * recruit_dist_l
    number0_ysal[y + 1, 1, 1, ] <- rec0 * recruit_dist_l
    NL_yl[y + 1, ] <- colSums(number_ysal[y + 1, 1, , ])
  }

  # Age marginals for downstream/reporting compatibility.
  number_ysa  <- array(0, dim = c(n_year + 1, n_season, n_age))
  number0_ysa <- array(0, dim = c(n_year + 1, n_season, n_age))
  for (y in seq_len(n_year + 1)) {
    for (s in seq_len(n_season)) {
      number_ysa[y, s, ]  <- rowSums(number_ysal[y, s, , ])
      number0_ysa[y, s, ] <- rowSums(number0_ysal[y, s, , ])
    }
  }
  dynamic_depletion_y <- spawning_biomass_y / spawning_biomass0_y

  REPORT(number_ysal)
  REPORT(number0_ysal)
  REPORT(NL_yl)
  REPORT(catch_pred_fyl)
  REPORT(catch_pred_fya)
  REPORT(catch_pred_ysf)
  REPORT(hrate_ysl)
  REPORT(spawning_biomass_y)
  REPORT(spawning_biomass0_y)
  REPORT(dynamic_depletion_y)
  RTMB::ADREPORT(spawning_biomass_y)
  RTMB::ADREPORT(spawning_biomass0_y)
  RTMB::ADREPORT(dynamic_depletion_y)

  return(list(number_ysal = number_ysal, number0_ysal = number0_ysal,
              number_ysa = number_ysa, number0_ysa = number0_ysa,
              NL_yl = NL_yl, lp_penalty = lp_penalty,
              catch_pred_fyl = catch_pred_fyl, catch_pred_fya = catch_pred_fya,
              catch_pred_ysf = catch_pred_ysf,
              spawning_biomass_y = spawning_biomass_y,
              spawning_biomass0_y = spawning_biomass0_y,
              dynamic_depletion_y = dynamic_depletion_y))
}

#' Harvest rate calculation
#'
#' Computes age-specific harvest rates by fishery for a single year-season
#' combination, using the Baranov catch equation.
#'
#' \code{weight_fya} is passed as an explicit argument (not read from
#' \code{data}) so that AD gradients propagate correctly if growth parameters
#' are estimated in the future.
#'
#' @param data A \code{list} of model data.  Must contain: \code{n_fishery},
#'   \code{n_age}, \code{catch_obs_ysf}, \code{catch_units_f}.
#' @param y Integer. Year index (1-based).
#' @param s Integer. Season index (1-based).
#' @param number_ysa Numeric array \code{[n_year+1, n_season, n_age]}.
#'   Current numbers-at-age.
#' @param sel_fya Numeric array \code{[n_fishery, n_year, n_age]}.
#'   Selectivity at age by fishery and year.
#' @param weight_fya Numeric array \code{[n_fishery, n_year, n_age]}. Mean
#'   weight at age by fishery and year.  Passed explicitly so AD gradients
#'   propagate if growth is ever estimated.
#' @return A named list with:
#' \describe{
#'   \item{h_rate_fa}{Harvest rate array \code{[n_fishery, n_age]}.}
#'   \item{h_rate_a}{Total harvest rate vector (length \code{n_age}).}
#'   \item{penalty}{Penalty from \code{\link{posfun}} for harvest-rate constraint.}
#' }
#' @importFrom RTMB ADoverload colSums
#' @export
#'
get_harvest_rate <- function(data, y, s, number_ysa, sel_fya, weight_fya) {
  "[<-" <- ADoverload("[<-")
  n_fishery <- data$n_fishery
  n_age <- data$n_age
  catch_obs_ysf <- data$catch_obs_ysf
  catch_units_f <- data$catch_units_f
  eps_denom <- 1e-6
  F_f <- numeric(n_fishery)
  h_rate_fa <- array(0, dim = c(n_fishery, n_age))
  for (f in seq_len(n_fishery)) {
    if (catch_obs_ysf[y, s, f] > 0) {
      if (catch_units_f[f] == 1) { # weight
        Nsum <- sum(number_ysa[y, s,] * sel_fya[f, y,] * weight_fya[f, y,]) + eps_denom
      } else if (catch_units_f[f] == 2) { # numbers
        Nsum <- sum(number_ysa[y, s,] * sel_fya[f, y,]) + eps_denom
      }
      F_f[f] <- catch_obs_ysf[y, s, f] / Nsum
      h_rate_fa[f,] <- F_f[f] * sel_fya[f, y,]
    }
  }
  sum_F <- sum(F_f)
  tmp <- posfun(x = 1 - sum_F, eps = 0.001)
  h_rate_a <- colSums(h_rate_fa)
  return(list(h_rate_fa = h_rate_fa, h_rate_a = h_rate_a, penalty = tmp$penalty))
}
