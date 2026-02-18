#' Initial numbers and Beverton-Holt parameters
#'
#' Computes the initial equilibrium numbers-at-age, unfished recruitment (R0), and Beverton-Holt stock-recruitment parameters.
#'
#' @param B0 Unfished spawning biomass.
#' @param h Beverton-Holt steepness parameter.
#' @param M_a a \code{vector} of natural mortality at age.
#' @param maturity_a a \code{vector} of maturity at age.
#' @return A list containing:
#' \describe{
#'   \item{Ninit}{Initial numbers-at-age (vector).}
#'   \item{R0}{Unfished recruitment (scalar).}
#'   \item{alpha}{BH alpha parameter.}
#'   \item{beta}{BH beta parameter.}
#' }
#' @importFrom RTMB ADoverload
#' @export
#'
get_initial_numbers <- function(B0, h, M_a, maturity_a) {
  "[<-" <- ADoverload("[<-")
  n_age <- length(M_a)
  rel_N <- numeric(n_age)
  rel_N[1] <- 1
  for (a in 2:n_age) rel_N[a] <- rel_N[a - 1] * exp(-M_a[a - 1])
  rel_N[n_age] <- rel_N[n_age] / (1 - exp(-M_a[n_age]))
  R0 <- B0 / sum(maturity_a * rel_N)
  alpha <- (4 * h * R0) / (5 * h - 1)
  beta  <- (B0 * (1 - h)) / (5 * h - 1)
  return(list(Ninit = R0 * rel_N, R0 = R0, alpha = alpha, beta = beta))
}

#' Population dynamics
#' 
#' Runs the core age- and season-structured population dynamics loop for Southern Bluefin Tuna.
#' Starts from initial equilibrium numbers (derived from B0 and h), applies seasonal harvest,
#' natural mortality, spawning, recruitment (Beverton-Holt + log-normal deviates), and
#' computes predicted catches and harvest rates. Handles both standard Baranov harvest
#' (removal_switch_f == 0) and direct age-sliced removals (removal_switch_f == 1).
#'
#' @param first_yr Integer. First model year (e.g. 1931).
#' @param first_yr_catch Integer. First year with catch data (determines when fishing begins).
#' @param B0 Numeric. Unfished equilibrium spawning biomass.
#' @param R0 Numeric. Unfished equilibrium recruitment (usually computed inside; kept for backward compatibility).
#' @param alpha Numeric. Beverton-Holt stock-recruitment alpha parameter.
#' @param beta Numeric. Beverton-Holt stock-recruitment beta parameter.
#' @param h Numeric (0.2–1). Steepness of the Beverton-Holt stock-recruitment relationship.
#' @param sigma_r Numeric > 0. Standard deviation of log recruitment deviations.
#' @param rdev_y Numeric vector of length n_year. Log-scale recruitment deviations (one per recruitment year).
#' @param M_a Numeric vector of length n_age. Natural mortality at age (age 0 to max_age).
#' @param phi_ya Numeric matrix (n_year + 1 rows, n_age columns). Spawning output (maturity × fecundity × weight) per recruit at age, by year (pre-computed via [get_phi]).
#' @param init_number_a 3D numeric array (n_season = 2, n_age). Initial numbers-at-age by season. Input provides initial conditions (year 1, season 1).
#' @param removal_switch_f Integer vector of length n_fishery (= 6). 0 = standard Baranov harvest rate, 1 = direct removal using sliced age frequencies (af_sliced_ysfa).
#' @param catch_obs_ysf 3D numeric array (n_catch_years, n_season = 2, n_fishery). Observed catch (weight) by year, season, and fishery.
#' @param sel_fya 3D numeric array (n_fishery, n_year, n_age). Fishery-specific selectivity at age by year (pre-computed via [get_selectivity]).
#' @param weight_fya 3D numeric array (n_fishery, n_year, n_age). Mean weight at age by fishery and year.
#' @return A named list containing:
#' \itemize{
#'   \item number_ysa – updated numbers-at-age array
#'   \item spawning_biomass_y – spawning biomass by year (length n_year + 1)
#'   \item hrate_ysfa – harvest rates by year, season, fishery, age
#'   \item hrate_ysa – total harvest rate by year, season, age
#'   \item catch_pred_fya – predicted catch numbers by fishery, year, age
#'   \item catch_pred_ysf – predicted catch weight by year, season, fishery
#'   \item lp_penalty – total penalty from posfun() (harvest rate constraints)
#'   \item (R0, alpha, beta also returned in some versions)
#' }
#' @importFrom RTMB ADoverload
#' @export
#' 
do_dynamics <- function(data, parameters,
                        B0, R0, alpha, beta, h = 0.95, sigma_r,
                        init_number_a, sel_fya) {
  
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  getAll(data, parameters, warn = FALSE)
  fy <- first_yr_catch - first_yr + 1
  n_age1 <- n_age - 1
  n_season1 <- n_season - 1
  S_a <- exp(-M_a / n_season)
  number_ysa <- array(0, dim = c(n_year + 1, n_season, n_age))
  number_ysa[1, 1,] <- init_number_a
  spawning_biomass_y <- numeric(n_year + 1)
  spawning_biomass_y[1] <- sum(number_ysa[1, 1,] * maturity_a)
  hrate_ysa  <- array(0, dim = c(n_year + 1, n_season, n_age))
  hrate_ysfa  <- array(0, dim = c(n_year + 1, n_season, n_fishery, n_age))
  catch_pred_fya <- array(0, dim = c(n_fishery, n_year, n_age))
  catch_pred_ysf <- array(0, dim = c(n_year, n_season, n_fishery))
  lp_penalty <- 0
  
  for (y in seq_len(n_year)) {
    for (s in seq_len(n_season1)) {
      if (y >= fy) {
        hr <- get_harvest_rate(data, y, s, number_ysa, sel_fya)
        hrate_ysfa[y, s,,] <- hr$h_rate_fa
        hrate_ysa[y, s,] <- hr$h_rate_a
        lp_penalty <- lp_penalty + hr$penalty
      }
      number_ysa[y, s + 1,] <- number_ysa[y, s,] * (1 - hrate_ysa[y, s,]) * S_a
    }
    # Last season
    if (y >= fy) {
      hr <- get_harvest_rate(data, y, n_season, number_ysa, sel_fya)
      hrate_ysfa[y, n_season,,] <- hr$h_rate_fa
      hrate_ysa[y, n_season,] <- hr$h_rate_a
      lp_penalty <- lp_penalty + hr$penalty
    }
    number_ysa[y + 1, 1, 2:n_age] <- number_ysa[y, n_season, 1:n_age1] * (1 - hrate_ysa[y, n_season, 1:n_age1]) * S_a[1:n_age1]
    number_ysa[y + 1, 1, n_age] <- number_ysa[y + 1, 1, n_age] + (number_ysa[y, n_season, n_age] * (1 - hrate_ysa[y, n_season, n_age]) * S_a[n_age])
    spawning_biomass_y[y + 1] <- sum(number_ysa[y + 1, 1,] * maturity_a)
    number_ysa[y + 1, 1, 1] <- get_recruitment(sbio = spawning_biomass_y[y + 1], rdev = rdev_y[y], B0 = B0, alpha = alpha, beta = beta, sigma_r = sigma_r)
    
    for (f in seq_len(n_fishery)) {
      for (s in seq_len(n_season)) {
        catch_pred_fya[f, y,] <- catch_pred_fya[f, y,] + hrate_ysfa[y, s, f,] * number_ysa[y, s,]
        for (a in seq_len(n_age)) {
          if (catch_units_f[f] == 1) { # weight
            catch_pred_ysf[y, s, f] <- catch_pred_ysf[y, s, f] + hrate_ysfa[y, s, f, a] * number_ysa[y, s, a] * weight_fya[f, y, a]
          } else { # numbers 
            catch_pred_ysf[y, s, f] <- catch_pred_ysf[y, s, f] + hrate_ysfa[y, s, f, a] * number_ysa[y, s, a]
          }
        }
      }
    }
  }
  
  REPORT(catch_pred_ysf)
  REPORT(catch_pred_fya)
  REPORT(hrate_ysa)
  REPORT(hrate_ysfa)
  REPORT(spawning_biomass_y)
  
  return(list(number_ysa = number_ysa, lp_penalty = lp_penalty))
}

#' Harvest rate calculation
#'
#' Computes age-specific harvest rates by fishery.
#'
#' @param y,s Year and season index.
#' @param first_yr,first_yr_catch Model and catch start years.
#' @param removal_switch_f Vector of slice-switch flags.
#' @param catch_obs_ysf Observed catch by year-season-fishery.
#' @param number_ysa 3D array of numbers-at-age.
#' @param sel_fya Selectivity array.
#' @param weight_fya Weight-at-age.
#' @return a \code{list} with age-specific harvest rates, fishery, and penalty term.
#' @importFrom RTMB ADoverload colSums
#' @export
#' 
get_harvest_rate <- function(data, y, s, number_ysa, sel_fya) {
  "[<-" <- ADoverload("[<-")
  getAll(data, warn = FALSE)
  F_f <- numeric(n_fishery)
  h_rate_fa <- array(0, dim = c(n_fishery, n_age))
  for (f in seq_len(n_fishery)) {
    if (catch_obs_ysf[y, s, f] > 0) {
      if (catch_units_f[f] == 1) { # weight
        Nsum <- sum(number_ysa[y, s,] * sel_fya[f, y,] * weight_fya[f, y,]) + 1e-6
      } else if (catch_units_f[f] == 2) { # numbers
        Nsum <- sum(number_ysa[y, s,] * sel_fya[f, y,]) + 1e-6
      }
      F_f[f] <- catch_obs_ysf[y, s, f] / Nsum
      h_rate_fa[f,] <- F_f[f] * sel_fya[f, y,]
    }
  }
  sum_F <- sum(F_f)
  tmp <- posfun(x = 1 - sum_F, eps = 0.001)
  # F_f <- F_f / sum_F * tmp$new
  # for (f in seq_len(n_fishery)) {
  #   if (catch_obs_ysf[yy, s, f] > 0) {
  #     if (removal_switch_f[f] == 0) {
  #       h_rate_fa[f,] <- F_f[f] * sel_fya[f,y,] # what about direct removals
  #     }
  #   }
  # }
  h_rate_a <- colSums(h_rate_fa)
  return(list(h_rate_fa = h_rate_fa, h_rate_a = h_rate_a, penalty = tmp$penalty))
}
