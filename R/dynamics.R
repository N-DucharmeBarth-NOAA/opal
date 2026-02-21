#' Initial numbers and Beverton-Holt parameters
#'
#' Computes the initial equilibrium numbers-at-age, unfished recruitment (R0), and Beverton-Holt stock-recruitment parameters.
#'
#' @param B0 Unfished spawning biomass.
#' @param h Beverton-Holt steepness parameter.
#' @param M_a a \code{vector} of natural mortality at age.
#' @param spawning_potential_a a \code{vector} of spawning potential at age
#'   (maturity × fecundity).
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
get_initial_numbers <- function(B0, h, M_a, spawning_potential_a) {
  "[<-" <- ADoverload("[<-")
  n_age <- length(M_a)
  rel_N <- numeric(n_age)
  rel_N[1] <- 1
  for (a in 2:n_age) rel_N[a] <- rel_N[a - 1] * exp(-M_a[a - 1])
  rel_N[n_age] <- rel_N[n_age] / (1 - exp(-M_a[n_age]))
  R0 <- B0 / sum(spawning_potential_a * rel_N)
  alpha <- (4 * h * R0) / (5 * h - 1)
  beta  <- (B0 * (1 - h)) / (5 * h - 1)
  return(list(Ninit = R0 * rel_N, R0 = R0, alpha = alpha, beta = beta))
}

#' Population dynamics
#'
#' Runs the core age- and season-structured population dynamics loop for bigeye
#' tuna. Starts from initial equilibrium numbers (derived from B0 and h),
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
#'   potential at age (maturity × fecundity). Passed explicitly so AD gradients
#'   propagate if growth is ever estimated.
#' @param weight_fya Numeric array \code{[n_fishery, n_year, n_age]}. Mean
#'   weight at age by fishery and year.  Passed explicitly so AD gradients
#'   propagate if growth is ever estimated.
#' @param init_number_a Numeric vector of length \code{n_age}. Initial
#'   equilibrium numbers-at-age (from \code{\link{get_initial_numbers}}).
#' @param sel_fya Numeric array \code{[n_fishery, n_year, n_age]}.
#'   Fishery-specific selectivity at age by year (from
#'   \code{\link{get_selectivity}}).
#' @return A named list with:
#' \describe{
#'   \item{number_ysa}{Numbers-at-age array \code{[n_year+1, n_season, n_age]}.}
#'   \item{lp_penalty}{Total penalty from \code{\link{posfun}} (harvest rate constraints).}
#' }
#' @importFrom RTMB ADoverload
#' @export
#'
do_dynamics <- function(data, parameters,
                        B0, R0, alpha, beta, h = 0.95, sigma_r = 0.6,
                        M_a, spawning_potential_a, weight_fya,
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
  spawning_biomass_y[1] <- sum(number_ysa[1, 1,] * spawning_potential_a)
  hrate_ysa  <- array(0, dim = c(n_year + 1, n_season, n_age))
  hrate_ysfa  <- array(0, dim = c(n_year + 1, n_season, n_fishery, n_age))
  catch_pred_fya <- array(0, dim = c(n_fishery, n_year, n_age))
  catch_pred_ysf <- array(0, dim = c(n_year, n_season, n_fishery))
  lp_penalty <- 0
  
  for (y in seq_len(n_year)) {
    for (s in seq_len(n_season1)) {
      if (y >= fy) {
        hr <- get_harvest_rate(data, y, s, number_ysa, sel_fya, weight_fya)
        hrate_ysfa[y, s,,] <- hr$h_rate_fa
        hrate_ysa[y, s,] <- hr$h_rate_a
        lp_penalty <- lp_penalty + hr$penalty
      }
      number_ysa[y, s + 1,] <- number_ysa[y, s,] * (1 - hrate_ysa[y, s,]) * S_a
    }
    # Last season
    if (y >= fy) {
      hr <- get_harvest_rate(data, y, n_season, number_ysa, sel_fya, weight_fya)
      hrate_ysfa[y, n_season,,] <- hr$h_rate_fa
      hrate_ysa[y, n_season,] <- hr$h_rate_a
      lp_penalty <- lp_penalty + hr$penalty
    }
    number_ysa[y + 1, 1, 2:n_age] <- number_ysa[y, n_season, 1:n_age1] * (1 - hrate_ysa[y, n_season, 1:n_age1]) * S_a[1:n_age1]
    number_ysa[y + 1, 1, n_age] <- number_ysa[y + 1, 1, n_age] + (number_ysa[y, n_season, n_age] * (1 - hrate_ysa[y, n_season, n_age]) * S_a[n_age])
    spawning_biomass_y[y + 1] <- sum(number_ysa[y + 1, 1,] * spawning_potential_a)
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
  h_rate_a <- colSums(h_rate_fa)
  return(list(h_rate_fa = h_rate_fa, h_rate_a = h_rate_a, penalty = tmp$penalty))
}
