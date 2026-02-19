#' Compute mean length-at-age using the Schnute parameterization of VB growth
#'
#' Uses the Schnute parameterization of the von Bertalanffy growth curve,
#' matching the SS3 formulation with \code{CV_Growth_Pattern = 2}.
#'
#' @param n_age Integer. Number of age classes.
#' @param A1 Integer. Reference age for L1 (data).
#' @param A2 Integer. Reference age for L2 (data).
#' @param L1 Numeric. Length at age A1 (may be AD).
#' @param L2 Numeric. Length at age A2 (may be AD).
#' @param k Numeric. VB growth coefficient (may be AD).
#' @return Numeric vector of length \code{n_age}: mean length at each age
#'   \code{a = 1, ..., n_age}.
#' @export
get_growth <- function(n_age, A1, A2, L1, L2, k) {
  ages <- 1:n_age
  mu_a <- L1 + (L2 - L1) * (1 - exp(-k * (ages - A1))) / (1 - exp(-k * (A2 - A1)))
  return(mu_a)
}

#' Compute SD of length-at-age from CV1 and CV2
#'
#' Linearly interpolates CV as a function of mean length, matching
#' SS3's \code{CV_Growth_Pattern = 2}.
#'
#' @param mu_a Numeric vector. Mean length at age (from \code{\link{get_growth}}).
#'   May be AD.
#' @param L1 Numeric. Length at age A1 (may be AD).
#' @param L2 Numeric. Length at age A2 (may be AD).
#' @param CV1 Numeric. CV at age A1 (may be AD).
#' @param CV2 Numeric. CV at age A2 (may be AD).
#' @return Numeric vector of SD at each age.
#' @export
get_sd_at_age <- function(mu_a, L1, L2, CV1, CV2) {
  cv_a <- CV1 + (mu_a - L1) / (L2 - L1) * (CV2 - CV1)
  sd_a <- mu_a * cv_a
  return(sd_a)
}

#' Compute weight at each length bin midpoint
#'
#' @param len_mid Numeric vector. Length bin midpoints (data).
#' @param lw_a Numeric. L-W scalar (data).
#' @param lw_b Numeric. L-W exponent (data).
#' @return Numeric vector of weight at each length bin.
#' @export
get_weight_at_length <- function(len_mid, lw_a, lw_b) {
  wt_at_len <- lw_a * len_mid^lw_b
  return(wt_at_len)
}

#' Convert maturity-at-length to maturity-at-age
#'
#' Uses the probability-of-length-at-age matrix (PLA) to convert a
#' maturity-at-length vector to maturity-at-age:
#' \code{mat_a = t(pla) \%*\% mat_l}
#'
#' @param pla Matrix (n_len x n_age). Probability of length at age (from
#'   \code{\link{get_pla}}).
#' @param maturity_at_length Numeric vector (length n_len). Maturity at each
#'   length bin (data).
#' @return Numeric vector (length n_age). Maturity at each age.
#' @export
get_maturity_at_age <- function(pla, maturity_at_length) {
  maturity_a <- as.vector(t(pla) %*% maturity_at_length)
  return(maturity_a)
}

#' Compute spawning potential at age
#'
#' Spawning potential is defined as maturity × fecundity at each age.
#' Both maturity-at-age and fecundity-at-age are derived from their
#' respective length-based vectors via the PLA (age-length key).
#'
#' @param pla Matrix (n_len x n_age). Probability of length at age (from \code{\link{get_pla}}).
#' @param maturity_at_length Numeric vector (length n_len). Maturity at each length bin (data).
#' @param fecundity_at_length Numeric vector (length n_len). Fecundity at each length bin (data).
#' @return Numeric vector (length n_age). Spawning potential at each age.
#' @export
get_spawning_potential <- function(pla, maturity_at_length, fecundity_at_length) {
  maturity_a <- as.vector(t(pla) %*% maturity_at_length)
  fecundity_a <- as.vector(t(pla) %*% fecundity_at_length)
  spawning_potential_a <- maturity_a * fecundity_a
  return(spawning_potential_a)
}
