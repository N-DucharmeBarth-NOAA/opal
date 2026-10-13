#' Compute mean length-at-age using the Schnute parameterization of VB growth
#'
#' Uses the Schnute parameterization of the von Bertalanffy growth curve.
#' Optionally reproduces Stock Synthesis (SS3) GrowthModel 3 (age-specific K),
#' including its cumulative multipliers and plus-group approximation.
#'
#' @param n_age Integer. Number of age classes.
#' @param A1 Integer. Reference age for L1 (data).
#' @param A2 Integer. Reference age for L2 (data).
#' @param L1 Numeric. Length at age A1 (may be AD).
#' @param L2 Numeric. Length at age A2 (may be AD).
#' @param log_k Numeric. VB growth coefficient (may be AD).
#' @param min_age Integer. Minimum age (default 1L).
#' @param log_k_mult Numeric vector (may be AD). Log K multipliers; \code{NULL}
#'   or length zero retains the single-K curve.
#' @param k_mult_ages Integer vector (data). Ages at which the multipliers apply,
#'   strictly increasing within the modeled ages and matching \code{log_k_mult}.
#' @param plus_group_adjust Logical (data). Apply SS3's oldest-age length
#'   approximation, independently of the multipliers (default \code{FALSE}).
#' @details Multipliers apply cumulatively: K at age a is the base K times
#'   the product of all multipliers at ages less than or equal to a. The final
#'   K carries forward. Linf is derived from L1, L2, A1 and A2 using the
#'   oldest-age K; growth from age a to a+1 uses K at a. Consequently, with
#'   multipliers the curve need not pass through L2 at A2. A1 must be a modeled
#'   age for this option. Below A1, opal retains base-K exponential growth
#'   toward this Linf; SS3 instead uses a linear segment.
#'
#'   The plus-group approximation averages a linear increase from the oldest
#'   length to Linf over a further \code{max_age} ages, with weights
#'   \code{exp(-0.2)^j}, j = 0, ..., max_age. The fixed 0.2 is not mortality;
#'   this option requires a positive integer oldest age.
#'
#'   \code{log_k} is the log rate per model age unit, whereas SS3 reads an
#'   annual rate and scales by season duration. For a quarterly age step,
#'   use \code{log(K_SS3 * 3/12)}; multipliers are dimensionless. This
#'   reimplementation is pinned to SS3 v3.24AB beginning-of-season lengths;
#'   later SS3 versions may differ. CV interpolation is handled separately by
#'   \code{\link{get_sd_at_age}}.
#' @references Methot, R.D. Jr. and Wetzel, C.R. (2013). Stock synthesis:
#'   A biological and statistical framework for fish stock assessment and
#'   fishery management. Fisheries Research 142: 86--99.
#'   \doi{10.1016/j.fishres.2012.10.012}.
#'
#'   \href{https://github.com/nmfs-ost/ss3-source-code/blob/d0da9a5/SS3trial.tpl}{Stock Synthesis source},
#'   sections SS_Label_Info_16.2.2--16.2.4.1.4.
#' @return Numeric vector of length \code{n_age}: mean length at each age
#'   \code{a = min_age, ..., min_age + n_age - 1}.
#' @export
get_growth <- function(n_age, A1, A2, L1, L2, log_k, min_age = 1L,
                       log_k_mult = NULL, k_mult_ages = NULL,
                       plus_group_adjust = FALSE) {
  "[<-" <- ADoverload("[<-")
  ages <- seq(min_age, by = 1, length.out = n_age)
  .validate_growth_options(n_age, A1, min_age, log_k_mult, k_mult_ages,
                           plus_group_adjust)
  k   <- exp(log_k)
  if (length(log_k_mult) == 0L) {
    mu_a <- L1 + (L2 - L1) * (1 - exp(-k * (ages - A1))) / (1 - exp(-k * (A2 - A1)))
    Linf <- L1 + (L2 - L1) / (1 - exp(-k * (A2 - A1)))
  } else {
    # Stock Synthesis GrowthModel 3 (Methot and Wetzel 2013): cumulative K,
    # oldest-age K for Linf, and the SS3 plus-group length approximation.
    k_a <- numeric(n_age) + k * 0
    for (i in seq_len(n_age)) {
      k_a[i] <- k * exp(sum(log_k_mult[k_mult_ages <= ages[i]]))
    }
    Linf <- L1 + (L2 - L1) / (1 - exp(-k_a[n_age] * (A2 - A1)))
    i1 <- match(A1, ages)
    mu_a <- numeric(n_age) + k * 0
    mu_a[i1] <- L1
    if (i1 < n_age) {
      for (i in seq.int(i1 + 1L, n_age)) {
        mu_a[i] <- Linf + (mu_a[i - 1L] - Linf) * exp(-k_a[i - 1L])
      }
    }
    if (i1 > 1L) {
      mu_a[seq_len(i1 - 1L)] <- Linf + (L1 - Linf) * exp(-k * (ages[seq_len(i1 - 1L)] - A1))
    }
  }
  if (plus_group_adjust) {
    max_age <- ages[n_age]
    j <- 0:max_age
    w <- exp(-0.2)^j
    L_n <- mu_a[n_age]
    mu_a[n_age] <- sum(w * (L_n + j / max_age * (Linf - L_n))) / sum(w)
  }
  return(mu_a)
}

.validate_growth_options <- function(n_age, A1, min_age = 1L,
                                     log_k_mult = NULL, k_mult_ages = NULL,
                                     plus_group_adjust = FALSE,
                                     check_length = TRUE) {
  if (check_length && length(log_k_mult) != length(k_mult_ages)) {
    stop("`log_k_mult` and `k_mult_ages` must have the same length.", call. = FALSE)
  }
  max_age <- min_age + n_age - 1L
  if (length(k_mult_ages)) {
    if (!is.numeric(k_mult_ages) || any(!is.finite(k_mult_ages)) ||
        any(k_mult_ages != floor(k_mult_ages)) ||
        any(k_mult_ages < min_age | k_mult_ages > max_age)) {
      stop("`k_mult_ages` must contain integer ages within the modeled age range.", call. = FALSE)
    }
    if (any(diff(k_mult_ages) <= 0)) {
      stop("`k_mult_ages` must be strictly increasing.", call. = FALSE)
    }
    if (length(A1) != 1L || !is.finite(A1) ||
        !A1 %in% seq(min_age, by = 1, length.out = n_age)) {
      stop("`A1` must be a modeled age when K multipliers are supplied.", call. = FALSE)
    }
  }
  if ((!is.logical(plus_group_adjust) && !is.numeric(plus_group_adjust)) ||
      length(plus_group_adjust) != 1L || is.na(plus_group_adjust) ||
      !plus_group_adjust %in% c(FALSE, TRUE)) {
    stop("`plus_group_adjust` must be a logical flag or 0/1.", call. = FALSE)
  }
  if (plus_group_adjust && (max_age <= 0 || max_age != floor(max_age))) {
    stop("Plus-group adjustment requires a positive integer oldest age.", call. = FALSE)
  }
  invisible(NULL)
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
#' @param log_CV1 Numeric. CV at age A1 (may be AD).
#' @param log_CV2 Numeric. CV at age A2 (may be AD).
#' @return Numeric vector of SD at each age.
#' @export
get_sd_at_age <- function(mu_a, L1, L2, log_CV1, log_CV2) {
  CV1 <- exp(log_CV1)
  CV2 <- exp(log_CV2)  
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
