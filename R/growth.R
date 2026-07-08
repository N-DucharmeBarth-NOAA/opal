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
#' @param log_k Numeric. VB growth coefficient (may be AD).
#' @param min_age Integer. Minimum age (default 1L).
#' @return Numeric vector of length \code{n_age}: mean length at each age
#'   \code{a = min_age, ..., min_age + n_age - 1}.
#' @export
get_growth <- function(n_age, A1, A2, L1, L2, log_k, min_age = 1L) {
  ages <- seq(min_age, by = 1, length.out = n_age)
  k   <- exp(log_k)
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

#' Recruit length distribution
#'
#' Distribution of recruits across length bins, assuming length-at-recruitment
#' is normal with mean \code{mu_r} and SD \code{sd_r}. Matches the first column
#' of \code{\link{get_pla}} when \code{mu_r = mu_a[1]} and \code{sd_r = sd_a[1]}.
#' AD-compatible (uses \code{RTMB::pnorm}) so gradients propagate to growth
#' parameters when they are estimated.
#'
#' @param len_lower Numeric vector of lower bounds of length bins (length n_len).
#' @param len_upper Numeric vector of upper bounds of length bins (length n_len).
#' @param mu_r Numeric. Mean length at recruitment (may be AD).
#' @param sd_r Numeric. SD of length at recruitment (may be AD).
#' @return Numeric vector (length n_len) summing to 1.
#' @importFrom RTMB ADoverload pnorm
#' @export
#'
get_recruit_length_dist <- function(len_lower, len_upper, mu_r, sd_r) {
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  n_len <- length(len_lower)
  edges <- c(len_lower, len_upper[n_len])
  p_edges <- pnorm((edges - mu_r) / sd_r)
  d <- p_edges[2:(n_len + 1)] - p_edges[1:n_len]
  d <- d / (sum(d) + 1e-12)
  return(d)
}

#' Growth transition matrix (length-bin to length-bin, per age step)
#'
#' Builds the growth transition operator used by the joint age-length dynamics
#' (\code{\link{do_dynamics_length}}). For each age step \code{a} (age \code{a}
#' to \code{a+1}) it returns a \code{[n_len, n_len]} matrix \code{G[a, dest, src]}
#' giving the probability that a fish in source length bin \code{src} moves to
#' destination bin \code{dest} over one year, discretised from a von Bertalanffy
#' one-step growth increment with normal spread. Columns (over \code{dest}) sum
#' to 1. This is the transition analogue of the static \code{\link{get_pla}}
#' age-length key: it advances an evolving length-at-age distribution rather than
#' assuming a fixed one.
#'
#' The one-step conditional mean is the Ford-Walford recursion
#' \eqn{E[L_{a+1} | L_a = l] = L_\infty (1 - \rho) + \rho l} with
#' \eqn{\rho = \exp(-k)} and \eqn{L_\infty} the asymptote of the Schnute growth
#' curve (plus a small per-age anchor correction so the transition preserves the
#' reference mean-at-age). Each source bin's next-year length is a binned normal
#' about that mean; tails beyond the length range fold into the first and
#' plus-length bins, and the column is normalised to conserve numbers. This is the
#' transition analogue of \code{\link{get_pla}} (the same normal discretisation),
#' so the propagated marginal reproduces the static age-length key under no
#' fishing. (A strict no-shrinkage upper-triangular form was avoided: folding all
#' downward mass onto the diagonal biases the mean upward and that bias compounds
#' over ages.)
#'
#' The transition SD is set so that, propagated forward, the marginal length-at-age
#' SD reproduces \code{sd_a} (the existing linear-CV SD-at-age). With
#' \code{sd_mode = "conditional"} (default) the conditional width is
#' \eqn{\sigma_{cond}[a] = \sqrt{sd_a[a+1]^2 - \rho^2 sd_a[a]^2}} (floored at
#' \code{sd_floor}); the length engine then reduces to \code{get_pla} under no
#' fishing. \code{sd_mode = "marginal"} uses the destination marginal SD directly
#' (faithful to \code{lbm}, but inflates the realised marginal SD).
#'
#' @param len_lower Numeric vector of lower bin bounds (length n_len). Data only.
#' @param len_upper Numeric vector of upper bin bounds (length n_len). Data only.
#' @param len_mid Numeric vector of bin midpoints (length n_len). Data only.
#' @param sd_a Numeric vector of SD of length-at-age (length n_age). May be AD.
#' @param L1 Numeric. Length at reference age A1 (may be AD).
#' @param L2 Numeric. Length at reference age A2 (may be AD).
#' @param log_k Numeric. Log VB growth coefficient (may be AD).
#' @param A1 Integer. Reference age for L1 (data).
#' @param A2 Integer. Reference age for L2 (data).
#' @param plus_group_growth Logical (or 0/1). If \code{TRUE} the plus-age-group
#'   keeps growing toward \eqn{L_\infty} via a self-transition each year; if
#'   \code{FALSE} its length distribution is frozen (identity), matching \code{lbm}.
#' @param sd_mode Either \code{"conditional"} (default) or \code{"marginal"}.
#' @param sd_floor Numeric lower bound on the transition SD. Default \code{1e-3}.
#' @return Numeric array \code{[n_age, n_len, n_len]}. Slices \code{1..n_age-1}
#'   are the age \code{a -> a+1} transitions; slice \code{n_age} is the
#'   plus-group self-transition.
#' @importFrom RTMB ADoverload pnorm
#' @export
#'
get_growth_matrix <- function(len_lower, len_upper, len_mid, sd_a,
                              L1, L2, log_k, A1, A2,
                              plus_group_growth = TRUE,
                              sd_mode = "conditional", sd_floor = 1e-3) {
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  n_len <- length(len_lower)
  n_age <- length(sd_a)
  edges <- c(len_lower, len_upper[n_len])          # n_len + 1 contiguous edges
  k    <- exp(log_k)
  rho  <- exp(-k)
  Linf <- L1 + (L2 - L1) / (1 - exp(-k * (A2 - A1)))   # a -> Inf limit of Schnute mu_a
  pg   <- as.logical(plus_group_growth)

  # Smooth positive floor ~ max(x, sd_floor^2), AD-safe (differentiable).
  floor2 <- sd_floor^2
  smax <- function(x) 0.5 * (x + floor2 + sqrt((x - floor2)^2 + 1e-12))

  # Per-step transition SD. Slot a = a->a+1; slot n_age = plus-group self-step.
  sigma_step <- numeric(n_age)
  if (identical(sd_mode, "conditional")) {
    for (a in seq_len(n_age - 1)) {
      sigma_step[a] <- sqrt(smax(sd_a[a + 1]^2 - rho^2 * sd_a[a]^2))
    }
    sigma_step[n_age] <- sqrt(smax(sd_a[n_age]^2 - rho^2 * sd_a[n_age]^2))
  } else {
    for (a in seq_len(n_age - 1)) sigma_step[a] <- sqrt(smax(sd_a[a + 1]^2))
    sigma_step[n_age] <- sqrt(smax(sd_a[n_age]^2))
  }

  # Per-age mean-preserving anchor. Anchoring the growth increment on raw bin
  # midpoints biases the first moment (a bin's mass is not centred on its
  # midpoint), and the bias accumulates over ages. Add a per-source-age shift
  # delta_a so the transition reproduces the reference mean-at-age exactly:
  # given the age-a marginal (binned normal, mean meanP_a), the output mean is
  # Linf(1-rho) + rho*meanP_a + delta_a, which we set equal to mu_a[a+1].
  mu_a <- get_growth(n_age, A1, A2, L1, L2, log_k)
  meanP <- numeric(n_age)
  for (a in seq_len(n_age)) {
    Fe <- pnorm((edges - mu_a[a]) / sd_a[a])
    p  <- Fe[2:(n_len + 1)] - Fe[1:n_len]
    p  <- p / (sum(p) + 1e-12)
    meanP[a] <- sum(p * len_mid)
  }
  delta <- numeric(n_age)
  for (a in seq_len(n_age - 1)) delta[a] <- mu_a[a + 1] - (Linf * (1 - rho) + rho * meanP[a])
  # plus-group self-step: hold the mean at the plus-group mean (no further ageing target)
  delta[n_age] <- mu_a[n_age] - (Linf * (1 - rho) + rho * meanP[n_age])

  # Each source bin's next-year length is a (mean-anchored) binned normal about
  # Mux[i]; the lower/upper tails beyond the length range fold into the first and
  # plus-length bins, and the column is normalised to conserve numbers. This is
  # the transition analogue of get_pla (same discretisation of a normal), so the
  # propagated marginal reproduces the static age-length key under no fishing. A
  # strict no-shrinkage (upper-triangular) construction was avoided: folding all
  # downward mass onto the diagonal biases the mean upward and that bias
  # compounds over ages (drifting length-at-age away from the PLA).
  build_step <- function(sigma, dlt) {
    Mux <- Linf * (1 - rho) + rho * len_mid + dlt  # mean-preserving source anchor
    Tm <- array(0, dim = c(n_len, n_len))          # [dest, src]
    for (i in seq_len(n_len)) {
      Fc  <- pnorm((edges - Mux[i]) / sigma)       # CDF at each edge, length n_len+1
      col <- Fc[2:(n_len + 1)] - Fc[1:n_len]       # probability in each destination bin
      col[1]     <- col[1] + Fc[1]                 # lower tail -> first bin
      col[n_len] <- col[n_len] + (1 - Fc[n_len + 1]) # upper tail -> plus-length bin
      col <- col / (sum(col) + 1e-12)              # conserve numbers (column sums to 1)
      Tm[, i] <- col
    }
    Tm
  }

  G <- array(0, dim = c(n_age, n_len, n_len))
  for (a in seq_len(n_age - 1)) G[a, , ] <- build_step(sigma_step[a], delta[a])
  if (pg) {
    G[n_age, , ] <- build_step(sigma_step[n_age], delta[n_age])
  } else {
    Ident <- array(0, dim = c(n_len, n_len))
    for (i in seq_len(n_len)) Ident[i, i] <- 1
    G[n_age, , ] <- Ident
  }
  return(G)
}
