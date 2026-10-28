#' Compute mean length-at-age
#'
#' Computes mean length at each modeled age under one of three growth curves,
#' selected by \code{growth_model}: the Schnute parameterization of the von
#' Bertalanffy (VB) curve, the Richards curve, or the growth cessation model.
#'
#' @param n_age Integer. Number of age classes.
#' @param A1 Integer. Reference age for L1 (data).
#' @param A2 Integer. Reference age for L2 (data).
#' @param L1 Numeric. Length at age A1 (may be AD). Used by \code{"vonB"}
#'   and \code{"richards"}.
#' @param L2 Numeric. Length at age A2 (may be AD). Used by \code{"vonB"}
#'   and \code{"richards"}.
#' @param log_k Numeric. Log growth coefficient per model age step (may be
#'   AD). Used by \code{"vonB"} and \code{"richards"}.
#' @param min_age Integer. Minimum age (default 1L).
#' @param log_k_mult Numeric vector (may be AD). Log K multipliers; \code{NULL}
#'   or length zero retains the single-K curve. \code{"vonB"} only.
#' @param k_mult_ages Integer vector (data). Ages at which the multipliers apply,
#'   strictly increasing within the modeled ages and matching \code{log_k_mult}.
#'   \code{"vonB"} only.
#' @param plus_group_adjust Logical (data). Apply SS3's oldest-age length
#'   approximation, independently of the multipliers (default \code{FALSE}).
#'   \code{"vonB"} only.
#' @param growth_model Character (data). Growth curve: \code{"vonB"}
#'   (default), \code{"richards"} or \code{"cessation"}. See Details.
#' @param b Numeric. Richards shape parameter on the natural scale (may be
#'   AD). Must be non-zero. Used by \code{"richards"}.
#' @param L0 Numeric. Length at absolute age 0, natural scale (may be AD).
#'   Used by \code{"cessation"}.
#' @param log_rmax Numeric. Log maximum growth rate, in length per model age
#'   step (may be AD). Used by \code{"cessation"}.
#' @param log_k_gc Numeric. Log steepness of the logistic decline in growth
#'   rate, per model age step (may be AD). Used by \code{"cessation"}.
#' @param t50 Numeric. Absolute age at the logistic midpoint (may be AD).
#'   Used by \code{"cessation"}.
#' @param age_offset Numeric (data). Added to the modeled ages to give the
#'   absolute age \code{t} used by \code{"cessation"} (default 0).
#' @details
#' \subsection{Growth models}{
#'   \tabular{lll}{
#'     \code{growth_model} \tab Curve \tab Arguments used \cr
#'     \code{"vonB"} (default) \tab Schnute VB \tab \code{L1}, \code{L2},
#'       \code{log_k} (optionally \code{log_k_mult}, \code{k_mult_ages},
#'       \code{plus_group_adjust}) \cr
#'     \code{"richards"} \tab Richards \tab \code{L1}, \code{L2},
#'       \code{log_k}, \code{b} \cr
#'     \code{"cessation"} \tab Growth cessation \tab \code{L0},
#'       \code{log_rmax}, \code{log_k_gc}, \code{t50}, \code{age_offset}
#'   }
#'   Stock Synthesis (SS3) GrowthModel 3 (age-specific K) is reproduced by
#'   \code{"vonB"} with \code{log_k_mult} and \code{k_mult_ages}.
#'
#'   All rates (\code{log_k}, \code{log_rmax} and \code{log_k_gc}) are per
#'   model age step.
#' }
#'
#' \subsection{von Bertalanffy (growth_model = "vonB")}{
#'   Uses the Schnute parameterization of the von Bertalanffy growth curve,
#'   \deqn{L_a = L_1 + (L_2 - L_1)\frac{1 - e^{-k(a - A_1)}}{1 - e^{-k(A_2 - A_1)}}.}{L_a = L1 + (L2 - L1) (1 - exp(-k (a - A1))) / (1 - exp(-k (A2 - A1))).}
#'   Optionally reproduces SS3 GrowthModel 3 (age-specific K), including its
#'   cumulative multipliers and plus-group approximation.
#'
#'   Multipliers apply cumulatively: K at age a is the base K times
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
#' }
#'
#' \subsection{Richards (growth_model = "richards")}{
#'   Uses the Richards (1959) curve as parameterized by Schnute (1981),
#'   \deqn{L_a = \left[L_1^b + (L_2^b - L_1^b)\frac{1 - e^{-k(a - A_1)}}{1 - e^{-k(A_2 - A_1)}}\right]^{1/b}.}{L_a = [L1^b + (L2^b - L1^b) (1 - exp(-k (a - A1))) / (1 - exp(-k (A2 - A1)))]^(1/b).}
#'   The parameterization and the guidance in this subsection follow the
#'   Stock Synthesis user manual (section 8.5.2.2.2, Richards growth
#'   function). See the Richards growth vignette on the Stock Synthesis
#'   website for further technical guidance.
#'
#'   \code{b = 1} gives the VB curve. As \code{b} approaches 0, the curve
#'   approaches the Gompertz curve; Gompertz growth can be obtained by fixing
#'   \code{b} at a small value such as \code{1e-4}. \code{b} is a
#'   floating-point parameter on the natural scale that may be positive or
#'   negative but not exactly 0. If its initial value is non-zero, there is
#'   effectively no risk of an estimated \code{b} landing exactly on 0.
#'
#'   When \code{A1} is greater than the youngest modeled age, some parameter
#'   combinations give undefined (\code{NaN}) predicted lengths at the younger
#'   ages. Below \code{A1} the fraction in the equation is negative, so with
#'   \code{b > 0} the bracket can become negative. With \code{b < 0} the
#'   asymptotic bracket can be zero or negative, in which case Linf is
#'   undefined. The choice of \code{A1} and \code{A2} also constrains the
#'   possible curve shapes. opal does not clamp the bracket.
#' }
#'
#' \subsection{Growth cessation (growth_model = "cessation")}{
#'   Uses the growth cessation model of Maunder et al. (2018), Eq. 2,
#'   \deqn{L_t = L_0 + r_{max}\left[\frac{\ln(e^{-k t_{50}} + 1) - \ln(e^{k(t - t_{50})} + 1)}{k} + t\right],}{L_t = L0 + rmax [(log(exp(-k t50) + 1) - log(exp(k (t - t50)) + 1)) / k + t],}
#'   with \code{rmax = exp(log_rmax)} and \code{k = exp(log_k_gc)}. The
#'   growth rate declines logistically from \code{rmax}, with midpoint age
#'   \code{t50}. \code{L(0) = L0}, and the asymptotic length is
#'   \eqn{L_0 + r_{max}[\ln(e^{-k t_{50}} + 1)/k + t_{50}]}{L0 + rmax [log(exp(-k t50) + 1) / k + t50]}.
#'   As \code{k} approaches 0, the curve approaches \code{L0 + rmax * t / 2}.
#'   Both log terms are computed as an overflow-safe \eqn{\log(1 + e^x)}{log(1 + exp(x))}
#'   (\code{RTMB::logspace_add} for AD inputs), which avoids overflow for
#'   large \code{k * (t - t50)}.
#'
#'   Unlike the VB and Richards curves, which depend only on \code{a - A1},
#'   this curve depends on absolute age \code{t}, so the age origin matters.
#'   \code{t} is the modeled age plus \code{age_offset}. For example, if
#'   model age 1 corresponds to biological age 0, use \code{age_offset = -1}.
#' }
#' @references Maunder, M.N., Deriso, R.B., Schaefer, K.M., Fuller, D.W.,
#'   Aires-da-Silva, A.M., Minte-Vera, C.V. and Campana, S.E. (2018). The
#'   growth cessation model: a growth model for species showing a near
#'   cessation in growth with application to bigeye tuna (\emph{Thunnus
#'   obesus}). Marine Biology 165: 76. \doi{10.1007/s00227-018-3336-9}.
#'
#'   Methot, R.D. Jr. and Wetzel, C.R. (2013). Stock synthesis:
#'   A biological and statistical framework for fish stock assessment and
#'   fishery management. Fisheries Research 142: 86--99.
#'   \doi{10.1016/j.fishres.2012.10.012}.
#'
#'   Methot, R.D. Jr., Wetzel, C.R., Taylor, I.G., Doering, K.L., Perl, E.F.
#'   and Johnson, K.F. Stock Synthesis User Manual, Version 3.30.25.1. NOAA
#'   Fisheries, Seattle, WA. Section 8.5.2.2.2, Richards growth function.
#'   \url{https://nmfs-ost.github.io/ss3-doc/SS330_User_Manual_release.html}.
#'
#'   Richards, F.J. (1959). A flexible growth function for empirical use.
#'   Journal of Experimental Botany 10: 290--301.
#'
#'   Schnute, J. (1981). A versatile growth model with statistically stable
#'   parameters. Canadian Journal of Fisheries and Aquatic Sciences 38:
#'   1128--1140.
#'
#'   \href{https://github.com/nmfs-ost/ss3-source-code/blob/d0da9a5/SS3trial.tpl}{Stock Synthesis source},
#'   sections SS_Label_Info_16.2.2--16.2.4.1.4.
#' @return Numeric vector of length \code{n_age}: mean length at each age
#'   \code{a = min_age, ..., min_age + n_age - 1}.
#' @examples
#' \dontrun{
#'   # Starting values from the IO yellowfin inputs (quarterly ages, SS3
#'   # GrowthModel 3 solution). In opal_model() the same options come from
#'   # data (growth_model, k_mult_ages, growth_plus_group_adjust,
#'   # growth_age_offset) and parameters (log_L1, log_L2, log_k, log_k_mult,
#'   # richards_b, log_L0, log_rmax, log_k_gc, t50_gc).
#'   data(io_yft_data)
#'   data(io_yft_parameters)
#'   d <- io_yft_data
#'   p <- io_yft_parameters
#'   L1 <- exp(p$log_L1)
#'   L2 <- exp(p$log_L2)
#'
#'   # von Bertalanffy (Schnute), the default: single K
#'   L_vb <- get_growth(d$n_age, d$A1, d$A2, L1, L2, p$log_k)
#'
#'   # SS3 GrowthModel 3: age-specific K multipliers and plus-group adjustment
#'   L_ss3 <- get_growth(d$n_age, d$A1, d$A2, L1, L2, p$log_k,
#'                       log_k_mult = p$log_k_mult,
#'                       k_mult_ages = d$k_mult_ages,
#'                       plus_group_adjust = TRUE)
#'
#'   # Richards (b = 1 recovers VB; b near 0 approximates Gompertz)
#'   L_rich <- get_growth(d$n_age, d$A1, d$A2, L1, L2, log_k = log(0.25),
#'                        growth_model = "richards", b = -1)
#'
#'   # Growth cessation: absolute age t = model age + age_offset
#'   L_gc <- get_growth(d$n_age, d$A1, d$A2, growth_model = "cessation",
#'                      L0 = 15, log_rmax = log(7), log_k_gc = log(0.6),
#'                      t50 = 18, age_offset = 0)
#'
#'   matplot(cbind(L_vb, L_ss3, L_rich, L_gc), type = "l", lty = 1,
#'           xlab = "Model age (quarters)", ylab = "Mean length (cm)")
#'   legend("bottomright", c("vonB", "SS3 GM3", "Richards", "Cessation"),
#'          col = 1:4, lty = 1)
#' }
#' @export
get_growth <- function(n_age, A1, A2, L1 = NULL, L2 = NULL, log_k = NULL,
                       min_age = 1L,
                       log_k_mult = NULL, k_mult_ages = NULL,
                       plus_group_adjust = FALSE,
                       growth_model = "vonB",
                       b = NULL,
                       L0 = NULL, log_rmax = NULL, log_k_gc = NULL, t50 = NULL,
                       age_offset = 0) {
  .validate_growth_options(n_age, A1, min_age, log_k_mult, k_mult_ages,
                           plus_group_adjust,
                           growth_model = growth_model, A2 = A2,
                           L1 = L1, L2 = L2, log_k = log_k, b = b,
                           L0 = L0, log_rmax = log_rmax,
                           log_k_gc = log_k_gc, t50 = t50,
                           age_offset = age_offset, check_required = TRUE)
  switch(growth_model,
    vonB = .growth_vb(n_age, A1, A2, L1, L2, log_k, min_age,
                      log_k_mult, k_mult_ages, plus_group_adjust),
    richards = .growth_richards(n_age, A1, A2, L1, L2, log_k, b, min_age),
    cessation = .growth_cessation(n_age, L0, log_rmax, log_k_gc, t50,
                                  min_age, age_offset)
  )
}

# Schnute von Bertalanffy growth, optionally with SS3 GrowthModel 3
# age-specific K and the SS3 plus-group length approximation. Inputs are
# validated by get_growth().
.growth_vb <- function(n_age, A1, A2, L1, L2, log_k, min_age = 1L,
                       log_k_mult = NULL, k_mult_ages = NULL,
                       plus_group_adjust = FALSE) {
  "[<-" <- ADoverload("[<-")
  ages <- seq(min_age, by = 1, length.out = n_age)
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

# Richards (1959) growth as parameterized by Schnute (1981), following the
# Stock Synthesis user manual (Methot et al.), section 8.5.2.2.2.
.growth_richards <- function(n_age, A1, A2, L1, L2, log_k, b, min_age) {
  ages <- seq(min_age, by = 1, length.out = n_age)
  k <- exp(log_k)
  frac <- (1 - exp(-k * (ages - A1))) / (1 - exp(-k * (A2 - A1)))
  inner <- L1^b + (L2^b - L1^b) * frac
  inner^(1 / b)
}

# Growth cessation model, Maunder et al. (2018) Mar Biol 165:76, Eq. 2.
.growth_cessation <- function(n_age, L0, log_rmax, log_k_gc, t50,
                              min_age, age_offset) {
  t <- seq(min_age, by = 1, length.out = n_age) + age_offset
  rmax <- exp(log_rmax)
  k <- exp(log_k_gc)
  # Overflow-safe softplus (RTMB::logspace_add on the AD path)
  sp0 <- log1pexp(-k * t50)        # log(exp(-k t50) + 1)
  spt <- log1pexp(k * (t - t50))   # log(exp(k(t - t50)) + 1)
  L0 + rmax * ((sp0 - spt) / k + t)
}

.is_growth_model <- function(growth_model) {
  is.character(growth_model) && length(growth_model) == 1L &&
    growth_model %in% c("vonB", "richards", "cessation")
}

# Structural checks only; these must not branch on AD values. Value checks on
# initial parameters are in .validate_growth_values().
.validate_growth_options <- function(n_age, A1, min_age = 1L,
                                     log_k_mult = NULL, k_mult_ages = NULL,
                                     plus_group_adjust = FALSE,
                                     check_length = TRUE,
                                     growth_model = "vonB", A2 = NULL,
                                     L1 = NULL, L2 = NULL, log_k = NULL,
                                     b = NULL, L0 = NULL, log_rmax = NULL,
                                     log_k_gc = NULL, t50 = NULL,
                                     age_offset = 0,
                                     check_required = FALSE) {
  if (!.is_growth_model(growth_model)) {
    stop("`growth_model` must be one of \"vonB\", \"richards\" or \"cessation\".",
         call. = FALSE)
  }
  if (check_required) {
    required <- switch(growth_model,
      vonB = list(L1 = L1, L2 = L2, log_k = log_k),
      richards = list(L1 = L1, L2 = L2, log_k = log_k, b = b),
      cessation = list(L0 = L0, log_rmax = log_rmax, log_k_gc = log_k_gc, t50 = t50)
    )
    missing <- names(required)[vapply(required, is.null, logical(1))]
    if (length(missing)) {
      stop(sprintf("`growth_model = \"%s\"` requires %s.", growth_model,
                   paste0("`", missing, "`", collapse = ", ")), call. = FALSE)
    }
  }
  if (growth_model != "vonB" && (length(log_k_mult) || length(k_mult_ages))) {
    stop("`log_k_mult` and `k_mult_ages` are supported only with `growth_model = \"vonB\"`.",
         call. = FALSE)
  }
  if (!is.numeric(age_offset) || length(age_offset) != 1L || !is.finite(age_offset)) {
    stop("`age_offset` must be a finite number.", call. = FALSE)
  }
  if (growth_model == "cessation") {
    modeled <- seq(min_age, by = 1, length.out = n_age)
    for (nm in c("A1", "A2")) {
      a <- if (nm == "A1") A1 else A2
      if (length(a) != 1L || !is.finite(a) || !a %in% modeled) {
        stop(sprintf("`%s` must be a modeled age when `growth_model = \"cessation\"`.", nm),
             call. = FALSE)
      }
    }
  }
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
  if (plus_group_adjust && growth_model != "vonB") {
    stop("`plus_group_adjust` is supported only with `growth_model = \"vonB\"`.", call. = FALSE)
  }
  if (plus_group_adjust && (max_age <= 0 || max_age != floor(max_age))) {
    stop("Plus-group adjustment requires a positive integer oldest age.", call. = FALSE)
  }
  invisible(NULL)
}

# Value checks on initial (non-AD) parameter values, run by opal_obj().
# opal_model() evaluates growth at ages 1, ..., n_age (the get_growth()
# default), so A1 is compared with 1.
.validate_growth_values <- function(data, parameters = NULL) {
  growth_model <- if (is.null(data$growth_model)) "vonB" else data$growth_model
  if (!.is_growth_model(growth_model)) {
    return(invisible(NULL))  # reported by .validate_growth_options()
  }
  if (!is.null(parameters)) {
    required <- switch(growth_model,
      vonB = character(), richards = "richards_b",
      cessation = c("log_L0", "log_rmax", "log_k_gc", "t50_gc"))
    missing <- setdiff(required, names(parameters))
    if (length(missing)) {
      stop(sprintf("`growth_model = \"%s\"` requires parameter(s) %s; see get_parameters().",
                   growth_model, paste0("`", missing, "`", collapse = ", ")),
           call. = FALSE)
    }
  }
  if (growth_model == "richards") {
    b <- parameters$richards_b
    if (!is.null(b) && (length(b) != 1L || !is.finite(b) || b == 0)) {
      stop("`richards_b` must be a finite, non-zero number.", call. = FALSE)
    }
    if (length(data$A1) == 1L && is.finite(data$A1) && data$A1 > 1) {
      warning("Richards growth with `A1` greater than the youngest modeled age: ",
              "some parameter combinations give NaN predicted lengths at young ages.",
              call. = FALSE)
    }
  }
  if (growth_model == "cessation" && is.null(data$growth_age_offset)) {
    message("Growth cessation uses absolute age t = model age + ",
            "`growth_age_offset`; `growth_age_offset` is not supplied and defaults to 0.")
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
