#' Logistic selectivity as a function of length
#'
#' @param len Numeric vector of length-bin midpoints, or of ages for an
#'   age-based fishery (`sel_basis_f = 2`); see [get_selectivity()].
#' @param par Numeric vector of length 6 containing selectivity parameters.
#'   Only the first two entries are used: `par[1]` (a) and `par[2]` (b).
#'   `a` is the inflection point on the real line, transformed via
#'   `mean(len) + a * sd(len)`, so `a = 0` gives inflection at `mean(len)`.
#'   `b` is log-scale 95% width, transformed via `exp(b) * sd(len)`,
#'   so `b = 0` gives a 95% width equal to `sd(len)`.
#' @return Numeric vector of selectivity values in (0, 1].
#' @export
#'
sel_logistic <- function(len, par) {
  "[<-" <- ADoverload("[<-")
  mu <- mean(len)
  sd <- sd(len)
  a <- par[1]
  b <- par[2]
  a50 <- mu + a * sd
  w95 <- exp(b) * sd
  sel <- 1 / (1 + exp(-log(19) * (len - a50) / w95))
  return(sel)
}

#' log(1 + exp(u)), overflow-safe
#'
#' The AD path uses RTMB's vectorized atomic `logspace_add`; the numeric path
#' uses the standard stable form.
#'
#' @importFrom RTMB logspace_add
#' @noRd
log1pexp <- function(u) {
  if (inherits(u, "advector")) {
    logspace_add(0 * u, u)
  } else {
    pmax(u, 0) + log1p(exp(-abs(u)))
  }
}

#' log(2^g - 1), overflow-safe for positive g
#'
#' @noRd
log_2pow_m1 <- function(g) {
  g * log(2) + log(1 - exp(-g * log(2)))
}

#' Double Richards selectivity as a function of length
#'
#' A smooth, six-parameter curve formed by multiplying an ascending
#' generalized-logistic (Richards) limb by one minus a descending Richards
#' limb:
#' \deqn{(1 + \exp(-z_1))^{-1/\gamma_1}
#' [1 - (1 + \exp(-z_2))^{-1/\gamma_2}].}
#' Parameters are on the real line and transformed as follows:
#' \describe{
#'   \item{`par[1]`}{Ascending 50% point: `mean(len) + par[1] * sd(len)`.}
#'   \item{`par[2]`}{Ascending width: `exp(par[2]) * sd(len)`.}
#'   \item{`par[3]`}{Ascending shape: `gamma1 = exp(par[3])`.}
#'   \item{`par[4]`}{Descending 50% point is ascending 50% point plus
#'     `exp(par[4]) * sd(len)`.}
#'   \item{`par[5]`}{Descending width: `exp(par[5]) * sd(len)`.}
#'   \item{`par[6]`}{Descending shape: `gamma2 = exp(par[6])`.}
#' }
#' The Maunder parameterization is recovered with
#' `alpha_k = L50_k + log(2^gamma_k - 1) / beta_k` and
#' `beta_k = log(19) / width_k`. The peak can be below one when the limbs
#' overlap; it is intentionally not rescaled. As either shape parameter tends
#' to negative infinity, its limb tends to a Gompertz curve with the same 50%
#' point. Large shape values flatten the limb toward 0.5 (at `par[3] =
#' par[6] = 7`, the curve is approximately 0.25). Mapping `par[3]` and
#' `par[6]` to zero gives a four-parameter double-logistic curve. Suggested
#' bounds for columns 3 and 6 are `c(-5, 5)`; use the existing `par_sel` bounds
#' for the remaining columns.
#'
#' @param len Numeric vector of length-bin midpoints, or of ages for an
#'   age-based fishery (`sel_basis_f = 2`); see [get_selectivity()].
#' @param par Numeric vector of length 6 containing selectivity parameters.
#' @return Numeric vector of selectivity values in [0, 1).
#' @references Maunder, M. (2025). Double Richards selectivity. Unpublished
#'   technical note, 19 November 2025.
#' @export
sel_double_richards <- function(len, par) {
  mu <- mean(len)
  sd <- sd(len)
  l50_1 <- mu + par[1] * sd
  w_1 <- exp(par[2]) * sd
  g_1 <- exp(par[3])
  l50_2 <- l50_1 + exp(par[4]) * sd
  w_2 <- exp(par[5]) * sd
  g_2 <- exp(par[6])
  b_1 <- log(19) / w_1
  b_2 <- log(19) / w_2
  z_1 <- b_1 * (len - l50_1) - log_2pow_m1(g_1)
  z_2 <- b_2 * (len - l50_2) - log_2pow_m1(g_2)
  asc <- exp(-log1pexp(-z_1) / g_1)
  dsc <- 1 - exp(-log1pexp(-z_2) / g_2)
  asc * dsc
}

#' Convert double Richards parameters to natural scale
#'
#' This numeric reporting helper returns the transformed parameter values and
#' Maunder's `alpha1`, `beta1`, `alpha2`, and `beta2`. Under the 50%-point
#' parameterization, `alpha2` can be non-positive when the shape parameters
#' differ.
#'
#' @param len Numeric vector of length-bin midpoints.
#' @param par Numeric vector of length 6 containing selectivity parameters.
#' @return Named numeric vector of natural-scale parameters.
#' @export
double_richards_natural <- function(len, par) {
  mu <- mean(len)
  sd <- sd(len)
  l50_asc <- mu + par[1] * sd
  width_asc <- exp(par[2]) * sd
  gamma_asc <- exp(par[3])
  l50_desc <- l50_asc + exp(par[4]) * sd
  width_desc <- exp(par[5]) * sd
  gamma_desc <- exp(par[6])
  beta1 <- log(19) / width_asc
  beta2 <- log(19) / width_desc
  alpha1 <- l50_asc + log_2pow_m1(gamma_asc) / beta1
  alpha_desc <- l50_desc + log_2pow_m1(gamma_desc) / beta2
  c(l50_asc = l50_asc, width_asc = width_asc, gamma_asc = gamma_asc,
    l50_desc = l50_desc, width_desc = width_desc, gamma_desc = gamma_desc,
    alpha1 = alpha1, beta1 = beta1, alpha2 = alpha_desc - alpha1,
    beta2 = beta2)
}

#' Selectivity-at-length for a single fishery, dispatched on type code
#'
#' @param len Numeric vector of length-bin midpoints, or of ages for an
#'   age-based fishery (`sel_basis_f = 2`); see [get_selectivity()].
#' @param par Numeric vector of length 6 containing selectivity parameters.
#' @param sel_type Selectivity type: 1 = logistic, 2 = double-normal, 3 =
#'   double Richards.
#' @return Numeric vector of selectivity values.
#' @export
sel_length <- function(len, par, sel_type) {
  switch(as.character(sel_type),
    "1" = sel_logistic(len, par),
    "2" = sel_double_normal(len, par),
    "3" = sel_double_richards(len, par),
    stop("Unknown sel_type: ", sel_type,
         ". Supported: 1 (logistic), 2 (double-normal), 3 (double Richards).")
  )
}

#' Double-normal selectivity as a function of length (SS3 pattern 24, full form)
#'
#' All parameters are on the real line, making this parameterization suitable
#' for unconstrained optimization and AD. The length-bin vector `x` is used
#' to define natural scales for position and width parameters.
#'
#' Only the full form is implemented (e and f both active). The simplified
#' forms where e <= -999 or f <= -999 are not supported.
#'
#' @param x Numeric vector of length-bin midpoints, or of ages for an
#'   age-based fishery (`sel_basis_f = 2`); see [get_selectivity()]. The
#'   values must be equally spaced, because `x[2] - x[1]` is used as the bin
#'   width.
#' @param par Numeric vector of length 6 containing selectivity parameters.
#'   \describe{
#'     \item{`par[1]` (a)}{Peak location (real line). Transformed via
#'       `mean(x) + a * sd(x)`, so `a = 0` places the peak at `mean(x)`.}
#'     \item{`par[2]` (b)}{Plateau width (real line). Controls the distance from peak to the
#'       start of the descending limb via logistic transform of the available
#'       range: `peak + bin_width + (0.99 * max(x) - peak - bin_width) / (1 + exp(-b))`.}
#'     \item{`par[3]` (c)}{Ascending width (real line, log-space). Actual denominator =
#'       `exp(c) * sd(x)^2`.}
#'     \item{`par[4]` (d)}{Descending width (real line, log-space). Actual denominator =
#'       `exp(d) * sd(x)^2`.}
#'     \item{`par[5]` (e)}{Initial selectivity (real line, logit-space). Transformed via
#'       `1 / (1 + exp(-e))`, so `e = 0` gives initial selectivity of 0.5,
#'       large negative values give ~0, large positive values give ~1.}
#'     \item{`par[6]` (f)}{Final selectivity (real line, logit-space). Same transform as `e`.}
#'   }
#' @return Numeric vector of selectivity values in \[0, 1\].
#' @export
#'
sel_double_normal <- function(x, par) {
  "[<-" <- ADoverload("[<-")
  mu <- mean(x)
  sd <- sd(x)
  # Extract parameters from vector
  a <- par[1]
  b <- par[2]
  c <- par[3]
  d <- par[4]
  e <- par[5]
  f <- par[6]
  # --- Transform parameters from real line to natural scale ---
  peak       <- mu + a * sd
  upselex    <- exp(c) * sd^2
  downselex  <- exp(d) * sd^2
  point1     <- 1 / (1 + exp(-e))
  point2     <- 1 / (1 + exp(-f))

  # --- Derived quantities ---
  startbin  <- 1
  j1        <- startbin - 1
  j2        <- length(x)
  bin_width <- x[2] - x[1]

  peak2  <- peak + bin_width + (0.99 * x[j2] - peak - bin_width) / (1 + exp(-b))

  t1min  <- exp(-(x[startbin] - peak)^2 / upselex)
  t2min  <- exp(-(x[j2] - peak2)^2 / downselex)

  # --- Compute selectivity ---
  t1 <- x - peak
  t2 <- x - peak2
  join1 <- 1 / (1 + exp(-(20 / (1 + abs(t1))) * t1))
  join2 <- 1 / (1 + exp(-(20 / (1 + abs(t2))) * t2))
  asc <- point1 + (1 - point1) * (exp(-t1^2 / upselex) - t1min) / (1 - t1min)
  dsc <- 1 + (point2 - 1) * (exp(-t2^2 / downselex) - 1) / (t2min - 1)
  sel <- rep(NA_real_, length(x))
  idx <- (j1 + 1):j2
  sel[idx] <- asc[idx] * (1 - join1[idx]) + join1[idx] * (1 - join2[idx] + dsc[idx] * join2[idx])
  if (startbin > 1) {
    sel[1:startbin] <- (x[1:startbin] / x[startbin])^2 * sel[startbin]
  }
  if (j2 < length(x)) {
    sel[(j2 + 1):length(x)] <- sel[j2]
  }
  return(sel)
}

#' Probability of length at age matrix (age-length key)
#'
#' Computes the probability distribution of fish lengths across different
#' age classes. AD-compatible: uses RTMB::pnorm and ADoverload so that
#' gradients propagate if growth parameters (mu_a, sd_a) are estimated.
#'
#' Bins are assumed contiguous and non-overlapping, i.e.
#' \code{len_upper[i] == len_lower[i + 1]} for all \code{i}. This allows a
#' single vectorized \code{pnorm} call per age on the \code{n_len + 1} unique
#' bin edges rather than two calls on the lower and upper bounds separately.
#'
#' @param len_lower Numeric vector of lower bounds of length bins (length L). Data only.
#' @param len_upper Numeric vector of upper bounds of length bins (length L). Data only.
#' @param mu_a Numeric vector of mean length at age (length A). May be AD.
#' @param sd_a Numeric vector of SD of length at age (length A). May be AD.
#' @return Matrix of dimensions L x A where columns sum to 1.
#' @importFrom RTMB ADoverload pnorm
#' @export
#'
get_pla <- function(len_lower, len_upper, mu_a, sd_a) {
  "[<-" <- ADoverload("[<-")
  n_len <- length(len_lower)
  n_age <- length(mu_a)
  edges <- c(len_lower, len_upper[n_len])
  pla <- array(0, dim = c(n_len, n_age))
  for (a in seq_len(n_age)) {
    p_edges <- pnorm((edges - mu_a[a]) / sd_a[a])
    pla[, a] <- p_edges[2:(n_len + 1)] - p_edges[1:n_len]
    col_sum <- sum(pla[, a])
    pla[, a] <- pla[, a] / (col_sum + 1e-12)
  }
  return(pla)
}

#' Validate the selectivity basis codes
#'
#' @param sel_basis_f Integer vector of selectivity basis codes, or `NULL`.
#' @param n_fishery Number of fisheries.
#' @return Integer vector of length `n_fishery` (all `1L` when `sel_basis_f`
#'   is `NULL`).
#' @noRd
.validate_sel_basis <- function(sel_basis_f, n_fishery) {
  if (is.null(sel_basis_f)) return(rep(1L, n_fishery))
  if (length(sel_basis_f) != n_fishery) {
    stop("'sel_basis_f' must have length n_fishery (", n_fishery,
         "), not ", length(sel_basis_f), ".", call. = FALSE)
  }
  if (!is.numeric(sel_basis_f) || anyNA(sel_basis_f) ||
      !all(sel_basis_f %in% c(1, 2))) {
    bad <- unique(sel_basis_f[is.na(sel_basis_f) | !sel_basis_f %in% c(1, 2)])
    stop("'sel_basis_f' must contain only 1 (length) or 2 (age); found: ",
         paste(bad, collapse = ", "), ".", call. = FALSE)
  }
  as.integer(sel_basis_f)
}

#' Compute selectivity-at-age from parametric selectivity curves
#'
#' Defines selectivity for each fishery as a parametric curve (logistic,
#' double-normal, or double Richards) evaluated either on length or directly
#' on age, chosen per fishery by `data$sel_basis_f`:
#' \describe{
#'   \item{1 = length (default)}{The curve is evaluated at `len_mid` and
#'     converted to selectivity-at-age by matrix-multiplying with the
#'     probability-of-length-at-age (`t(pla) %*% sel_at_length`).}
#'   \item{2 = age}{The curve is evaluated directly at the ages `age_a`, and
#'     `pla` is not used, so the selectivity does not depend on growth.
#'     Parameters are scaled by `mean(age_a)` and `sd(age_a)` rather than by
#'     the length bins.}
#' }
#' If `data$sel_basis_f` is absent or `NULL`, every fishery is length-based.
#'
#' Selectivity is currently time-invariant within each fishery (constant
#' across years).
#'
#' @param data A list containing model data. Required elements:
#'   n_fishery, n_year, n_age, sel_type_f. Optional: `sel_basis_f`, an integer
#'   vector of length n_fishery (1 = length, 2 = age), and `min_age` (default
#'   1), used to build `age_a` when it is not supplied.
#' @param par_sel Numeric matrix of dimensions `[n_fishery, 6]`. Each row is
#'   a real-line parameter vector. For logistic (sel_type_f == 1), only
#'   columns 1:2 are used. For double-normal (sel_type_f == 2) and double
#'   Richards (sel_type_f == 3), all 6 are used.
#' @param pla Numeric matrix of dimensions `[n_len, n_age]` containing
#'   probability-of-length-at-age (age-length key). Computed via
#'   \code{get_pla()}. Used only for length-based fisheries.
#' @param len_mid Numeric vector (length n_len) of length-bin midpoints.
#' @param age_a Numeric vector (length n_age) of ages at which age-based
#'   curves are evaluated. Defaults to
#'   `seq(min_age, by = 1, length.out = n_age)`, as in \code{get_growth()}.
#' @return 3D array `sel_fya` of dimensions `[n_fishery, n_year, n_age]`.
#' @importFrom RTMB ADoverload
#' @export
#'
get_selectivity <- function(data, par_sel, pla, len_mid, age_a = NULL) {
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  n_fishery <- data$n_fishery
  n_year <- data$n_year
  n_age <- data$n_age
  sel_basis_f <- .validate_sel_basis(data$sel_basis_f, n_fishery)
  if (is.null(age_a)) {
    min_age <- if (is.null(data$min_age)) 1L else data$min_age
    age_a <- seq(min_age, by = 1, length.out = n_age)
  } else if (length(age_a) != n_age) {
    stop("'age_a' must have length n_age (", n_age, ").", call. = FALSE)
  }
  sel_fya <- array(0, dim = c(n_fishery, n_year, n_age))
  for (f in seq_len(n_fishery)) {
    par_f <- par_sel[f, ]
    if (sel_basis_f[f] == 2L) {
      # Age-based: evaluate the curve directly on age
      sel_at_age <- sel_length(age_a, par_f, data$sel_type_f[f])
    } else {
      sel_at_length <- sel_length(len_mid, par_f, data$sel_type_f[f])
      # Convert to selectivity-at-age via PLA
      sel_at_age <- c(t(pla) %*% sel_at_length)
    }
    # Time-invariant: replicate across years
    for (y in seq_len(n_year)) {
      sel_fya[f, y, ] <- sel_at_age
    }
  }
  return(sel_fya)
}

#' Selectivity scaling vector for each fishery
#'
#' @param basis Character vector of length 1 or `n_fishery`, each `"length"`
#'   or `"age"`.
#' @param n_fishery Number of fisheries.
#' @param sel_lengths Length-bin midpoints.
#' @param sel_ages Ages, required when any fishery is age-based.
#' @return List of length `n_fishery` holding the x vector for each fishery.
#' @noRd
.selex_scale_vectors <- function(basis, n_fishery, sel_lengths, sel_ages) {
  if (!length(basis) %in% c(1L, n_fishery)) {
    stop("'basis' must have length 1 or n_fishery (", n_fishery, ").",
         call. = FALSE)
  }
  bad <- setdiff(basis, c("length", "age"))
  if (length(bad) > 0L) {
    stop("'basis' must be \"length\" or \"age\"; found: ",
         paste(bad, collapse = ", "), ".", call. = FALSE)
  }
  basis <- rep_len(basis, n_fishery)
  if (any(basis == "age") && is.null(sel_ages)) {
    stop("'sel_ages' must be supplied for age-based fisheries.", call. = FALSE)
  }
  if (any(basis == "length") && is.null(sel_lengths)) {
    stop("'sel_lengths' must be supplied for length-based fisheries.",
         call. = FALSE)
  }
  lapply(basis, function(b) if (b == "age") sel_ages else sel_lengths)
}

#' Convert SS3 selectivity parameters to RTMB real-line parameterization
#'
#' Takes SS3 natural-scale selectivity parameters and converts them to the
#' centered real-line parameterization used by the RTMB \code{sel_logistic}
#' and \code{sel_double_normal} functions.
#'
#' Length-based fisheries (`basis = "length"`) take SS3 `size_selex`
#' patterns 1 (logistic) and 24 (double normal), scaled by `sel_lengths`.
#' Age-based fisheries (`basis = "age"`, used with `sel_basis_f = 2`) take
#' SS3 `age_selex` patterns 12 (logistic) and 20 (double normal), which share
#' the parameter definitions of patterns 1 and 24, scaled by `sel_ages`.
#'
#' SS3's `-999` codes for the initial and final double-normal selectivity
#' ("no adjustment") have no exact equivalent in opal. They are converted to
#' `-9`, which sets the initial or final selectivity close to zero. With a
#' final value of `-999`, SS3 leaves the descending limb unadjusted at the
#' oldest age or length bin, so the two curves can differ there.
#'
#' @param ss3_pars A data.frame or matrix with one row per fishery. For
#'   double-normal (pattern 24 or 20): columns are peak, top_logit,
#'   ascend_se, descend_se, start_logit, end_logit. For logistic (pattern 1
#'   or 12): columns are inflection, width.
#' @param sel_type_f Integer vector (length n_fishery). 1 = logistic, 2 =
#'   double-normal, 3 = double Richards.
#' @param sel_lengths Numeric vector of selectivity length-bin midpoints
#'   (same vector that will be passed to sel_logistic/sel_double_normal).
#'   Used for length-based fisheries.
#' @param basis Character vector of length 1 or n_fishery: `"length"`
#'   (default) or `"age"`, giving the selectivity basis of each fishery.
#' @param sel_ages Numeric vector of ages at which age-based curves are
#'   evaluated, i.e. `seq(min_age, by = 1, length.out = n_age)`. Required
#'   when any fishery has `basis = "age"`.
#' @return Numeric matrix `[n_fishery, 6]` of RTMB real-line parameters.
#' @export
#'
convert_ss3_selex_to_rtmb <- function(ss3_pars, sel_type_f, sel_lengths,
                                      basis = "length", sel_ages = NULL) {
  n_fishery <- length(sel_type_f)
  x_f <- .selex_scale_vectors(basis, n_fishery, sel_lengths, sel_ages)
  par_sel <- matrix(0, nrow = n_fishery, ncol = 6)

  for (f in seq_len(n_fishery)) {
    mu_len <- mean(x_f[[f]])
    sd_len <- sd(x_f[[f]])
    if (sel_type_f[f] == 1L) {
      # Logistic (SS3 pattern 1 or 12): inflection, width
      ss3_inflection <- ss3_pars[f, 1]
      ss3_width <- ss3_pars[f, 2]
      par_sel[f, 1] <- (ss3_inflection - mu_len) / sd_len
      par_sel[f, 2] <- log(ss3_width / sd_len)
      # par_sel[f, 3:6] remain 0 (unused)
    } else if (sel_type_f[f] == 2L) {
      # Double-normal (SS3 pattern 24 or 20): peak, top_logit, ascend_se, descend_se, start_logit, end_logit
      ss3_peak        <- ss3_pars[f, 1]
      ss3_top_logit   <- ss3_pars[f, 2]
      ss3_ascend_se   <- ss3_pars[f, 3]
      ss3_descend_se  <- ss3_pars[f, 4]
      ss3_start_logit <- ss3_pars[f, 5]
      ss3_end_logit   <- ss3_pars[f, 6]

      par_sel[f, 1] <- (ss3_peak - mu_len) / sd_len        # a: peak location
      par_sel[f, 2] <- ss3_top_logit                         # b: plateau (same space)
      par_sel[f, 3] <- ss3_ascend_se - 2 * log(sd_len)           # c: ascending width
      par_sel[f, 4] <- ss3_descend_se - 2 * log(sd_len)          # d: descending width

      # Handle start_logit = -999 (SS3 convention for "fix at 0")
      if (ss3_start_logit <= -999) {
        par_sel[f, 5] <- -9
      } else {
        par_sel[f, 5] <- ss3_start_logit                     # e: initial selectivity
      }

      # Handle end_logit = -999
      if (ss3_end_logit <= -999) {
        par_sel[f, 6] <- -9
      } else {
        par_sel[f, 6] <- ss3_end_logit                       # f: final selectivity
      }
    } else if (sel_type_f[f] == 3L) {
      stop("sel_type_f = 3 (double Richards) has no SS3 equivalent; supply par_sel directly for these fisheries.")
    } else {
      stop("Unknown sel_type_f: ", sel_type_f[f])
    }
  }
  return(par_sel)
}

#' Convert RTMB real-line selectivity parameters back to SS3 natural scale
#'
#' Inverse of \code{convert_ss3_selex_to_rtmb}. Useful for verifying
#' round-trip conversion and reporting parameter values in natural units.
#'
#' @param par_sel Numeric matrix `[n_fishery, 6]` of RTMB real-line parameters.
#' @param sel_type_f Integer vector (1 = logistic, 2 = double-normal, 3 =
#'   double Richards).
#' @param sel_lengths Numeric vector of selectivity length-bin midpoints.
#'   Used for length-based fisheries.
#' @inheritParams convert_ss3_selex_to_rtmb
#' @return Data.frame with SS3-scale parameter values per fishery. Double
#'   Richards rows retain `NA` parameter values because there is no SS3
#'   equivalent. Age-based rows correspond to SS3 `age_selex` patterns 12
#'   and 20.
#' @export
#'
convert_rtmb_selex_to_ss3 <- function(par_sel, sel_type_f, sel_lengths,
                                      basis = "length", sel_ages = NULL) {
  n_fishery <- nrow(par_sel)
  x_f <- .selex_scale_vectors(basis, n_fishery, sel_lengths, sel_ages)
  result <- data.frame(
    fleet = 1:n_fishery,
    type = sel_type_f,
    peak_or_inflection = NA_real_,
    top_logit_or_width = NA_real_,
    ascend_se = NA_real_,
    descend_se = NA_real_,
    start_logit = NA_real_,
    end_logit = NA_real_
  )
  for (f in seq_len(n_fishery)) {
    mu_len <- mean(x_f[[f]])
    sd_len <- sd(x_f[[f]])
    if (sel_type_f[f] == 1L) {
      result$peak_or_inflection[f] <- mu_len + par_sel[f, 1] * sd_len
      result$top_logit_or_width[f] <- exp(par_sel[f, 2]) * sd_len
    } else if (sel_type_f[f] == 2L) {
      result$peak_or_inflection[f] <- mu_len + par_sel[f, 1] * sd_len
      result$top_logit_or_width[f] <- par_sel[f, 2]
      result$ascend_se[f]  <- par_sel[f, 3] + 2 * log(sd_len)
      result$descend_se[f] <- par_sel[f, 4] + 2 * log(sd_len)
      result$start_logit[f] <- par_sel[f, 5]
      result$end_logit[f] <- par_sel[f, 6]
    } else if (sel_type_f[f] != 3L) {
      stop("Unknown sel_type_f: ", sel_type_f[f])
    }
  }
  return(result)
}
