#' Get priors
#' 
#' @param parameters A \code{list} specifying the parameters to be passed to \code{MakeADFun}. Can be generated using the `get_parameters()` function. Vector parameters (e.g., \code{log_cpue_q}) are supported.
#' @param data A \code{list} of data inputs (optional). Used to retrieve prior
#'   center values for growth/variability parameters (e.g.,
#'   \code{data$prior_log_L1_mean}).
#' @return A \code{list} of priors.
#' @export
#' 
get_priors <- function(parameters, data = NULL) {
  priors <- list()
  log_B0_mean <- if (!is.null(data) && !is.null(data$prior_log_B0_mean)) data$prior_log_B0_mean else parameters$log_B0
  priors[["log_B0"]] <- list(type = "normal", par1 = log_B0_mean, par2 = 1.0, index = which("log_B0" == names(parameters)))
  log_cpue_q_mean <- if (!is.null(data) && !is.null(data$prior_log_cpue_q_mean)) data$prior_log_cpue_q_mean else parameters$log_cpue_q
  priors[["log_cpue_q"]] <- list(type = "normal", par1 = log_cpue_q_mean, par2 = 1.5, index = which("log_cpue_q" == names(parameters)))
  
  # Broad normal priors for selectivity parameters (all on real line)
  # par_sel is a matrix [n_fishery, 6] — treat as a single block
  # Normal(0, 2) is vague: allows peak to shift ~2 SD from mean length,
  # widths to vary by exp(±2) ≈ 0.14x to 7.4x the length SD, etc.
  # Imported SS3-style matrices can use large negative sentinel values for
  # inactive parameters. Center those entries on their fixed value so they do
  # not add a misleading constant to the objective.
  if ("par_sel" %in% names(parameters)) {
    par_sel_mean <- ifelse(parameters$par_sel <= -100, parameters$par_sel, 0)
    priors[["par_sel"]] <- list(type = "normal", par1 = par_sel_mean, par2 = 2, index = which("par_sel" == names(parameters)))
  }

  if ("log_lf_tau" %in% names(parameters)) {
    priors[["log_lf_tau"]] <- list(type = "normal", par1 = 0, par2 = 2, index = which("log_lf_tau" == names(parameters)))
  }

  # Growth priors (on log scale — Normal priors on log-parameters are lognormal on natural scale)
  # Prior centers are stored in data to match initialization values
  if ("log_L1" %in% names(parameters)) {
    log_L1_mean <- if (!is.null(data) && !is.null(data$prior_log_L1_mean)) data$prior_log_L1_mean else parameters$log_L1
    priors[["log_L1"]] <- list(type = "normal", par1 = log_L1_mean, par2 = 0.1,
                                index = which("log_L1" == names(parameters)))
  }
  if ("log_L2" %in% names(parameters)) {
    log_L2_mean <- if (!is.null(data) && !is.null(data$prior_log_L2_mean)) data$prior_log_L2_mean else parameters$log_L2
    priors[["log_L2"]] <- list(type = "normal", par1 = log_L2_mean, par2 = 0.2,
                                index = which("log_L2" == names(parameters)))
  }
  if ("log_k" %in% names(parameters)) {
    log_k_mean <- if (!is.null(data) && !is.null(data$prior_log_k_mean)) data$prior_log_k_mean else parameters$log_k
    priors[["log_k"]] <- list(type = "normal", par1 = log_k_mean, par2 = 0.3,
                               index = which("log_k" == names(parameters)))
  }
  if ("log_k_mult" %in% names(parameters) && length(parameters$log_k_mult)) {
    log_k_mult_mean <- if (!is.null(data) && !is.null(data$prior_log_k_mult_mean)) data$prior_log_k_mult_mean else parameters$log_k_mult
    priors[["log_k_mult"]] <- list(type = "normal", par1 = log_k_mult_mean, par2 = 0.3,
                                 index = which("log_k_mult" == names(parameters)))
  }
  # Alternative growth curves: richards_b and t50_gc are on the natural scale
  if ("richards_b" %in% names(parameters)) {
    richards_b_mean <- if (!is.null(data) && !is.null(data$prior_richards_b_mean)) data$prior_richards_b_mean else parameters$richards_b
    priors[["richards_b"]] <- list(type = "normal", par1 = richards_b_mean, par2 = 0.5,
                                   index = which("richards_b" == names(parameters)))
  }
  if ("log_L0" %in% names(parameters)) {
    log_L0_mean <- if (!is.null(data) && !is.null(data$prior_log_L0_mean)) data$prior_log_L0_mean else parameters$log_L0
    priors[["log_L0"]] <- list(type = "normal", par1 = log_L0_mean, par2 = 0.1,
                               index = which("log_L0" == names(parameters)))
  }
  if ("log_rmax" %in% names(parameters)) {
    log_rmax_mean <- if (!is.null(data) && !is.null(data$prior_log_rmax_mean)) data$prior_log_rmax_mean else parameters$log_rmax
    priors[["log_rmax"]] <- list(type = "normal", par1 = log_rmax_mean, par2 = 0.3,
                                 index = which("log_rmax" == names(parameters)))
  }
  if ("log_k_gc" %in% names(parameters)) {
    log_k_gc_mean <- if (!is.null(data) && !is.null(data$prior_log_k_gc_mean)) data$prior_log_k_gc_mean else parameters$log_k_gc
    priors[["log_k_gc"]] <- list(type = "normal", par1 = log_k_gc_mean, par2 = 0.5,
                                 index = which("log_k_gc" == names(parameters)))
  }
  if ("t50_gc" %in% names(parameters)) {
    t50_gc_mean <- if (!is.null(data) && !is.null(data$prior_t50_gc_mean)) data$prior_t50_gc_mean else parameters$t50_gc
    priors[["t50_gc"]] <- list(type = "normal", par1 = t50_gc_mean, par2 = 2,
                               index = which("t50_gc" == names(parameters)))
  }

  # Variability priors
  if ("log_CV1" %in% names(parameters)) {
    log_CV1_mean <- if (!is.null(data) && !is.null(data$prior_log_CV1_mean)) data$prior_log_CV1_mean else parameters$log_CV1
    priors[["log_CV1"]] <- list(type = "normal", par1 = log_CV1_mean, par2 = 0.3,
                                 index = which("log_CV1" == names(parameters)))
  }
  if ("log_CV2" %in% names(parameters)) {
    log_CV2_mean <- if (!is.null(data) && !is.null(data$prior_log_CV2_mean)) data$prior_log_CV2_mean else parameters$log_CV2
    priors[["log_CV2"]] <- list(type = "normal", par1 = log_CV2_mean, par2 = 0.3,
                                 index = which("log_CV2" == names(parameters)))
  }

  evaluate_priors(parameters, priors)
  return(priors)
}

#' Evaluate priors
#' 
#' @param parameters A \code{list} specifying the parameters to be passed to \code{MakeADFun}. Can be generated using the `get_parameters()` function.
#' @param priors A \code{list} of named \code{list}s specifying priors for the parameters. Can be generated using the `get_priors()` function.
#' @return A \code{numeric} value.
#' @details Supported distributions are `normal`, `student` (three degrees of
#'   freedom), `lognormal`, and `beta` (mean and precision). Priors act on the
#'   stored parameter scale; for example, a normal prior on `log_B0` is a normal
#'   density on log spawning output. Locations and scales must be scalars or
#'   match the parameter block length. An optional `lambda` field applies a
#'   non-negative scalar or per-element multiplier to the log density; it
#'   defaults to 1. Indices identify blocks in `parameters`.
#' @importFrom RTMB dnorm dlnorm
#' @importFrom RTMBdist dbeta2 dt2
#' @export
#' @examples
#' {
#'   parameters <- list(log_B0 = log(1e6))
#'   priors <- list(
#'     log_B0 = list(
#'       type = "normal", par1 = log(1e6), par2 = 1,
#'       index = which("log_B0" == names(parameters))
#'     )
#'   )
#'   evaluate_priors(parameters, priors)
#' }
#' 
evaluate_priors <- function(parameters, priors) {
  n <- length(priors)
  if (n == 0L) return(0)
  .opal_validate_priors(parameters, priors)
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  lp <- numeric(n)
  for (i in 1:n) {
    type <- priors[[i]]$type
    par1 <- priors[[i]]$par1
    par2 <- priors[[i]]$par2
    lambda <- priors[[i]]$lambda
    idx <- priors[[i]]$index
    x <- parameters[[idx]]
    if (is.null(lambda)) lambda <- 1
    if (type == "normal") lp[i] <- sum(.opal_weight_likelihood(
      dnorm(x = x, mean = par1, sd = par2, log = TRUE), lambda))
    if (type == "student") lp[i] <- sum(.opal_weight_likelihood(
      dt2(x = x, mu = par1, sigma = par2, df = 3, log = TRUE), lambda))
    # if (type == "cauchy") lp[i] <- sum(dcauchy(x = x, location = par1, scale = par2, log = TRUE))
    if (type == "lognormal") lp[i] <- sum(.opal_weight_likelihood(
      dlnorm(x = x, meanlog = par1, sdlog = par2, log = TRUE), lambda))
    if (type == "beta") lp[i] <- sum(.opal_weight_likelihood(
      dbeta2(x = x, mu = par1, phi = par2, log = TRUE), lambda))
  }
  return(-sum(lp))
}

.opal_validate_priors <- function(parameters, priors) {
  for (i in seq_along(priors)) {
    p <- priors[[i]]
    if (!is.list(p) || length(p$type) != 1L || is.na(p$type) ||
        !p$type %in% c("normal", "student", "lognormal", "beta")) {
      stop("Unsupported prior distribution in prior ", i, ".", call. = FALSE)
    }
    idx <- p$index
    if (!is.numeric(idx) || length(idx) != 1L || !is.finite(idx) ||
        idx != floor(idx) || idx < 1L || idx > length(parameters)) {
      stop("Invalid parameter index in prior ", i, ".", call. = FALSE)
    }
    label <- names(priors)[i]
    if (length(label) && label %in% names(parameters) && label != names(parameters)[idx]) {
      stop("Prior name and parameter index disagree in prior ", i, ".", call. = FALSE)
    }
    for (field in c("par1", "par2", "lambda")) {
      if (field == "lambda" && is.null(p$lambda)) next
      value <- p[[field]]
      if (!is.numeric(value) || !length(value) || any(!is.finite(value)) ||
          !length(value) %in% c(1L, length(parameters[[idx]]))) {
        stop("Prior ", field, " must be finite with length one or the parameter length.",
             call. = FALSE)
      }
      if (field == "lambda" && any(value < 0)) {
        stop("Prior lambda must be non-negative.", call. = FALSE)
      }
    }
    if (any(p$par2 <= 0)) stop("Prior scale or precision must be positive.", call. = FALSE)
    if (p$type == "beta" && any(p$par1 <= 0 | p$par1 >= 1)) {
      stop("Beta prior mean must be between zero and one.", call. = FALSE)
    }
  }
  invisible(priors)
}
