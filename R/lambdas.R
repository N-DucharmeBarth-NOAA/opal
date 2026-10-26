#' Resolve likelihood-component lambdas
#'
#' Expands the optional `data$lambdas` list to one multiplier per CPUE index or
#' composition fleet. Scalar values are recycled; unspecified components use 1.
#' @details The supported entries are `cpue` (one or `n_index` values), `lf`
#'   and `wf` (one or `n_fishery` values each), and scalar `rec`, `init_rec`,
#'   `prior`, and `penalty` multipliers. `prior_mapped_only` is a logical flag
#'   that removes prior contributions for elements fixed by the parameter map.
#' @param data A model data list.
#' @return A fully expanded named list of lambda values.
#' @export
get_lambdas <- function(data) {
  cpue_data <- data$cpue_data
  n_index <- data$n_index
  if (is.null(n_index)) {
    n_index <- if (is.data.frame(cpue_data) && nrow(cpue_data)) {
      if ("index" %in% names(cpue_data)) max(cpue_data$index, na.rm = TRUE) else 1L
    } else 1L
  }
  n_fishery <- data$n_fishery
  if (is.null(n_fishery)) {
    fisheries <- unlist(data[c("removal_switch_f", "lf_fishery_f", "wf_fishery_f",
                               "lf_fishery", "wf_fishery")], use.names = FALSE)
    n_fishery <- if (length(fisheries)) max(c(length(data$removal_switch_f), fisheries),
      na.rm = TRUE) else 1L
  }
  for (name in c("n_index", "n_fishery")) {
    size <- if (name == "n_index") n_index else n_fishery
    if (!is.numeric(size) || length(size) != 1L || !is.finite(size) ||
        size < 1L || size != floor(size)) {
      stop("`", name, "` must be a positive integer.", call. = FALSE)
    }
  }
  if (!is.null(cpue_data) && is.data.frame(cpue_data) &&
      "index" %in% names(cpue_data) &&
      (any(!is.finite(cpue_data$index)) || any(cpue_data$index != floor(cpue_data$index)) ||
       any(cpue_data$index < 1L | cpue_data$index > n_index))) {
    stop("CPUE indices must be integers between 1 and `n_index`.", call. = FALSE)
  }
  lambdas <- data$lambdas
  if (is.null(lambdas)) lambdas <- list()
  .opal_validate_lambdas(lambdas, n_index, n_fishery)
  defaults <- list(cpue = 1, lf = 1, wf = 1, rec = 1, init_rec = 1,
                   prior = 1, penalty = 1, prior_mapped_only = FALSE)
  result <- defaults
  for (name in names(lambdas)) {
    if (!is.null(lambdas[[name]])) result[[name]] <- lambdas[[name]]
  }
  result$cpue <- rep(result$cpue, length.out = n_index)
  result$lf <- rep(result$lf, length.out = n_fishery)
  result$wf <- rep(result$wf, length.out = n_fishery)
  result
}

.opal_validate_lambdas <- function(lambdas, n_index, n_fishery) {
  if (!is.list(lambdas) || (length(lambdas) &&
      (is.null(names(lambdas)) || anyNA(names(lambdas)) ||
       any(!nzchar(names(lambdas))) || anyDuplicated(names(lambdas))))) {
    stop("`lambdas` must be a named list.", call. = FALSE)
  }
  allowed <- c("cpue", "lf", "wf", "rec", "init_rec", "prior", "penalty",
               "prior_mapped_only")
  unknown <- setdiff(names(lambdas), allowed)
  if (length(unknown)) {
    stop("Unknown lambda name: ", paste(unknown, collapse = ", "), ".", call. = FALSE)
  }
  sizes <- c(cpue = n_index, lf = n_fishery, wf = n_fishery,
             rec = 1L, init_rec = 1L, prior = 1L, penalty = 1L)
  for (name in intersect(names(lambdas), names(sizes))) {
    value <- lambdas[[name]]
    if (is.null(value)) next
    if (!is.numeric(value) || !length(value) || any(!is.finite(value)) ||
        any(value < 0) || !length(value) %in% c(1L, sizes[[name]])) {
      stop("Lambda `", name, "` must be finite, non-negative, and have length one or ",
           sizes[[name]], ".", call. = FALSE)
    }
  }
  if ("prior_mapped_only" %in% names(lambdas)) {
    value <- lambdas$prior_mapped_only
    if (!is.null(value) &&
        (!is.logical(value) || length(value) != 1L || is.na(value))) {
      stop("Lambda `prior_mapped_only` must be TRUE or FALSE.", call. = FALSE)
    }
  }
  invisible(lambdas)
}

.opal_weight_likelihood <- function(lp, lambda) {
  if (!length(lp)) return(lp)
  if (length(lambda) == 1L) {
    weights <- rep(lambda, length(lp))
  } else if (length(lambda) == length(lp)) {
    weights <- lambda
  } else {
    return(lp)
  }
  if (all(weights == 1)) return(lp)
  if (all(weights == 0)) return(rep(0, length(lp)))
  "[<-" <- RTMB::ADoverload("[<-")
  out <- lp
  zero <- which(weights == 0)
  if (length(zero)) out[zero] <- 0
  scaled <- which(weights != 0 & weights != 1)
  if (length(scaled)) out[scaled] <- weights[scaled] * lp[scaled]
  out
}

.opal_resolve_prior_lambdas <- function(priors, map, parameters) {
  for (i in seq_along(priors)) {
    idx <- priors[[i]]$index
    name <- names(parameters)[idx]
    mapped <- map[[name]]
    if (is.null(mapped)) next
    active <- !is.na(as.vector(mapped))
    lambda <- priors[[i]]$lambda
    if (is.null(lambda)) lambda <- rep(1, length(active))
    priors[[i]]$lambda <- rep(lambda, length.out = length(active)) * active
  }
  priors
}

.opal_lambda_tempered <- function(data, include_prior = FALSE) {
  lambdas <- get_lambdas(data)
  numeric_lambdas <- lambdas[setdiff(names(lambdas), "prior_mapped_only")]
  any(vapply(numeric_lambdas, function(x) any(x != 1), logical(1))) ||
    (include_prior && any(vapply(data$priors %||% list(), function(p) {
      !is.null(p$lambda) && any(p$lambda != 1)
    }, logical(1))))
}

.opal_warn_tempered <- function(data, include_prior = FALSE, context = "objective") {
  if (.opal_lambda_tempered(data, include_prior)) {
    warning("Likelihood/prior lambdas other than 1 temper the ", context, ".", call. = FALSE)
  }
  invisible(NULL)
}

.opal_warn_lambda_identifiability <- function(data, parameters, map) {
  lambdas <- get_lambdas(data)
  issues <- character()
  active <- function(name, i) {
    if (is.null(parameters[[name]])) return(FALSE)
    j <- if (length(parameters[[name]]) == 1L) 1L else i
    if (j > length(parameters[[name]])) return(FALSE)
    mapped <- map[[name]]
    is.null(mapped) || (j <= length(mapped) && !is.na(as.vector(mapped)[j]))
  }
  for (i in which(lambdas$cpue == 0)) {
    if (active("log_cpue_q", i) || active("log_cpue_tau", i)) {
      issues <- c(issues, paste0("CPUE index ", i))
    }
  }
  for (i in which(lambdas$lf == 0)) {
    if (active("log_lf_tau", i)) issues <- c(issues, paste0("LF fleet ", i))
  }
  for (i in which(lambdas$wf == 0)) {
    if (active("log_wf_tau", i)) issues <- c(issues, paste0("WF fleet ", i))
  }
  if (length(issues)) {
    warning("Zero likelihood lambdas may leave estimated parameters unidentified (",
            paste(unique(issues), collapse = ", "),
            "); map the corresponding parameters off.", call. = FALSE)
  }
  invisible(issues)
}

#' Set a prior lambda
#'
#' Sets the contribution weight for all elements of a prior block or for one
#' selected element. For an `opal_obj`, the updated object is returned through
#' [opal_update()] so dependent results are invalidated.
#' @param x An `opal_obj` or a named priors list.
#' @param parameter Name of the parameter block.
#' @param value Non-negative finite lambda value(s).
#' @param element Optional one-based element index within the block.
#' @return The updated `opal_obj` or priors list.
#' @export
opal_set_prior_lambda <- function(x, parameter, value, element = NULL) {
  object <- inherits(x, "opal_obj")
  priors <- if (object) x$data$priors else x
  if (!is.list(priors) || is.null(names(priors)) || !parameter %in% names(priors)) {
    stop("No prior is defined for `", parameter, "`.", call. = FALSE)
  }
  if (!is.numeric(value) || !length(value) || any(!is.finite(value)) || any(value < 0)) {
    stop("`value` must be finite and non-negative.", call. = FALSE)
  }
  prior <- priors[[parameter]]
  known_length <- max(length(prior$par1), length(prior$par2),
                      length(prior$lambda), 1L)
  parameter_value <- if (object) x$parameters[[parameter]] else NULL
  parameter_length_known <- length(parameter_value) > 0L
  block_length <- if (length(parameter_value)) length(parameter_value) else known_length
  if (is.null(element)) {
    if (!length(value) %in% c(1L, block_length) &&
        (parameter_length_known || known_length > 1L)) {
      stop("`value` must have length one or the prior block length.", call. = FALSE)
    }
    prior$lambda <- value
  } else {
    if (!is.numeric(element) || length(element) != 1L || !is.finite(element) ||
        element != floor(element) || element < 1L ||
        ((parameter_length_known || known_length > 1L) && element > block_length) ||
        length(value) != 1L) {
      stop("`element` must select one prior element and `value` must be scalar.", call. = FALSE)
    }
    lambda <- if (is.null(prior$lambda)) rep(1, max(block_length, element)) else
      rep(prior$lambda, length.out = max(block_length, element))
    lambda[element] <- value
    prior$lambda <- lambda
  }
  priors[[parameter]] <- prior
  if (object) opal_update(x, priors = priors) else priors
}
