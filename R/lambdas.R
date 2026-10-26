#' Resolve likelihood-component lambdas
#'
#' Expands the optional `data$lambdas` list to one multiplier per CPUE index or
#' composition fleet. Scalar values are recycled; unspecified components use 1.
#' @details The supported entries are `cpue` (one or `n_index` values), `lf`
#'   and `wf` (one or `n_fishery` values each), and scalar `rec`, `init_rec`,
#'   `prior`, and `penalty` multipliers. `prior_mapped_only` is a logical flag
#'   that removes prior contributions for elements fixed by the parameter map.
#'
#'   `opal_model()` cannot see the parameter map, so `prior_mapped_only` is
#'   applied only when the objective is built from an `opal_obj` (through
#'   [opal_build()], [opal_rtmb()], [opal_fit()] or [opal_mcmc()]). An
#'   objective built directly with `RTMB::MakeADFun(cmb(opal_model, data), ...)`,
#'   or rebuilt from a legacy `opal_fit`, still includes priors on mapped
#'   elements, and `opal_model()` warns that the flag was not applied.
#' @param data A model data list.
#' @return A fully expanded named list of lambda values.
#' @export
get_lambdas <- function(data) {
  cpue_data <- data$cpue_data
  index <- if (is.data.frame(cpue_data) && "index" %in% names(cpue_data)) {
    cpue_data$index
  } else numeric()
  if (isTRUE(data$cpue_switch > 0) && length(index) &&
      (any(!is.finite(index)) || any(index != floor(index)) || any(index < 1L))) {
    # get_cpue_like() would silently drop these rows from the likelihood.
    stop("CPUE indices must be positive integers.", call. = FALSE)
  }
  # As in get_cpue_like(), an index above `n_index` extends the index count.
  valid <- is.finite(index) & index == floor(index) & index >= 1L
  n_index <- max(c(data$n_index, index[valid], 1L))
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
    # A mismatch means the observation order changed; never drop a lambda silently.
    stop("Lambda length (", length(lambda), ") does not match the number of ",
         "likelihood terms (", length(lp), ").", call. = FALSE)
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
  # opal_model() warns when prior_mapped_only is set without this marker.
  attr(priors, "opal_mapped_only_resolved") <- TRUE
  priors
}

# Data passed to MakeADFun by the staged workflow, with prior_mapped_only
# applied. The model function cannot see the map, so this is resolved here.
.opal_runtime_data <- function(data, map, parameters) {
  if (isTRUE(get_lambdas(data)$prior_mapped_only) && length(data$priors)) {
    data$priors <- .opal_resolve_prior_lambdas(data$priors, map, parameters)
  }
  data
}

# TRUE when prior_mapped_only would remove a prior contribution, that is when
# a prior covers at least one element fixed by the map.
.opal_prior_mapped_only_changes <- function(data, map, parameters) {
  if (!isTRUE(get_lambdas(data)$prior_mapped_only) || !length(data$priors)) {
    return(FALSE)
  }
  resolved <- .opal_resolve_prior_lambdas(data$priors, map, parameters)
  any(vapply(seq_along(resolved), function(i) {
    after <- resolved[[i]]$lambda
    if (is.null(after)) return(FALSE)
    before <- data$priors[[i]]$lambda
    if (is.null(before)) before <- 1
    any(rep(before, length.out = length(after)) != after)
  }, logical(1)))
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
