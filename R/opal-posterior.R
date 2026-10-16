.opal_store_derived <- function(x, name, result, scope, settings) {
  if (!is.character(name) || length(name) != 1L || is.na(name) || !nzchar(name)) stop("Supply one result name.")
  record <- list(result = result, scope = scope, identity = .opal_check_identity(x, scope),
                 settings = settings, version = 1L)
  record$payload_id <- .opal_obj_hash(record)
  x$derived[[name]] <- record
  .opal_obj_seal(x)
}

#' Retrieve a stored assessment analysis
#' @param x An `opal_obj`.
#' @param name Name of the analysis under `x$derived`.
#' @return The stored result after checking source identity and integrity.
#' @family assessment workflow
#' @export
opal_derived <- function(x, name) {
  validate_opal_obj(x, results = TRUE)
  record <- x$derived[[name]]
  if (is.null(record)) stop("No stored analysis named '", name, "'.", call. = FALSE)
  if (is.null(record$scope) || !identical(record$identity, .opal_check_identity(x, record$scope)) ||
      !identical(record$payload_id, .opal_obj_hash(record[setdiff(names(record), "payload_id")]))) {
    stop("Stored analysis is stale or has been modified; recompute it.", call. = FALSE)
  }
  record$result
}

#' Summarise posterior parameters and derived model quantities
#'
#' Evaluates each selected joint draw using its own biological inputs and
#' reports. Summaries are stored in the portable assessment object.
#' @param x An `opal_obj` with posterior draws.
#' @param quantities Numeric report names to summarise.
#' @param probs Distinct, increasing quantile probabilities between zero and one.
#' @param draws Optional retained draw indices, in iteration-within-chain order.
#' @param name Name of the stored result.
#' @return An updated object. [opal_derived()] returns long-form `parameters`
#'   and `reports` summaries, report dimensions, draw identifiers, and the
#'   source MCMC validation status. Element indices follow R's column order.
#' @details Summaries never establish posterior acceptance. Non-finite report
#'   values cause an error, rather than silently discarding draws. Marginal
#'   posterior samples lacking random effects cannot be used. Use
#'   [opal_check()] to assess mixing, sampler diagnostics, bounds, and biology.
#' @family assessment workflow
#' @export
opal_posterior <- function(x,
    quantities = c("spawning_biomass_y", "static_depletion_y", "dynamic_depletion_y"),
    probs = c(0.025, 0.5, 0.975), draws = NULL, name = "posterior") {
  if (!is.numeric(probs) || !length(probs) || any(!is.finite(probs)) ||
      any(probs < 0 | probs > 1) || is.unsorted(probs, strictly = TRUE)) stop("Invalid quantile probabilities.")
  if (!is.character(quantities) || !length(quantities) || anyDuplicated(quantities)) stop("Supply distinct report names.")
  context <- .opal_posterior_context(x, draws)
  first <- context$object$report(context$values[1L, ])
  if (!all(quantities %in% names(first)) || any(!vapply(first[quantities], is.numeric, logical(1)))) {
    stop("All requested quantities must be numeric model reports.", call. = FALSE)
  }
  dimensions <- lapply(first[quantities], function(z) if (is.null(dim(z))) length(z) else dim(z))
  reports <- lapply(first[quantities], function(z) matrix(NA_real_, nrow(context$values), length(z)))
  for (i in seq_len(nrow(context$values))) {
    r <- if (i == 1L) first else context$object$report(context$values[i, ])
    for (quantity in quantities) {
      z <- r[[quantity]]
      if (length(z) != ncol(reports[[quantity]]) || any(!is.finite(z))) {
        stop("Invalid report at retained draw ", context$ids$draw[i], ": ", quantity, call. = FALSE)
      }
      reports[[quantity]][i, ] <- as.numeric(z)
    }
  }
  summarise <- function(values, quantity) {
    out <- data.frame(quantity = quantity, element = seq_len(ncol(values)),
      mean = colMeans(values), sd = apply(values, 2L, stats::sd))
    q <- vapply(seq_len(ncol(values)), function(i) stats::quantile(values[, i], probs,
      names = FALSE), numeric(length(probs)))
    q <- t(matrix(q, nrow = length(probs)))
    colnames(q) <- paste0("q", format(probs, trim = TRUE, scientific = FALSE))
    out <- cbind(out, q)
    rownames(out) <- NULL
    out
  }
  parameter_summary <- summarise(context$values, "parameter")
  parameter_summary$parameter <- colnames(context$values)
  result <- list(parameters = parameter_summary,
    reports = do.call(rbind, Map(summarise, reports, names(reports))),
    dimensions = dimensions, draws = context$ids, probs = probs,
    validation = summary(x)$checks$mcmc)
  .opal_store_derived(x, name, result, "mcmc", list(quantities = quantities, probs = probs, draws = draws))
}

#' Plot observed and fitted compositions
#' @param x A fitted `opal_obj`.
#' @param type Length (`"lf"`) or weight (`"wf"`) compositions.
#' @param fishery Optional fishery indices to include.
#' @return A ggplot with observation proportions and fitted probabilities.
#'   Removed and zero-size rows are excluded. Bin numbers refer to the prepared
#'   data, including any tail aggregation.
#' @family assessment diagnostics
#' @export
plot_composition <- function(x, type = c("lf", "wf"), fishery = NULL) {
  type <- match.arg(type)
  object <- opal_rtmb(x, fresh = TRUE)
  p <- object$env$last.par.best
  rows <- .opal_composition_rows(x$data, object$env$parList(par = p), object$report(p), type)
  rows <- Filter(function(z) z$include && (is.null(fishery) || z$fishery %in% fishery), rows)
  if (!length(rows)) stop("No included compositions match this selection.")
  df <- do.call(rbind, lapply(rows, function(z) data.frame(fishery = z$fishery,
    year = z$year, row = z$row, bin = z$bin, observed = z$observed / sum(z$observed), predicted = z$predicted)))
  ggplot2::ggplot(df, ggplot2::aes(x = .data$bin, y = .data$observed)) +
    ggplot2::geom_point() + ggplot2::geom_line(ggplot2::aes(y = .data$predicted)) +
    ggplot2::facet_wrap(~fishery + year + row) +
    ggplot2::labs(x = paste(if (type == "lf") "Length" else "Weight", "bin"),
      y = "Proportion", caption = "Observed points; fitted line") + ggplot2::theme_bw()
}

#' Compare parameter priors and posterior distributions
#' @param x An `opal_obj` with joint posterior samples and `data$priors`.
#' @param parameters Optional expanded active parameter names to plot.
#' @return A ggplot with marginal posterior histograms and prior densities on
#'   each parameter's stored scale. Fixed and shared mapped elements are
#'   omitted when there is no unambiguous scalar prior correspondence.
#' @details The displayed priors are the specified, untruncated densities.
#'   Optimiser/sampler bounds may truncate the realised target. For shared
#'   parameter maps, multiple prior contributions need not equal one density;
#'   such parameters and blocks with multiple priors are excluded with a warning.
#' @family assessment diagnostics
#' @export
plot_prior_posterior <- function(x, parameters = NULL) {
  context <- .opal_posterior_context(x)
  priors <- x$data$priors
  if (!length(priors)) stop("No parameter priors are specified.")
  .opal_validate_priors(x$parameters, priors)
  frames <- curves <- list()
  skipped <- character()
  for (prior in priors) {
    block <- names(x$parameters)[prior$index]
    if (sum(vapply(priors, function(z) z$index == prior$index, logical(1))) > 1L) {
      skipped <- c(skipped, block)
      next
    }
    base <- x$parameters[[block]]
    map <- x$map[[block]]
    if (is.null(map)) map <- factor(seq_along(base))
    active <- which(!is.na(map))
    if (anyDuplicated(map[active])) { skipped <- c(skipped, block); next }
    # RTMB orders active factor levels, then expands repeated block names.
    active <- active[order(as.integer(map[active]))]
    names <- .opal_expand_parameter_names(rep(block, length(active)))
    for (j in seq_along(active)) {
      label <- names[j]
      if (!label %in% colnames(context$values) || (!is.null(parameters) && !label %in% parameters)) next
      i <- active[j]
      par1 <- rep(prior$par1, length.out = length(base))[i]
      par2 <- rep(prior$par2, length.out = length(base))[i]
      values <- context$values[, label]
      lim <- range(values)
      pad <- max(diff(lim) * 0.1, .Machine$double.eps^0.5)
      grid <- seq(lim[1] - pad, lim[2] + pad, length.out = 300)
      density <- switch(prior$type,
        normal = stats::dnorm(grid, par1, par2),
        student = stats::dt((grid - par1) / par2, 3) / par2,
        lognormal = stats::dlnorm(grid, par1, par2),
        beta = stats::dbeta(grid, par1 * par2, (1 - par1) * par2))
      frames[[label]] <- data.frame(parameter = label, value = values)
      curves[[label]] <- data.frame(parameter = label, value = grid, density = density)
    }
  }
  if (length(skipped)) warning("Ambiguous prior blocks omitted: ", paste(unique(skipped), collapse = ", "), call. = FALSE)
  if (!length(frames)) stop("No unambiguous sampled parameters match the priors and selection.")
  ggplot2::ggplot(do.call(rbind, frames), ggplot2::aes(x = .data$value)) +
    ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)), bins = 30,
      fill = "grey80", colour = "white") +
    ggplot2::geom_line(data = do.call(rbind, curves), ggplot2::aes(y = .data$density), colour = "#0072B2") +
    ggplot2::facet_wrap(~parameter, scales = "free") +
    ggplot2::labs(x = "Stored parameter scale", y = "Density", caption = "Posterior histogram; prior line") +
    ggplot2::theme_bw()
}
