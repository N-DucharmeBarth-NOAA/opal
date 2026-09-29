#' Profile an assessment parameter
#'
#' Fixes a scalar parameter or a mapped parameter group at each supplied value,
#' then re-optimises the remaining parameters using the usual two-pass fit.
#' @param x A fitted `opal_obj`.
#' @param parameter Name of a parameter block, on its stored scale.
#' @param values Distinct finite profile values on that scale.
#' @param element One-based element within the block, in R column order.
#' @param fit_args Named arguments passed to [opal_fit()].
#' @param name Name of the stored result.
#' @return An updated object containing a table, component contributions,
#'   conditional fits, and failed-point messages, accessible with [opal_derived()].
#' @details This profiles the complete penalised objective, including priors
#'   and process penalties; it is a likelihood profile only when those terms
#'   are absent. Fixed effects can be profiled in a Laplace model, but random
#'   effect blocks cannot. Shared map elements are fixed together. The
#'   original fit is the reference; a negative objective difference indicates
#'   that a conditional fit improved on it. Failed points are never interpolated.
#' @family assessment sensitivity
#' @export
opal_profile <- function(x, parameter, values, element = 1L, fit_args = list(), name = "profile") {
  validate_opal_obj(x)
  if (is.null(x$fit$opt)) stop("Profiling requires a fitted model.")
  if (length(parameter) != 1L || !parameter %in% names(x$parameters) || parameter %in% x$random) {
    stop("Choose one fixed-effect parameter block.")
  }
  if (length(element) != 1L || !is.finite(element) || element != floor(element) ||
      element < 1 || element > length(x$parameters[[parameter]])) stop("Invalid parameter element.")
  if (!is.numeric(values) || !length(values) || any(!is.finite(values)) || anyDuplicated(values)) stop("Supply distinct finite profile values.")
  .opal_validate_named_list(fit_args, "fit_args")
  if (length(intersect(names(fit_args), c("data", "x", "check")))) stop("Profile controls cannot override the model or checks.")
  mapping <- x$map[[parameter]]
  if (is.null(mapping)) mapping <- factor(seq_along(x$parameters[[parameter]]))
  group <- as.character(mapping[element])
  if (is.na(group)) stop("The profiled element is already fixed.")
  members <- which(as.character(mapping) == group)
  active <- which(names(x$fit$opt$par) == parameter)
  removed <- active[match(group, levels(mapping))]
  if (any(values < x$bounds$lower[removed] | values > x$bounds$upper[removed])) stop("Profile values lie outside the configured bounds.")
  fits <- vector("list", length(values)); table <- components <- vector("list", length(values))
  for (i in seq_along(values)) {
    attempt <- tryCatch({
      p <- x$fit$parameters; p[[parameter]][members] <- values[i]
      map <- x$map; mapping[members] <- NA; map[[parameter]] <- droplevels(mapping)
      candidate <- opal_update(x, parameters = p, map = map,
        bounds = list(lower = unname(x$bounds$lower[-removed]), upper = unname(x$bounds$upper[-removed])))
      candidate <- opal_build(candidate)
      object <- opal_rtmb(candidate)
      if (!length(object$par)) {
        candidate <- opal_attach_fit(candidate, list(par = object$par,
          objective = object$fn(object$par), convergence = 0L), check = FALSE)
      } else candidate <- do.call(opal_fit, c(list(candidate, check = FALSE), fit_args))
      object <- opal_rtmb(candidate)
      gradient <- if (length(object$par)) max(abs(object$gr(candidate$fit$opt$par))) else 0
      biology <- opal_diagnose(candidate)
      report <- opal_report(candidate)
      terms <- intersect(c("lp_prior", "lp_penalty", "lp_rec", "lp_init_rec", "lp_cpue", "lp_lf", "lp_wf"), names(report))
      components[[i]] <- data.frame(value = values[i], component = terms,
        objective = vapply(report[terms], sum, numeric(1)))
      fits[[i]] <- candidate
      data.frame(value = values[i], objective = candidate$fit$opt$objective,
        delta = candidate$fit$opt$objective - x$fit$opt$objective,
        max_gradient = gradient, passes = candidate$fit$opt$convergence == 0 &&
          is.finite(gradient) && gradient <= 1e-3 && biology$passes, error = NA_character_)
    }, error = function(e) data.frame(value = values[i], objective = NA_real_, delta = NA_real_,
      max_gradient = NA_real_, passes = FALSE, error = conditionMessage(e)))
    table[[i]] <- attempt
  }
  result <- list(table = do.call(rbind, table), components = do.call(rbind, components),
    fits = fits, parameter = parameter, element = element, mapped_elements = members,
    reference_objective = x$fit$opt$objective)
  .opal_store_derived(x, name, result, "fit", list(values = values, fit_args = fit_args))
}

#' Plot an objective profile
#' @param x An Opal object with a stored profile.
#' @param name Name of the stored profile.
#' @return A ggplot. Failed conditional fits are shown as crosses when their
#'   objectives are available; errors remain in the stored table.
#' @family assessment sensitivity
#' @export
plot_opal_profile <- function(x, name = "profile") {
  result <- opal_derived(x, name)
  ggplot2::ggplot(result$table, ggplot2::aes(x = .data$value, y = .data$delta)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey60") +
    ggplot2::geom_point(ggplot2::aes(shape = .data$passes), size = 2.5) +
    ggplot2::scale_shape_manual(values = c(`FALSE` = 4, `TRUE` = 16)) +
    ggplot2::labs(x = paste0(result$parameter, "[", result$element, "]"),
      y = "Change in penalised objective", shape = "Checks pass") + ggplot2::theme_bw()
}

#' Fit a reproducible grid of assessment scenarios
#' @param x A configured or fitted `opal_obj`.
#' @param scenarios Named list of argument lists passed to [opal_update()].
#' @param directory Optional directory for per-scenario portable checkpoints.
#' @param resume Reuse checkpoints only when target and fitting settings match.
#' @param fit_args Named arguments passed to [opal_fit()].
#' @return An `opal_grid` list with `models`, `summary`, and `settings`.
#'   Every scenario is retained, including failures. `passes` uses the full
#'   fit check, not just the optimiser code. Checkpoint identities include the
#'   scientific configuration and fitting controls.
#' @family assessment sensitivity
#' @export
opal_grid <- function(x, scenarios, directory = NULL, resume = TRUE, fit_args = list()) {
  validate_opal_obj(x)
  .opal_validate_named_list(scenarios, "scenarios")
  .opal_validate_named_list(fit_args, "fit_args")
  if (!length(scenarios) || anyDuplicated(names(scenarios))) stop("Supply distinct named scenarios.")
  if (length(intersect(names(fit_args), c("data", "x", "check")))) stop("Grid controls cannot override model inputs or checks.")
  .opal_obj_flag(resume, "resume")
  if (!is.null(directory)) dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  models <- stats::setNames(vector("list", length(scenarios)), names(scenarios))
  rows <- vector("list", length(scenarios))
  for (i in seq_along(scenarios)) {
    label <- names(scenarios)[i]
    reused <- FALSE
    attempt <- tryCatch({
      args <- scenarios[[i]]
      .opal_validate_named_list(args, "scenario")
      if ("x" %in% names(args)) stop("A scenario cannot replace x.")
      candidate <- do.call(opal_update, c(list(x), args))
      candidate <- opal_build(candidate)
      key <- .opal_obj_hash(list(identity = candidate$identity, control = candidate$control,
        fit_args = fit_args, validation = .opal_validation_version))
      path <- if (!is.null(directory)) file.path(directory, paste0(key, ".rds")) else NULL
      if (resume && !is.null(path) && file.exists(path)) {
        saved <- readRDS(path)
        if (!identical(saved$key, key) || !identical(saved$payload_id, .opal_obj_hash(saved$model))) stop("Invalid grid checkpoint.")
        validate_opal_obj(saved$model, results = TRUE)
        candidate <- saved$model
        reused <- TRUE
      } else {
        candidate <- do.call(opal_fit, c(list(candidate, check = TRUE), fit_args))
        if (!is.null(path)) {
          saved <- list(key = key, model = candidate, payload_id = .opal_obj_hash(candidate))
          temporary <- tempfile(tmpdir = directory)
          on.exit(unlink(temporary), add = TRUE)
          saveRDS(saved, temporary)
          if (!file.rename(temporary, path)) stop("Could not write grid checkpoint.")
        }
      }
      models[[i]] <- candidate
      data.frame(model = label, objective = candidate$fit$opt$objective,
        passes = identical(summary(candidate)$checks$fit, "passed"), reused = reused, error = NA_character_)
    }, error = function(e) data.frame(model = label, objective = NA_real_, passes = FALSE,
      reused = FALSE, error = conditionMessage(e)))
    rows[[i]] <- attempt
  }
  structure(list(models = models, summary = do.call(rbind, rows),
    settings = list(scenarios = scenarios, fit_args = fit_args, directory = directory)), class = "opal_grid")
}

#' Sample accepted members of an assessment grid
#' @param grid An `opal_grid` returned by [opal_grid()].
#' @param ... Arguments passed to [opal_mcmc()]. Use an explicit seed.
#' @return The grid with updated models and a `mcmc_passes` summary column.
#'   Failed fits are skipped, and failed sampling attempts remain in each model.
#' @family assessment sensitivity
#' @export
opal_grid_mcmc <- function(grid, ...) {
  if (!inherits(grid, "opal_grid")) stop("Supply an opal_grid.")
  grid$summary$mcmc_passes <- FALSE
  for (i in seq_along(grid$models)) {
    model <- grid$models[[i]]
    if (is.null(model) || !identical(summary(model)$checks$fit, "passed")) next
    grid$models[[i]] <- opal_mcmc(model, ...)
    grid$summary$mcmc_passes[i] <- identical(summary(grid$models[[i]])$checks$mcmc, "passed")
  }
  grid
}

#' Select balanced posterior draws from a model grid
#' @param grid An `opal_grid` with accepted posteriors in every member.
#' @param n_per_model Number of draws selected without replacement per model.
#' @param seed Reproducible selection seed, restored on exit.
#' @return A table of model, chain, iteration, retained draw index, and equal
#'   model weights. Selection is balanced across chains within each model,
#'   with any remainder allocated to the first chains. Models are not silently
#'   dropped, and equal model weighting is an explicit assumption.
#' @family assessment sensitivity
#' @export
opal_grid_draws <- function(grid, n_per_model, seed = 123L) {
  if (!inherits(grid, "opal_grid") || !length(grid$models)) stop("Supply a non-empty opal_grid.")
  if (length(n_per_model) != 1 || !is.finite(n_per_model) || n_per_model < 1 || n_per_model != floor(n_per_model)) stop("Supply a positive integer draw count.")
  .opal_with_seed(seed, {
    out <- lapply(seq_along(grid$models), function(i) {
      x <- grid$models[[i]]
      if (is.null(x) || !identical(summary(x)$checks$mcmc, "passed")) stop("Every grid member requires a current passing posterior check.")
      ids <- .opal_posterior_context(x)$ids
      count <- rep(n_per_model %/% x$mcmc$chains, x$mcmc$chains)
      if (n_per_model %% x$mcmc$chains) count[seq_len(n_per_model %% x$mcmc$chains)] <- count[seq_len(n_per_model %% x$mcmc$chains)] + 1L
      selected <- unlist(lapply(seq_along(count), function(chain) {
        available <- which(ids$chain == chain)
        if (count[chain] > length(available)) stop("Insufficient retained draws for balanced selection.")
        available[sample.int(length(available), count[chain])]
      }))
      cbind(model = names(grid$models)[i], ids[selected, ], weight = 1 / (length(grid$models) * n_per_model))
    })
    do.call(rbind, out)
  })
}
