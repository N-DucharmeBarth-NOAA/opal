#' Update an Opal model configuration
#'
#' Changes to the scientific configuration clear dependent fits, samples, and
#' reports. Metadata and optimiser controls do not alter the scientific target.
#' Defaults are resolved again when their upstream configuration changes.
#' @param x An `opal_obj`.
#' @param data,parameters,map,random,bounds,control,makeadfun_args As in [opal_obj()].
#' @param priors Optional replacement for `data$priors`.
#' @param metadata Named metadata to merge.
#' @return An updated `opal_obj`; the input object is not modified.
#' @details
#' Changes to data, priors, parameters, maps, random effects, bounds, or RTMB
#' construction settings clear dependent fits, posteriors, diagnostics, and
#' projections. Metadata is merged, and optimiser controls can change without
#' invalidating results. A supplied map is retained unless explicitly replaced;
#' check its dimensions when changing parameter structure.
#' @family assessment workflow
#' @examples
#' inputs <- opaka_quickstart_inputs()
#' assessment <- opal_obj(inputs$data, inputs$parameters, inputs$map)
#' labelled <- opal_update(assessment, metadata = list(label = "Baseline"))
#' summary(labelled)
#' @export
opal_update <- function(x, data, parameters, map, random, bounds, control,
                        makeadfun_args, priors, metadata) {
  validate_opal_obj(x)
  old <- x
  supplied <- as.list(match.call())[-1L]
  supplied$x <- NULL
  configuration <- c("data", "parameters", "map", "random", "bounds", "makeadfun_args")
  for (name in intersect(names(supplied), configuration)) {
    value <- get(name, inherits = FALSE)
    if (name == "data") value <- .opal_sanitize_tables(value)
    x[name] <- list(value)
    if (name %in% c("parameters", "map", "bounds")) {
      x$configuration$origins[[name]] <- if (is.null(value)) "unresolved" else "supplied"
    }
  }
  if (!missing(priors)) x$data["priors"] <- list(priors)
  if (!missing(control)) x["control"] <- list(control)
  if (!missing(metadata)) {
    .opal_validate_named_list(metadata, "metadata")
    x$provenance$metadata <- utils::modifyList(x$provenance$metadata, metadata)
  }
  upstream <- !identical(old$data, x$data) || !identical(old$parameters, x$parameters)
  if (upstream && missing(map) && old$configuration$origins$map != "supplied") {
    x["map"] <- list(NULL)
    x$configuration$origins$map <- "unresolved"
  }
  if (!identical(old$data, x$data) && missing(parameters) &&
      old$configuration$origins$parameters == "default") {
    x["parameters"] <- list(NULL)
    x$configuration$origins$parameters <- "unresolved"
  }
  if ((upstream || !identical(old$map, x$map) || !identical(old$random, x$random)) &&
      missing(bounds) && old$configuration$origins$bounds != "supplied") {
    x["bounds"] <- list(NULL)
    x$configuration$origins$bounds <- "unresolved"
  }
  if (!identical(.opal_obj_identity(x)$target, old$identity$target)) {
    x["build"] <- list(NULL)
    x$fit <- list(opt = NULL, parameters = NULL, last_par_best = NULL,
                  diagnostics = list(), estimability = NULL)
    x["mcmc"] <- list(NULL)
    x$mcmc_history <- list()
    x$derived <- list()
    x$validation <- list()
  }
  .opal_obj_portable(unclass(x), "opal_obj")
  .opal_obj_seal(x)
}

#' Attach an externally optimised Opal fit
#' @param x A configured `opal_obj`.
#' @param opt Optimiser result with named `par`, `objective`, and `convergence`.
#' @param check Run fitting diagnostics after attachment.
#' @param check_args Named arguments passed to [opal_check()].
#' @return A fitted `opal_obj`.
#' @description
#' Attach a result from an external optimiser to the configuration that produced
#' it. The objective, parameter layout, and bounds are verified before retaining
#' the point. For ordinary optimisation, use [opal_fit()].
#' @family assessment workflow
#' @examples
#' assessment <- opal_read(system.file("extdata", "opaka_quickstart_fit.rds",
#'                                    package = "opal"))
#' # Reattach the already verified optimum; no optimisation is performed.
#' assessment <- opal_attach_fit(assessment, assessment$fit$opt, check = FALSE)
#' @export
opal_attach_fit <- function(x, opt, check = TRUE, check_args = list()) {
  x <- opal_build(x)
  .opal_obj_flag(check, "check")
  object <- opal_rtmb(x, fresh = TRUE)
  # Capture exactly the supplied point, even when it is worse than an earlier
  # evaluation. RTMB's best-so-far state must not silently substitute a point.
  object$env$value.best <- Inf
  # Reuse the legacy constructor's layout, bounds, and objective checks.
  captured <- opal_fit(data = x$data, obj = object, opt = opt,
                       bounds = x$bounds, control = x$control,
                       makeadfun_args = x$makeadfun_args)
  x$fit <- list(opt = captured$fit$opt, parameters = captured$parameters,
                last_par_best = object$env$last.par.best,
                diagnostics = list(), estimability = NULL)
  object$par <- captured$fit$opt$par
  object$opt <- captured$fit$opt
  x$derived <- list()
  x$validation$fit <- NULL
  x <- .opal_obj_seal(x)
  assign(.opal_obj_runtime_id(x), object, .opal_obj_cache)
  if (check) x <- .opal_run_check(x, "fit", check_args)
  x
}

.opal_optimise <- function(x, n_passes = 2L, control = NULL,
                           check = TRUE, check_args = list()) {
  validate_opal_obj(x)
  .opal_obj_flag(check, "check")
  if (length(n_passes) != 1L || !is.numeric(n_passes) ||
      !is.finite(n_passes) || n_passes < 1 || n_passes != as.integer(n_passes)) {
    stop("`n_passes` must be a positive integer.", call. = FALSE)
  }
  if (!is.null(control)) x <- opal_update(x, control = control)
  x <- opal_build(x)
  object <- opal_rtmb(x, fresh = TRUE)
  start <- object$par
  effective_control <- .opal_fit_or(x$control, list(eval.max = 10000L, iter.max = 10000L))
  passes <- vector("list", n_passes)
  for (pass in seq_len(n_passes)) {
    opt <- stats::nlminb(start, object$fn, object$gr,
                         lower = x$bounds$lower, upper = x$bounds$upper,
                         control = effective_control)
    passes[[pass]] <- opt
    start <- opt$par
  }
  x$control <- effective_control
  x <- opal_attach_fit(x, opt, check = FALSE)
  x$fit$diagnostics$optimisation <- list(n_passes = n_passes, passes = passes)
  if (check) x <- .opal_run_check(x, "fit", check_args)
  .opal_obj_seal(x)
}

.opal_check_identity <- function(x, scope) {
  if (scope == "fit") return(x$identity)
  list(target = x$identity$target, posterior = x$mcmc$payload_id)
}

.opal_run_check <- function(x, scope, args) {
  .opal_validate_named_list(args, "check_args")
  if (length(intersect(names(args), c("x", "scope")))) {
    stop("`check_args` cannot override x or scope.", call. = FALSE)
  }
  do.call(opal_check, c(list(x = x, scope = scope), args))
}

#' Check fitted or sampled Opal results
#'
#' Checks retain the object even when diagnostics fail. Validation records are
#' tied to the exact fitted or sampled state and are separate from lifecycle
#' stage. Short MCMC smoke tests are not evidence of convergence.
#' @param x An `opal_obj`.
#' @param scope Check the fitted optimum or MCMC.
#' @param gradient_tolerance Maximum absolute fitting gradient.
#' @param max_rhat Maximum rank-normalised R-hat.
#' @param min_ess Minimum bulk and tail effective sample sizes.
#' @param stop_on_failure Stop rather than warn when a check fails.
#' @param penalty_tolerance,catch_tolerance Biological thresholds as in [opal_diagnose()].
#' @return The object with an attached validation record.
#' @details
#' Fit checks cover optimiser convergence, maximum absolute gradient, a
#' positive-definite Hessian, parameter bounds, and finite non-negative numbers
#' at age, plus the biological checks in [opal_diagnose()]. Every retained
#' joint posterior draw is checked for bounds and biology. MCMC checks also
#' require at least two chains, finite R-hat and effective
#' sample sizes within the thresholds, known sampler diagnostics, no
#' divergences, and no maximum-tree-depth hits. Missing sampler diagnostics
#' prevent a passing MCMC check, even if imported parameter draws are usable.
#' The record is stored in `x$validation[[scope]]`, including metrics, settings,
#' and the identity of the checked result. These checks do not establish
#' scientific adequacy of a model or projection scenario.
#' @family assessment workflow
#' @examples
#' assessment <- opal_read(system.file("extdata", "opaka_quickstart_fit.rds",
#'                                    package = "opal"))
#' assessment <- opal_check(assessment, scope = "fit")
#' assessment$validation$fit$metrics
#' @export
opal_check <- function(x, scope = c("fit", "mcmc"), gradient_tolerance = 1e-3,
                       max_rhat = 1.01, min_ess = 100, stop_on_failure = FALSE,
                       penalty_tolerance = 1e-10, catch_tolerance = 1e-6) {
  scope <- match.arg(scope)
  validate_opal_obj(x, results = scope == "mcmc")
  .opal_obj_flag(stop_on_failure, "stop_on_failure")
  thresholds <- c(gradient_tolerance, max_rhat, min_ess)
  if (length(thresholds) != 3L || any(!is.finite(thresholds)) || any(thresholds <= 0)) {
    stop("Diagnostic thresholds must be positive finite scalars.", call. = FALSE)
  }
  if (scope == "fit" && is.null(x$fit$opt)) stop("No fitted optimum to check.", call. = FALSE)
  if (scope == "mcmc" && is.null(x$mcmc)) stop("No posterior to check.", call. = FALSE)
  result <- tryCatch({
    if (scope == "fit") {
      object <- opal_rtmb(x, fresh = TRUE)
      par <- x$fit$opt$par
      gradient <- object$gr(par)
      hessian <- if (is.function(object$he)) object$he(par) else stats::optimHess(par, object$fn, object$gr)
      positive <- !inherits(try(chol(hessian), silent = TRUE), "try-error")
      report <- object$report(x$fit$last_par_best)
      object$env$last.par.best <- x$fit$last_par_best
      estimability <- tryCatch(check_estimability(object, hessian), error = identity)
      x$fit$estimability <- .opal_compact_estimability(estimability,
        .opal_expand_parameter_names(names(par)), names(par))
      metrics <- list(convergence = x$fit$opt$convergence,
        max_gradient = max(abs(gradient)), positive_hessian = positive,
        inside_bounds = all(par >= x$bounds$lower & par <= x$bounds$upper),
        valid_population = all(is.finite(report$number_ysa)) && all(report$number_ysa >= 0))
      metrics$biology <- .opal_diagnose_report(x$data, x$fit$parameters, report,
        penalty_tolerance, catch_tolerance)
      passes <- identical(as.integer(metrics$convergence), 0L) &&
        is.finite(metrics$max_gradient) && metrics$max_gradient <= gradient_tolerance &&
        positive && metrics$inside_bounds && metrics$valid_population && metrics$biology$passes
    } else {
      if (!requireNamespace("posterior", quietly = TRUE)) {
        stop("Install the suggested posterior package to check MCMC.")
      }
      m <- x$mcmc
      keep <- seq.int(m$warmup + 1L, dim(m$samples)[1L])
      draws <- m$samples[keep, , m$par_names, drop = FALSE]
      diagnostics <- posterior::summarise_draws(posterior::as_draws_array(draws),
        "rhat", "ess_bulk", "ess_tail")
      divergent <- depth_hits <- 0L
      sampler_known <- is.list(m$sampler_params) && length(m$sampler_params) == m$chains &&
        length(m$max_treedepth) == 1L && is.finite(m$max_treedepth) && m$max_treedepth > 0
      for (params in m$sampler_params) {
        params <- as.matrix(params)
        if (nrow(params) == dim(m$samples)[1L]) params <- params[keep, , drop = FALSE]
        required <- c("divergent__", "treedepth__")
        if (nrow(params) != m$iter || !all(required %in% colnames(params))) sampler_known <- FALSE
        if (all(required %in% colnames(params)) &&
            any(!is.finite(params[, required, drop = FALSE]))) sampler_known <- FALSE
        if ("divergent__" %in% colnames(params)) divergent <- divergent + sum(params[, "divergent__"])
        if ("treedepth__" %in% colnames(params) && length(m$max_treedepth) == 1L && is.finite(m$max_treedepth)) {
          depth_hits <- depth_hits + sum(params[, "treedepth__"] >= m$max_treedepth)
        }
      }
      metrics <- list(parameters = as.data.frame(diagnostics), divergences = divergent,
                      treedepth_hits = depth_hits, sampler_diagnostics_known = sampler_known)
      metrics$biology <- .opal_check_draws(x, penalty_tolerance, catch_tolerance)
      passes <- m$chains >= 2L && all(is.finite(diagnostics$rhat)) &&
        all(diagnostics$rhat <= max_rhat) && all(is.finite(diagnostics$ess_bulk)) &&
        all(is.finite(diagnostics$ess_tail)) && all(diagnostics$ess_bulk >= min_ess) &&
        all(diagnostics$ess_tail >= min_ess) && sampler_known && divergent == 0 && depth_hits == 0 &&
        metrics$biology$passes
    }
    list(passes = isTRUE(passes), metrics = metrics, error = NULL)
  }, error = function(e) list(passes = FALSE, metrics = list(), error = conditionMessage(e)))
  record <- c(result, list(identity = .opal_check_identity(x, scope),
    version = .opal_validation_version,
    settings = list(gradient_tolerance = gradient_tolerance, max_rhat = max_rhat, min_ess = min_ess,
                    penalty_tolerance = penalty_tolerance, catch_tolerance = catch_tolerance),
    checked_at = format(Sys.time(), tz = "UTC", usetz = TRUE)))
  record$payload_id <- .opal_obj_hash(record)
  x$validation[[scope]] <- record
  x <- .opal_obj_seal(x)
  if (!record$passes) {
    message <- paste0("The ", scope, " check failed; results are retained.",
                      if (!is.null(record$error)) paste0(" ", record$error))
    if (stop_on_failure) stop(message, call. = FALSE)
    warning(message, call. = FALSE)
  }
  x
}
