# Portable staged assessment objects -----------------------------------------

.opal_obj_schema <- 1L
.opal_obj_cache <- new.env(parent = emptyenv())

.opal_obj_hash <- function(x) {
  # Version 2 avoids ALTREP encodings. The writer-version header is not part
  # of the payload identity: the same data must hash equally across R versions.
  bytes <- serialize(x, NULL, version = 2L)
  bytes[7:10] <- as.raw(0L)
  path <- tempfile()
  on.exit(unlink(path), add = TRUE)
  writeBin(bytes, path)
  unname(tools::md5sum(path))
}

.opal_obj_identity <- function(x) {
  target <- list(
    model = x$model[c("name", "schema_version", "scientific_version", "signature")],
    data = x$data, parameters = x$parameters, map = x$map,
    random = x$random, makeadfun_args = x$makeadfun_args, bounds = x$bounds
  )
  list(target = .opal_obj_hash(target), fit = .opal_obj_hash(x$fit[c(
    "opt", "parameters", "last_par_best"
  )]))
}

.opal_obj_runtime_id <- function(x) .opal_obj_hash(x$identity)

.opal_obj_stage <- function(x) {
  if (!is.null(x$mcmc)) return("sampled")
  if (!is.null(x$fit$opt)) return("fitted")
  if (!is.null(x$build)) return("built")
  if (!is.null(x$parameters) && !is.null(x$map)) return("configured")
  "data"
}

.opal_obj_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop("`", name, "` must be TRUE or FALSE.", call. = FALSE)
  }
}

.opal_obj_portable <- function(x, name) {
  if (.opal_contains_nonportable(x)) {
    stop("`", name, "` contains non-portable state.", call. = FALSE)
  }
}

.opal_obj_parameters <- function(parameters, map, random) {
  if (!is.character(random) || anyNA(random) || any(!nzchar(random)) || anyDuplicated(random)) {
    stop("Invalid random-effect names.", call. = FALSE)
  }
  if (is.null(parameters)) {
    if (!is.null(map)) {
      .opal_validate_named_list(map, "map")
      if (any(!vapply(map, is.factor, logical(1L)))) stop("Map entries must be factors.", call. = FALSE)
    }
    # Parameter names and dimensions can only be checked after defaults resolve.
    return(invisible(NULL))
  }
  .opal_validate_named_list(parameters, "parameters")
  if (!length(parameters) || any(vapply(parameters, function(p) {
    !is.numeric(p) || !length(p) || any(!is.finite(p))
  }, logical(1L)))) stop("Parameters must be finite numeric arrays.", call. = FALSE)
  if (!is.null(map)) {
    .opal_validate_named_list(map, "map")
    if (length(setdiff(names(map), names(parameters)))) {
      stop("The map contains unknown parameters.", call. = FALSE)
    }
    for (name in names(map)) {
      if (!is.factor(map[[name]]) ||
          length(map[[name]]) != length(parameters[[name]]) ||
          (!is.null(dim(map[[name]])) &&
           !identical(dim(map[[name]]), dim(parameters[[name]])))) {
        stop("Invalid map for `", name, "`.", call. = FALSE)
      }
    }
  }
  if (!is.character(random) || anyNA(random) || anyDuplicated(random) ||
      length(setdiff(random, names(parameters)))) {
    stop("Invalid random-effect names.", call. = FALSE)
  }
}

#' Validate a staged Opal object
#'
#' Checks structural and model-state integrity without building an objective.
#' Full result verification additionally checks posterior payload identities.
#' @param x An `opal_obj`.
#' @param results Verify stored posterior draws and attempt history.
#' @return `x`, invisibly.
#' @export
validate_opal_obj <- function(x, results = FALSE) {
  .opal_obj_flag(results, "results")
  required <- c("schema_version", "model", "data", "parameters", "map",
                "random", "bounds", "control", "makeadfun_args", "configuration",
                "build", "fit", "mcmc", "mcmc_history", "derived", "validation",
                "provenance", "identity")
  if (!inherits(x, "opal_obj") || !is.list(x) ||
      length(setdiff(required, names(x))) ||
      !identical(x$schema_version, .opal_obj_schema)) {
    stop("Invalid or unsupported opal_obj schema.", call. = FALSE)
  }
  for (name in c("data", "makeadfun_args", "configuration", "fit", "derived",
                 "validation", "provenance")) {
    .opal_validate_named_list(x[[name]], name)
  }
  .opal_obj_parameters(x$parameters, x$map, x$random)
  if (!is.null(x$control)) .opal_validate_named_list(x$control, "control")
  if (length(intersect(names(x$makeadfun_args),
                       c("func", "parameters", "map", "random", "silent")))) {
    stop("Core MakeADFun arguments cannot be overridden.", call. = FALSE)
  }
  if (!identical(x$model$name, "opal_model") ||
      !is.character(x$model$scientific_version)) {
    stop("Invalid model metadata.", call. = FALSE)
  }
  if (!is.null(x$fit$opt)) {
    opt <- x$fit$opt
    if (!is.list(opt) || !is.numeric(opt$par) || is.null(names(opt$par)) ||
        any(!is.finite(opt$par)) || length(opt$objective) != 1L ||
        !is.finite(opt$objective) || is.null(x$fit$parameters) ||
        !is.numeric(x$fit$last_par_best) || any(!is.finite(x$fit$last_par_best))) {
      stop("Invalid fitted state.", call. = FALSE)
    }
    .opal_obj_parameters(x$fit$parameters, x$map, x$random)
    if (!identical(names(x$fit$parameters), names(x$parameters)) ||
        !identical(lapply(x$fit$parameters, dim), lapply(x$parameters, dim)) ||
        !identical(lengths(x$fit$parameters), lengths(x$parameters))) {
      stop("Fitted parameter structure differs from the configuration.", call. = FALSE)
    }
    if (!is.null(x$bounds)) .opal_normalize_bounds(x$bounds, opt$par)
  }
  if (!identical(x$identity, .opal_obj_identity(x))) {
    stop("Model state was modified directly; use opal_update().", call. = FALSE)
  }
  if (!is.list(x$mcmc_history)) stop("Invalid MCMC history.", call. = FALSE)
  if (!is.null(x$mcmc) &&
      (!inherits(x$mcmc, "opal_mcmc") ||
       !identical(x$mcmc$target_id, x$identity$target))) {
    stop("Posterior does not belong to this model configuration.", call. = FALSE)
  }
  if (results) {
    .opal_obj_portable(unclass(x), "opal_obj")
    payloads <- c(list(x$mcmc), lapply(x$mcmc_history, `[[`, "mcmc"))
    for (m in payloads) {
      if (is.null(m)) next
      .opal_portable_mcmc(m, m$settings)
      if (!identical(m$payload_id, .opal_obj_hash(m[setdiff(names(m), "payload_id")]))) {
        stop("Stored posterior payload was modified.", call. = FALSE)
      }
    }
  }
  invisible(x)
}

.opal_obj_seal <- function(x) {
  x$identity <- .opal_obj_identity(x)
  x$provenance$updated_at <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  validate_opal_obj(x)
  x
}

.opal_obj_compatible <- function(x) {
  current <- .opal_model_metadata()
  if (!identical(x$model$scientific_version, current$scientific_version) ||
      !identical(x$model$schema_version, current$schema_version)) {
    stop("This object uses a different scientific model contract.", call. = FALSE)
  }
  if (!identical(x$model$signature, current$signature)) {
    warning("The model implementation checksum differs from the saved object.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Create a portable Opal assessment object
#'
#' Creates an S3 object before fitting. Configuration, fitted results, and
#' posterior samples remain ordinary R data. RTMB objectives are cached only
#' within the current session. Explicit inputs are recommended for custom
#' stocks; `model` identifies a bundled parameter configuration when needed.
#' @param data Named model data list.
#' @param parameters Named initial parameter list, or NULL for bundled defaults.
#' @param map RTMB map. NULL resolves defaults; `list()` explicitly frees all parameters.
#' @param random Names of random-effect parameters.
#' @param bounds Bounds accepted by [get_bounds()], or NULL to derive them.
#' @param control Optimiser controls.
#' @param model Optional bundled configuration name for [get_parameters()].
#' @param makeadfun_args Additional portable arguments for `RTMB::MakeADFun()`.
#' @param metadata Named user metadata.
#' @return An `opal_obj`. Construction never optimises or samples.
#' @export
#' @examples
#' inputs <- opaka_quickstart_inputs()
#' x <- opal_obj(inputs$data, inputs$parameters, inputs$map)
#' summary(x)
opal_obj <- function(data, parameters = NULL, map = NULL, random = character(),
                     bounds = NULL, control = NULL, model = NULL,
                     makeadfun_args = list(), metadata = list()) {
  data <- .opal_sanitize_tables(data)
  .opal_validate_named_list(metadata, "metadata")
  x <- structure(list(
    schema_version = .opal_obj_schema, model = .opal_model_metadata(),
    data = data, parameters = parameters, map = map, random = random,
    bounds = bounds, control = control, makeadfun_args = makeadfun_args,
    configuration = list(model = model, origins = list(
      parameters = if (is.null(parameters)) "unresolved" else "supplied",
      map = if (is.null(map)) "unresolved" else "supplied",
      bounds = if (is.null(bounds)) "unresolved" else "supplied")),
    build = NULL, fit = list(opt = NULL, parameters = NULL, last_par_best = NULL,
                            diagnostics = list(), estimability = NULL),
    mcmc = NULL, mcmc_history = list(), derived = list(), validation = list(),
    provenance = list(created_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
                      metadata = metadata), identity = NULL
  ), class = c("opal_obj", "list"))
  .opal_obj_portable(unclass(x), "opal_obj")
  .opal_obj_seal(x)
}

.opal_obj_make_runtime <- function(x, silent = TRUE) {
  .opal_obj_compatible(x)
  if (is.null(x$parameters) || is.null(x$map)) {
    stop("Build the configuration with opal_build() first.", call. = FALSE)
  }
  parameters <- .opal_fit_or(x$fit$parameters, x$parameters)
  object <- do.call(RTMB::MakeADFun, c(list(
    func = cmb(opal_model, x$data), parameters = parameters, map = x$map,
    random = x$random, silent = silent), x$makeadfun_args))
  if (!is.null(x$fit$opt)) {
    if (!identical(names(object$par), names(x$fit$opt$par))) {
      stop("Rebuilt active parameter layout differs.", call. = FALSE)
    }
    object$env$value.best <- Inf
    value <- object$fn(x$fit$opt$par)
    if (!isTRUE(all.equal(as.numeric(value), as.numeric(x$fit$opt$objective),
                         tolerance = 1e-6))) {
      stop("Rebuilt objective differs from the saved fitted objective.", call. = FALSE)
    }
    if (!identical(names(object$env$last.par.best), names(x$fit$last_par_best))) {
      stop("Rebuilt complete parameter layout differs.", call. = FALSE)
    }
    if (!isTRUE(all.equal(as.numeric(object$env$last.par.best),
                          as.numeric(x$fit$last_par_best), tolerance = 1e-6))) {
      stop("Rebuilt complete fitted state differs.", call. = FALSE)
    }
    object$par <- x$fit$opt$par
    object$opt <- x$fit$opt
    object$env$last.par.best <- x$fit$last_par_best
  } else {
    if (!is.null(x$build) && !identical(x$build$active_names,
          .opal_expand_parameter_names(names(object$par)))) {
      stop("Rebuilt active parameter layout differs.", call. = FALSE)
    }
    value <- object$fn(object$par)
    if (!is.finite(value) || any(!is.finite(object$gr(object$par)))) {
      stop("Initial objective or gradient is not finite.", call. = FALSE)
    }
    if (!is.null(x$build) &&
        !isTRUE(all.equal(as.numeric(value), x$build$objective, tolerance = 1e-6))) {
      stop("Rebuilt initial objective differs from the saved build.", call. = FALSE)
    }
  }
  object
}

#' Build or access an Opal runtime objective
#' @param x An `opal_obj`.
#' @param silent Silence RTMB construction messages.
#' @return `opal_build()` returns the updated object; `opal_rtmb()` returns a
#'   transient RTMB objective. Use `fresh = TRUE` for isolated mutable work.
#' @export
opal_build <- function(x, silent = TRUE) {
  validate_opal_obj(x)
  .opal_obj_flag(silent, "silent")
  if (!is.null(x$build)) {
    opal_rtmb(x)
    return(x)
  }
  if (is.null(x$parameters)) {
    if (is.null(x$configuration$model)) {
      stop("Supply parameters or an explicit bundled `model` configuration.", call. = FALSE)
    }
    x$parameters <- get_parameters(model = x$configuration$model)
    x$configuration$origins$parameters <- "default"
  }
  if (is.null(x$map)) {
    x$map <- get_map(x$parameters)
    x$configuration$origins$map <- "default"
  }
  .opal_obj_parameters(x$parameters, x$map, x$random)
  object <- .opal_obj_make_runtime(x, silent)
  if (is.null(x$bounds)) {
    x$bounds <- get_bounds(object, x$parameters)
    x$configuration$origins$bounds <- "default"
  }
  x$bounds <- .opal_normalize_bounds(x$bounds, object$par)
  if (any(object$par < x$bounds$lower | object$par > x$bounds$upper)) {
    stop("Initial parameters are outside the bounds.", call. = FALSE)
  }
  x$build <- list(active_names = .opal_expand_parameter_names(names(object$par)),
                  objective = as.numeric(object$fn(object$par)))
  x <- .opal_obj_seal(x)
  assign(.opal_obj_runtime_id(x), object, envir = .opal_obj_cache)
  x
}

#' @rdname opal_build
#' @param fresh Build an isolated objective instead of accessing the cache.
#' @export
opal_rtmb <- function(x, fresh = FALSE) {
  validate_opal_obj(x)
  .opal_obj_flag(fresh, "fresh")
  .opal_obj_compatible(x)
  id <- .opal_obj_runtime_id(x)
  if (!fresh && exists(id, .opal_obj_cache, inherits = FALSE)) {
    return(get(id, .opal_obj_cache, inherits = FALSE))
  }
  object <- .opal_obj_make_runtime(x)
  if (!fresh) assign(id, object, envir = .opal_obj_cache)
  object
}

#' Report a configured or fitted Opal model
#' @param x An `opal_obj`.
#' @return The named model report at the stored parameters.
#' @export
opal_report <- function(x) {
  object <- opal_rtmb(x, fresh = TRUE)
  object$report(object$env$last.par.best)
}

#' Summarise a staged Opal assessment
#' @param x,object An `opal_obj` or its summary.
#' @param ... Unused arguments.
#' @return `summary()` returns a compact list; print methods return invisibly.
#' @export
summary.opal_obj <- function(object, ...) {
  validate_opal_obj(object)
  checks <- lapply(c("fit", "mcmc"), function(scope) {
    record <- object$validation[[scope]]
    if (is.null(record)) return("not run")
    if (!identical(record$identity, .opal_check_identity(object, scope))) return("stale")
    if (isTRUE(record$passes)) "passed" else "failed"
  })
  structure(list(stage = .opal_obj_stage(object),
    scientific_version = object$model$scientific_version,
    parameters = if (is.null(object$build)) NA_integer_ else length(object$build$active_names),
    objective = object$fit$opt$objective, checks = stats::setNames(checks, c("fit", "mcmc")),
    mcmc = if (is.null(object$mcmc)) NULL else list(chains = object$mcmc$chains,
      retained = object$mcmc$iter,
      predates_fit = !identical(object$mcmc$fit_id, object$identity$fit)),
    attempts = length(object$mcmc_history), metadata = object$provenance$metadata),
    class = c("summary.opal_obj", "list"))
}

#' @rdname summary.opal_obj
#' @export
print.summary.opal_obj <- function(x, ...) {
  cat("<opal_obj> ", x$stage, "\n", sep = "")
  cat("Active parameters:", x$parameters, "\n")
  if (!is.null(x$objective)) cat("Objective:", format(x$objective), "\n")
  cat("Fit check:", x$checks$fit, " | MCMC check:", x$checks$mcmc, "\n")
  if (!is.null(x$mcmc)) {
    cat("MCMC:", x$mcmc$chains, "chains;", x$mcmc$retained, "retained iterations\n")
    if (x$mcmc$predates_fit) cat("MCMC predates the current optimisation.\n")
  }
  invisible(x)
}

#' @rdname summary.opal_obj
#' @export
print.opal_obj <- function(x, ...) {
  print(summary(x), ...)
  invisible(x)
}
