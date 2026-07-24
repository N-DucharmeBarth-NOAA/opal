# Portable fitted-model objects ------------------------------------------------

.opal_fit_schema_version <- 1L
.opal_model_schema_version <- 1L
.opal_model_scientific_version <- "opal_model_contract_v1"
.opal_fit_runtime_cache <- new.env(parent = emptyenv())

.opal_fit_or <- function(x, y) {
  if (is.null(x)) y else x
}

.opal_expand_parameter_names <- function(x) {
  if (is.null(x)) return(character())
  index <- ave(seq_along(x), x, FUN = seq_along)
  count <- ave(seq_along(x), x, FUN = length)
  ifelse(count == 1L, x, paste0(x, "[", index, "]"))
}

.opal_contains_nonportable <- function(x) {
  if (typeof(x) %in% c("closure", "environment", "externalptr", "weakref")) {
    return(TRUE)
  }
  if (is.list(x) && any(vapply(x, .opal_contains_nonportable, logical(1L)))) {
    return(TRUE)
  }
  x_attributes <- attributes(x)
  if (!is.null(x_attributes)) {
    return(any(vapply(
      x_attributes,
      .opal_contains_nonportable,
      logical(1L)
    )))
  }
  FALSE
}

.opal_sanitize_tables <- function(x) {
  if (inherits(x, "spec_tbl_df")) {
    attr(x, "spec") <- NULL
    attr(x, "problems") <- NULL
    class(x) <- setdiff(class(x), "spec_tbl_df")
  }
  if (is.list(x)) {
    for (i in seq_along(x)) {
      x[i] <- list(.opal_sanitize_tables(x[[i]]))
    }
  }
  x
}

.opal_validate_named_list <- function(x, label) {
  if (!is.list(x)) {
    stop("`", label, "` must be a list.", call. = FALSE)
  }
  if (length(x) &&
      (is.null(names(x)) || any(!nzchar(names(x))) || anyDuplicated(names(x)))) {
    stop(
      "`", label, "` must have non-empty, unique names.",
      call. = FALSE
    )
  }
  invisible(x)
}

.opal_package_version <- function(package) {
  tryCatch(
    as.character(utils::packageVersion(package)),
    error = function(e) NA_character_
  )
}

.opal_function_text <- function(fun) {
  paste(
    c(
      deparse(formals(fun), width.cutoff = 500L),
      deparse(body(fun), width.cutoff = 500L)
    ),
    collapse = "\n"
  )
}

.opal_text_checksum <- function(text) {
  file <- tempfile(fileext = ".txt")
  on.exit(unlink(file), add = TRUE)
  writeLines(enc2utf8(paste(text, collapse = "\n")), file, useBytes = TRUE)
  unname(tools::md5sum(file))
}

.opal_object_checksum <- function(x) {
  file <- tempfile(fileext = ".rds")
  on.exit(unlink(file), add = TRUE)
  saveRDS(x, file, compress = FALSE, version = 3L)
  unname(tools::md5sum(file))
}

.opal_model_metadata <- function() {
  dependencies <- c("RTMB", "RTMBdist", "SparseNUTS")
  list(
    name = "opal_model",
    schema_version = .opal_model_schema_version,
    scientific_version = .opal_model_scientific_version,
    signature = .opal_text_checksum(.opal_function_text(opal_model)),
    package = "opal",
    package_version = .opal_package_version("opal"),
    dependency_versions = stats::setNames(
      vapply(dependencies, .opal_package_version, character(1L)),
      dependencies
    ),
    r_version = as.character(getRversion())
  )
}

.opal_normalize_bounds <- function(bounds, opt_par) {
  if (is.null(bounds)) return(NULL)

  if (is.data.frame(bounds)) {
    if (!all(c("lower", "upper") %in% names(bounds))) {
      stop("A bounds data frame must contain `lower` and `upper`.", call. = FALSE)
    }
    supplied_names <- bounds$active_name
    if (is.null(supplied_names)) supplied_names <- bounds$parameter
    bounds <- list(
      lower = bounds$lower,
      upper = bounds$upper,
      parameter = supplied_names
    )
  }
  if (!is.list(bounds) || !all(c("lower", "upper") %in% names(bounds))) {
    stop(
      "`bounds` must be a data frame or list containing `lower` and `upper`.",
      call. = FALSE
    )
  }

  lower <- as.numeric(bounds$lower)
  upper <- as.numeric(bounds$upper)
  if (length(lower) != length(opt_par) ||
      length(upper) != length(opt_par)) {
    stop("Stored bounds do not match `opt$par`.", call. = FALSE)
  }
  if (anyNA(lower) || anyNA(upper) || any(lower > upper)) {
    stop("Stored bounds are invalid.", call. = FALSE)
  }

  raw_names <- names(opt_par)
  active_names <- .opal_expand_parameter_names(raw_names)
  supplied_names <- bounds$active_name
  if (is.null(supplied_names)) supplied_names <- bounds$parameter
  if (is.null(supplied_names) && !is.null(names(bounds$lower))) {
    supplied_names <- names(bounds$lower)
  }
  if (!is.null(supplied_names) &&
      !identical(as.character(supplied_names), raw_names) &&
      !identical(as.character(supplied_names), active_names)) {
    stop("Stored bounds are not ordered like `opt$par`.", call. = FALSE)
  }

  list(
    lower = stats::setNames(lower, active_names),
    upper = stats::setNames(upper, active_names)
  )
}

.opal_compact_estimability <- function(x, active_names, rtmb_names) {
  if (is.null(x)) return(NULL)

  make_summary <- function(status, message, n_bad = NA_integer_,
                           implicated = character()) {
    list(
      status = status,
      message = as.character(message),
      n_parameters = as.integer(length(active_names)),
      n_non_estimable_combinations = as.integer(n_bad),
      implicated_parameters = unique(as.character(implicated))
    )
  }

  if (inherits(x, "condition")) {
    return(make_summary("error", conditionMessage(x)))
  }
  if (!is.list(x)) {
    stop(
      "`estimability` must be a check result, compact summary, condition, or NULL.",
      call. = FALSE
    )
  }
  if (!is.null(x$status)) {
    if (!x$status %in% c("estimable", "non_estimable", "error")) {
      stop("Unknown compact estimability status.", call. = FALSE)
    }
    implicated <- .opal_fit_or(x$implicated_parameters, character())
    unknown <- setdiff(implicated, active_names)
    if (length(unknown)) {
      stop("Estimability output contains unknown parameters.", call. = FALSE)
    }
    return(make_summary(
      x$status,
      .opal_fit_or(x$message, x$status),
      .opal_fit_or(x$n_non_estimable_combinations, NA_integer_),
      implicated
    ))
  }
  if (!"WhichBad" %in% names(x)) {
    stop("Full estimability results must contain `WhichBad`.", call. = FALSE)
  }

  which_bad <- as.integer(x$WhichBad)
  implicated <- character()
  if (is.data.frame(x$BadParams) &&
      all(c("Param", "Param_check") %in% names(x$BadParams)) &&
      identical(as.character(x$BadParams$Param), rtmb_names)) {
    bad <- is.na(x$BadParams$Param_check) |
      as.character(x$BadParams$Param_check) != "OK"
    implicated <- active_names[which(bad)]
  }
  if (!length(which_bad)) {
    return(make_summary(
      "estimable",
      paste(
        "All", format(length(active_names), big.mark = ","),
        "active fixed-effect parameters are estimable."
      ),
      0L
    ))
  }
  make_summary(
    "non_estimable",
    paste(
      format(length(which_bad), big.mark = ","),
      "non-estimable parameter combination(s) detected."
    ),
    length(which_bad),
    implicated
  )
}

.opal_mcmc_sample_names <- function(samples) {
  sample_dimnames <- dimnames(samples)
  if (length(sample_dimnames) >= 3L) sample_dimnames[[3L]] else NULL
}

.opal_normalize_timing <- function(mcmc, chains) {
  timing <- list()
  known <- TRUE
  for (name in c("time.warmup", "time.sampling", "time.total")) {
    value <- mcmc[[name]]
    if (is.null(value) && is.list(mcmc$timing)) {
      value <- mcmc$timing[[name]]
    }
    if (is.null(value)) {
      value <- rep(0, chains)
      known <- FALSE
    }
    value <- as.numeric(value)
    if (length(value) == 1L && chains > 1L) value <- rep(value, chains)
    if (length(value) != chains || any(!is.finite(value)) || any(value < 0)) {
      value <- rep(0, chains)
      known <- FALSE
    }
    timing[[name]] <- value
  }
  list(values = timing, known = known)
}

.opal_portable_or_null <- function(x) {
  if (is.null(x) || .opal_contains_nonportable(x)) NULL else x
}

.opal_portable_mcmc <- function(mcmc, settings = list()) {
  if (is.null(mcmc)) return(NULL)
  .opal_validate_named_list(settings, "mcmc_settings")

  if (is.matrix(mcmc)) {
    sample_names <- colnames(mcmc)
    samples <- array(
      mcmc,
      dim = c(nrow(mcmc), 1L, ncol(mcmc)),
      dimnames = list(
        iteration = rownames(mcmc),
        chain = "1",
        variable = sample_names
      )
    )
    mcmc <- list(samples = samples, warmup = 0L, iter = nrow(samples))
  } else if (is.array(mcmc) && length(dim(mcmc)) == 3L) {
    mcmc <- list(samples = mcmc)
  }
  if (!is.list(mcmc) || is.null(mcmc$samples)) {
    stop(
      "`mcmc` must be a SparseNUTS-style fit, matrix, or three-dimensional array.",
      call. = FALSE
    )
  }

  samples <- as.array(mcmc$samples)
  if (length(dim(samples)) != 3L || any(dim(samples) < 1L) ||
      !is.numeric(samples) || any(!is.finite(samples))) {
    stop(
      "Stored MCMC samples must be a finite numeric iteration-chain-variable array.",
      call. = FALSE
    )
  }
  sample_names <- .opal_mcmc_sample_names(samples)
  if (is.null(sample_names) || any(!nzchar(sample_names)) ||
      anyDuplicated(sample_names)) {
    stop("Stored MCMC variables must have non-empty, unique names.", call. = FALSE)
  }
  par_names <- .opal_fit_or(mcmc$par_names, setdiff(sample_names, "lp__"))
  par_names <- as.character(par_names)
  if (!identical(par_names, setdiff(sample_names, "lp__"))) {
    stop("MCMC parameter names do not match the sample array.", call. = FALSE)
  }

  warmup <- as.integer(.opal_fit_or(mcmc$warmup, 0L))
  iter <- as.integer(.opal_fit_or(mcmc$iter, dim(samples)[1L] - warmup))
  thin <- as.integer(.opal_fit_or(mcmc$thin, 1L))
  chains <- as.integer(dim(samples)[2L])
  if (length(warmup) != 1L || is.na(warmup) || warmup < 0L ||
      length(iter) != 1L || is.na(iter) || iter < 1L ||
      warmup + iter != dim(samples)[1L] ||
      length(thin) != 1L || is.na(thin) || thin < 1L) {
    stop("MCMC warmup, iteration, or thinning metadata is invalid.", call. = FALSE)
  }

  timing <- .opal_normalize_timing(mcmc, chains)
  monitor <- if (is.null(mcmc$monitor)) {
    NULL
  } else {
    as.data.frame(mcmc$monitor)
  }
  mle <- mcmc$mle
  if (is.list(mle)) {
    mle <- mle[intersect(c("est", "se", "cor"), names(mle))]
  }

  out <- list(
    samples = samples,
    sampler_params = .opal_portable_or_null(mcmc$sampler_params),
    par_names = par_names,
    sample_names = sample_names,
    monitor = monitor,
    warmup = warmup,
    iter = iter,
    thin = thin,
    chains = chains,
    algorithm = as.character(.opal_fit_or(mcmc$algorithm, "stored")),
    metric = as.character(.opal_fit_or(mcmc$metric, NA_character_)),
    max_treedepth = suppressWarnings(as.integer(.opal_fit_or(
      mcmc$max_treedepth,
      NA_integer_
    ))),
    model = as.character(.opal_fit_or(mcmc$model, "RTMB")),
    mle = .opal_portable_or_null(mle),
    timing = timing$values,
    timing_known = timing$known,
    inits = .opal_portable_or_null(mcmc$inits),
    samples_unbounded = .opal_portable_or_null(mcmc$samples_unbounded),
    settings = utils::modifyList(
      list(timing_known = timing$known),
      settings
    ),
    source_class = class(mcmc),
    parameter_scope = NULL
  )
  class(out) <- c("opal_mcmc", "list")
  out
}

.opal_attach_mcmc <- function(mcmc, obj, settings = list()) {
  portable <- .opal_portable_mcmc(mcmc, settings = settings)
  if (is.null(portable)) return(NULL)

  active_names <- .opal_expand_parameter_names(names(obj$par))
  complete_names <- .opal_expand_parameter_names(names(obj$env$last.par.best))
  portable$parameter_scope <- if (identical(portable$par_names, active_names)) {
    "active"
  } else if (identical(portable$par_names, complete_names)) {
    "complete"
  } else {
    stop(
      "Stored MCMC parameter names do not match the active or complete model layout.",
      call. = FALSE
    )
  }
  portable
}

.opal_fit_runtime_payload <- function(x) {
  x[setdiff(names(x), "runtime_id")]
}

.opal_fit_runtime_id <- function(x) {
  paste0("opal-fit-", .opal_object_checksum(.opal_fit_runtime_payload(x)))
}

.opal_cache_runtime_object <- function(x, obj) {
  assign(x$runtime_id, obj, envir = .opal_fit_runtime_cache)
  invisible(x)
}

.opal_cached_runtime_object <- function(x) {
  if (exists(x$runtime_id, envir = .opal_fit_runtime_cache, inherits = FALSE)) {
    return(get(x$runtime_id, envir = .opal_fit_runtime_cache, inherits = FALSE))
  }
  NULL
}

#' Create a portable fitted opal model object
#'
#' Captures the plain R state needed to reproduce a fitted opal model without
#' serializing the transient RTMB objective. The objective is retained only in
#' a session cache and can be rebuilt with [rebuild_opal_object()]. Optional
#' SparseNUTS output is normalized to a package-owned `opal_mcmc` payload.
#'
#' @param data Named model data list used to construct `obj`.
#' @param obj Fitted RTMB objective created with `opal_model`.
#' @param opt Optimizer result, normally returned by [stats::nlminb()].
#' @param bounds Optional bounds data frame from [get_bounds()] or a list with
#'   `lower` and `upper`.
#' @param control Optional optimizer control list.
#' @param estimability Optional output from [check_estimability()]. A compact
#'   summary is retained.
#' @param diagnostics Optional named list of fit diagnostics.
#' @param metadata Optional named list of user metadata.
#' @param mcmc Optional SparseNUTS-style fit, posterior matrix, or
#'   iteration-chain-variable array.
#' @param mcmc_settings Optional named list of sampler settings not already
#'   present in `mcmc`.
#' @param derived Optional named list of portable derived results, such as
#'   projections or retrospective summaries.
#' @param makeadfun_args Optional named list of additional arguments needed to
#'   rebuild `RTMB::MakeADFun()`. Core arguments are reserved.
#' @param optimizer Non-empty character name of the optimizer used.
#'
#' @return An object inheriting from `opal_fit`.
#' @export
#'
opal_fit <- function(data, obj, opt, bounds = NULL, control = NULL,
                     estimability = NULL, diagnostics = list(),
                     metadata = list(), mcmc = NULL, mcmc_settings = list(),
                     derived = list(), makeadfun_args = list(),
                     optimizer = "nlminb") {
  if (!is.list(data)) stop("`data` must be a list.", call. = FALSE)
  data <- .opal_sanitize_tables(data)
  .opal_validate_named_list(data, "data")
  .opal_validate_named_list(diagnostics, "diagnostics")
  .opal_validate_named_list(metadata, "metadata")
  .opal_validate_named_list(mcmc_settings, "mcmc_settings")
  .opal_validate_named_list(derived, "derived")
  .opal_validate_named_list(makeadfun_args, "makeadfun_args")

  if (!is.list(obj) || is.null(obj$env) || !is.function(obj$fn) ||
      !is.function(obj$report)) {
    stop(
      "`obj` must be a fitted RTMB objective with `env`, `fn`, and `report`.",
      call. = FALSE
    )
  }
  if (!is.list(opt) || !is.numeric(opt$par) || is.null(names(opt$par)) ||
      any(!is.finite(opt$par))) {
    stop("`opt$par` must be a named, finite numeric vector.", call. = FALSE)
  }
  if (!is.numeric(opt$objective) || length(opt$objective) != 1L ||
      !is.finite(opt$objective)) {
    stop("`opt$objective` must be one finite numeric value.", call. = FALSE)
  }
  if (!is.numeric(opt$convergence) || length(opt$convergence) != 1L ||
      !is.finite(opt$convergence) ||
      opt$convergence != as.integer(opt$convergence)) {
    stop("`opt$convergence` must be one integer-like value.", call. = FALSE)
  }
  opt$convergence <- as.integer(opt$convergence)
  if (!is.null(control) && !is.list(control)) {
    stop("`control` must be a list or `NULL`.", call. = FALSE)
  }
  if (!is.character(optimizer) || length(optimizer) != 1L ||
      is.na(optimizer) || !nzchar(optimizer)) {
    stop("`optimizer` must be one non-empty character value.", call. = FALSE)
  }

  reserved <- intersect(
    names(makeadfun_args),
    c("func", "parameters", "map", "random", "silent")
  )
  if (length(reserved)) {
    stop(
      "`makeadfun_args` must not contain: ",
      paste(reserved, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  map <- .opal_fit_or(obj$env$map, list())
  random <- as.character(.opal_fit_or(obj$env$.random, character()))
  .opal_validate_named_list(map, "obj$env$map")
  if (any(!nzchar(random)) || anyDuplicated(random)) {
    stop("The random-effect specification in `obj` is invalid.", call. = FALSE)
  }

  active_names <- .opal_expand_parameter_names(names(obj$par))
  opt_names <- .opal_expand_parameter_names(names(opt$par))
  if (!identical(active_names, opt_names)) {
    stop("`opt$par` does not match the active parameter layout in `obj`.",
         call. = FALSE)
  }

  objective <- as.numeric(obj$fn(opt$par))
  if (length(objective) != 1L || !is.finite(objective) ||
      !isTRUE(all.equal(objective, as.numeric(opt$objective), tolerance = 1e-6))) {
    stop("`obj` does not reproduce `opt$objective` at `opt$par`.",
         call. = FALSE)
  }
  last_par_best <- obj$env$last.par.best
  if (!is.numeric(last_par_best) || is.null(names(last_par_best)) ||
      any(!is.finite(last_par_best))) {
    stop("`obj` did not retain a complete finite best parameter vector.",
         call. = FALSE)
  }
  parameters <- obj$env$parList(opt$par)
  .opal_validate_named_list(parameters, "fitted parameters")

  unknown_map <- setdiff(names(map), names(parameters))
  if (length(unknown_map)) {
    stop("The parameter map contains unknown parameters.", call. = FALSE)
  }
  for (name in names(map)) {
    if (!is.factor(map[[name]]) ||
        length(map[[name]]) != length(parameters[[name]]) ||
        (!is.null(dim(map[[name]])) &&
         !identical(
           as.integer(dim(map[[name]])),
           as.integer(dim(parameters[[name]]))
         ))) {
      stop("The map entry for `", name, "` is invalid.", call. = FALSE)
    }
  }
  if (length(setdiff(random, names(parameters)))) {
    stop("The random-effect specification contains unknown parameters.",
         call. = FALSE)
  }

  bounds <- .opal_normalize_bounds(bounds, opt$par)
  if (!is.null(bounds) &&
      any(opt$par < bounds$lower - 1e-8 | opt$par > bounds$upper + 1e-8)) {
    stop("`opt$par` is outside the supplied bounds.", call. = FALSE)
  }

  portable_inputs <- list(
    data = data,
    parameters = parameters,
    map = map,
    makeadfun_args = makeadfun_args,
    diagnostics = diagnostics,
    metadata = metadata,
    mcmc_settings = mcmc_settings,
    derived = derived
  )
  nonportable <- names(portable_inputs)[vapply(
    portable_inputs,
    .opal_contains_nonportable,
    logical(1L)
  )]
  if (length(nonportable)) {
    stop(
      "Portable fit state contains non-portable values: ",
      paste(nonportable, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  portable_mcmc <- .opal_attach_mcmc(
    mcmc,
    obj = obj,
    settings = mcmc_settings
  )
  out <- structure(
    list(
      schema_version = .opal_fit_schema_version,
      model = .opal_model_metadata(),
      data = data,
      parameters = parameters,
      map = map,
      random = random,
      makeadfun_args = makeadfun_args,
      bounds = bounds,
      control = control,
      optimization = list(method = optimizer),
      fit = list(
        opt = opt,
        estimability = .opal_compact_estimability(
          estimability,
          active_names = active_names,
          rtmb_names = names(opt$par)
        ),
        diagnostics = diagnostics
      ),
      mcmc = portable_mcmc,
      derived = derived,
      provenance = list(
        created_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
        updated_at = NULL,
        metadata = metadata
      ),
      runtime_id = NULL
    ),
    class = c("opal_fit", "list")
  )
  out$runtime_id <- .opal_fit_runtime_id(out)
  validate_opal_fit(out)

  obj$par <- opt$par
  obj$opt <- opt
  .opal_cache_runtime_object(out, obj)
  out
}

#' Validate a portable opal fit
#'
#' Performs structural validation without rebuilding the RTMB objective.
#'
#' @param x Object expected to inherit from `opal_fit`.
#' @return `x`, invisibly, when valid.
#' @export
#'
validate_opal_fit <- function(x) {
  if (!inherits(x, "opal_fit")) {
    stop("`x` must inherit from `opal_fit`.", call. = FALSE)
  }
  if (!identical(as.integer(x$schema_version), .opal_fit_schema_version)) {
    stop(
      "Unsupported `opal_fit` schema version: ", x$schema_version, ".",
      call. = FALSE
    )
  }
  required <- c(
    "schema_version", "model", "data", "parameters", "map", "random",
    "makeadfun_args", "bounds", "control", "optimization", "fit", "mcmc",
    "derived", "provenance", "runtime_id"
  )
  missing <- setdiff(required, names(x))
  if (length(missing)) {
    stop("`opal_fit` is missing: ", paste(missing, collapse = ", "), ".",
         call. = FALSE)
  }
  .opal_validate_named_list(x$data, "data")
  .opal_validate_named_list(x$parameters, "parameters")
  .opal_validate_named_list(x$map, "map")
  .opal_validate_named_list(x$makeadfun_args, "makeadfun_args")
  .opal_validate_named_list(x$derived, "derived")

  if (.opal_contains_nonportable(.opal_fit_runtime_payload(x))) {
    stop("`opal_fit` contains non-portable state.", call. = FALSE)
  }
  if (!is.character(x$random) || any(!nzchar(x$random)) ||
      anyDuplicated(x$random) ||
      length(setdiff(x$random, names(x$parameters)))) {
    stop("Invalid random-effect specification in `opal_fit`.", call. = FALSE)
  }
  if (length(setdiff(names(x$map), names(x$parameters)))) {
    stop("The stored map contains unknown parameters.", call. = FALSE)
  }
  for (name in names(x$map)) {
    if (!is.factor(x$map[[name]]) ||
        length(x$map[[name]]) != length(x$parameters[[name]]) ||
        (!is.null(dim(x$map[[name]])) &&
         !identical(
           as.integer(dim(x$map[[name]])),
           as.integer(dim(x$parameters[[name]]))
         ))) {
      stop("The stored map entry for `", name, "` is invalid.", call. = FALSE)
    }
  }

  opt <- x$fit$opt
  if (!is.list(x$fit) || !is.list(opt) || !is.numeric(opt$par) ||
      is.null(names(opt$par)) || any(!is.finite(opt$par)) ||
      !is.numeric(opt$objective) || length(opt$objective) != 1L ||
      !is.finite(opt$objective) ||
      !is.numeric(opt$convergence) || length(opt$convergence) != 1L ||
      !is.finite(opt$convergence)) {
    stop("Invalid optimizer result in `opal_fit`.", call. = FALSE)
  }
  .opal_validate_named_list(x$fit$diagnostics, "fit diagnostics")

  active_names <- .opal_expand_parameter_names(names(opt$par))
  if (!is.null(x$bounds)) {
    if (!is.list(x$bounds) ||
        !identical(names(x$bounds), c("lower", "upper")) ||
        length(x$bounds$lower) != length(opt$par) ||
        length(x$bounds$upper) != length(opt$par) ||
        !identical(names(x$bounds$lower), active_names) ||
        !identical(names(x$bounds$upper), active_names) ||
        anyNA(x$bounds$lower) || anyNA(x$bounds$upper) ||
        any(x$bounds$lower > x$bounds$upper) ||
        any(opt$par < x$bounds$lower - 1e-8 |
            opt$par > x$bounds$upper + 1e-8)) {
      stop("Invalid bounds in `opal_fit`.", call. = FALSE)
    }
  }

  if (!is.null(x$mcmc)) {
    if (!inherits(x$mcmc, "opal_mcmc") ||
        !is.array(x$mcmc$samples) ||
        length(dim(x$mcmc$samples)) != 3L ||
        !identical(
          .opal_mcmc_sample_names(x$mcmc$samples),
          x$mcmc$sample_names
        ) ||
        !identical(
          x$mcmc$par_names,
          setdiff(x$mcmc$sample_names, "lp__")
        ) ||
        !x$mcmc$parameter_scope %in% c("active", "complete") ||
        x$mcmc$warmup + x$mcmc$iter != dim(x$mcmc$samples)[1L] ||
        x$mcmc$chains != dim(x$mcmc$samples)[2L]) {
      stop("Invalid MCMC payload in `opal_fit`.", call. = FALSE)
    }
    if (identical(x$mcmc$parameter_scope, "active") &&
        !identical(x$mcmc$par_names, active_names)) {
      stop("Stored MCMC parameters do not match the saved active layout.",
           call. = FALSE)
    }
  }

  if (!is.list(x$model) ||
      !identical(x$model$name, "opal_model") ||
      length(x$model$schema_version) != 1L ||
      length(x$model$scientific_version) != 1L ||
      length(x$model$signature) != 1L) {
    stop("Invalid model metadata in `opal_fit`.", call. = FALSE)
  }
  if (!is.list(x$optimization) ||
      !is.character(x$optimization$method) ||
      length(x$optimization$method) != 1L ||
      !nzchar(x$optimization$method) ||
      (!is.null(x$control) && !is.list(x$control)) ||
      !is.list(x$provenance) ||
      !is.list(x$provenance$metadata) ||
      !is.character(x$runtime_id) ||
      length(x$runtime_id) != 1L ||
      !grepl("^opal-fit-[0-9a-f]{32}$", x$runtime_id)) {
    stop("Invalid optimization or provenance metadata in `opal_fit`.",
         call. = FALSE)
  }
  invisible(x)
}

#' Check compatibility of a saved opal fit
#'
#' @param x An `opal_fit` object.
#' @return A list containing `compatible`, `errors`, and `warnings`.
#' @export
#'
opal_fit_compatibility <- function(x) {
  validate_opal_fit(x)
  current <- .opal_model_metadata()
  errors <- character()
  warnings <- character()

  if (!identical(x$model$name, current$name) ||
      !identical(
        as.integer(x$model$schema_version),
        as.integer(current$schema_version)
      ) ||
      !identical(x$model$scientific_version, current$scientific_version)) {
    errors <- c(
      errors,
      "The saved fit uses a different opal scientific model contract."
    )
  }
  if (!identical(x$model$signature, current$signature)) {
    warnings <- c(
      warnings,
      "The opal_model implementation checksum differs from the saved fit."
    )
  }
  if (!identical(x$model$package_version, current$package_version)) {
    warnings <- c(
      warnings,
      paste0(
        "The saved fit used opal ", x$model$package_version,
        "; the current version is ", current$package_version, "."
      )
    )
  }
  list(
    compatible = !length(errors),
    errors = errors,
    warnings = warnings,
    saved = x$model,
    current = current
  )
}

#' Rebuild the RTMB objective for a portable opal fit
#'
#' Recreates the objective from stored data, fitted parameters, map, and random
#' effects, then verifies the active parameter layout and saved objective.
#'
#' @param x An `opal_fit` object.
#' @param strict Stop on compatibility or parameter-order differences.
#' @param check_objective Compare the rebuilt and saved objective values.
#' @param tolerance Relative tolerance passed to [all.equal()].
#' @param silent Passed to `RTMB::MakeADFun()`.
#' @param cache Cache the rebuilt objective for this R session.
#' @return A newly constructed RTMB objective.
#' @export
#'
rebuild_opal_object <- function(x, strict = TRUE, check_objective = TRUE,
                                tolerance = 1e-6, silent = FALSE,
                                cache = TRUE) {
  validate_opal_fit(x)
  compatibility <- opal_fit_compatibility(x)
  if (length(compatibility$errors)) {
    message <- paste(compatibility$errors, collapse = " ")
    if (isTRUE(strict)) stop(message, call. = FALSE)
    warning(message, call. = FALSE)
  }
  if (length(compatibility$warnings)) {
    warning(paste(compatibility$warnings, collapse = " "), call. = FALSE)
  }

  makeadfun_args <- c(
    list(
      func = cmb(opal_model, x$data),
      parameters = x$parameters,
      map = x$map,
      random = x$random
    ),
    x$makeadfun_args
  )
  makeadfun_args$silent <- silent
  obj <- do.call(RTMB::MakeADFun, makeadfun_args)

  current_names <- .opal_expand_parameter_names(names(obj$par))
  saved_names <- .opal_expand_parameter_names(names(x$fit$opt$par))
  par <- x$fit$opt$par
  if (!identical(current_names, saved_names)) {
    if (isTRUE(strict)) {
      stop("Active parameter names or ordering differ from the saved fit.",
           call. = FALSE)
    }
    if (length(current_names) == length(saved_names) &&
        setequal(current_names, saved_names)) {
      warning("Active parameter order differs; aligning by canonical names.",
              call. = FALSE)
      par <- par[match(current_names, saved_names)]
    } else {
      stop("The saved active parameter layout cannot be aligned.",
           call. = FALSE)
    }
  }
  names(par) <- names(obj$par)

  rebuilt_objective <- as.numeric(obj$fn(par))
  if (isTRUE(check_objective) &&
      !isTRUE(all.equal(
        rebuilt_objective,
        as.numeric(x$fit$opt$objective),
        tolerance = tolerance
      ))) {
    message <- paste0(
      "Rebuilt objective (", signif(rebuilt_objective, 10L),
      ") differs from saved objective (",
      signif(x$fit$opt$objective, 10L), ")."
    )
    if (isTRUE(strict)) stop(message, call. = FALSE)
    warning(message, call. = FALSE)
  }

  obj$par <- par
  obj$opt <- x$fit$opt
  obj$opt$par <- par
  if (!is.null(x$mcmc) &&
      identical(x$mcmc$parameter_scope, "complete") &&
      !identical(
        x$mcmc$par_names,
        .opal_expand_parameter_names(names(obj$env$last.par.best))
      )) {
    stop("Stored MCMC parameters do not match the rebuilt complete layout.",
         call. = FALSE)
  }
  attr(obj, "opal_fit_compatibility") <- compatibility
  if (isTRUE(cache)) .opal_cache_runtime_object(x, obj)
  obj
}

#' Access the runtime objective for an opal fit
#'
#' @param x An `opal_fit` object.
#' @param fresh Construct an isolated objective instead of using the session
#'   cache.
#' @return A fitted RTMB objective.
#' @export
#'
opal_fit_object <- function(x, fresh = FALSE) {
  validate_opal_fit(x)
  if (!is.logical(fresh) || length(fresh) != 1L || is.na(fresh)) {
    stop("`fresh` must be `TRUE` or `FALSE`.", call. = FALSE)
  }
  if (!identical(x$runtime_id, .opal_fit_runtime_id(x))) {
    stop("The fit payload was modified without refreshing its runtime identity.",
         call. = FALSE)
  }
  if (isTRUE(fresh)) {
    return(rebuild_opal_object(x, strict = TRUE, silent = TRUE, cache = FALSE))
  }
  obj <- .opal_cached_runtime_object(x)
  if (is.null(obj)) {
    obj <- rebuild_opal_object(x, strict = TRUE, silent = TRUE, cache = TRUE)
  }
  obj
}

#' Recreate the fitted model report
#'
#' @param x An `opal_fit` object.
#' @return The named report list returned by `opal_model`.
#' @export
#'
opal_fit_report <- function(x) {
  obj <- opal_fit_object(x)
  report <- obj$report(obj$env$last.par.best)
  if (!is.list(report) || is.null(names(report))) {
    stop("The rebuilt model did not return a named report list.",
         call. = FALSE)
  }
  report
}

#' Update portable results attached to an opal fit
#'
#' Adds or replaces normalized MCMC output and merges portable diagnostics,
#' derived results, or metadata. The fitted model state is unchanged.
#'
#' @param x An `opal_fit` object.
#' @param mcmc Optional replacement MCMC output. If omitted, existing output is
#'   retained; explicitly supply `NULL` to remove it.
#' @param mcmc_settings Optional settings merged into a replacement MCMC
#'   payload.
#' @param diagnostics Optional named list merged into fit diagnostics.
#' @param derived Optional named list merged into derived results.
#' @param metadata Optional named list merged into user metadata.
#' @return An updated `opal_fit`.
#' @export
#'
update_opal_fit <- function(x, mcmc, mcmc_settings = list(),
                            diagnostics = list(), derived = list(),
                            metadata = list()) {
  validate_opal_fit(x)
  .opal_validate_named_list(mcmc_settings, "mcmc_settings")
  .opal_validate_named_list(diagnostics, "diagnostics")
  .opal_validate_named_list(derived, "derived")
  .opal_validate_named_list(metadata, "metadata")

  obj <- opal_fit_object(x)
  if (!missing(mcmc)) {
    x["mcmc"] <- list(
      .opal_attach_mcmc(
        mcmc,
        obj = obj,
        settings = mcmc_settings
      )
    )
  }
  x$fit$diagnostics <- utils::modifyList(x$fit$diagnostics, diagnostics)
  x$derived <- utils::modifyList(x$derived, derived)
  x$provenance$metadata <- utils::modifyList(
    x$provenance$metadata,
    metadata
  )
  x$provenance$updated_at <- format(Sys.time(), tz = "UTC", usetz = TRUE)

  if (.opal_contains_nonportable(.opal_fit_runtime_payload(x))) {
    stop("Updated fit state contains non-portable values.", call. = FALSE)
  }
  x$runtime_id <- .opal_fit_runtime_id(x)
  validate_opal_fit(x)
  .opal_cache_runtime_object(x, obj)
  x
}

#' Convert normalized opal posterior draws to a tmbfit
#'
#' @param x An `opal_fit` with MCMC output or an `opal_mcmc` object.
#' @return A list inheriting from `tmbfit`.
#' @export
#'
opal_as_tmbfit <- function(x) {
  if (inherits(x, "opal_fit")) {
    validate_opal_fit(x)
    if (is.null(x$mcmc)) {
      stop("This `opal_fit` has no stored MCMC output.", call. = FALSE)
    }
    x <- x$mcmc
  }
  if (!inherits(x, "opal_mcmc")) {
    stop("`x` must be an `opal_fit` or `opal_mcmc`.", call. = FALSE)
  }

  samples <- x$samples
  sampler_params <- x$sampler_params
  samples_unbounded <- x$samples_unbounded
  warmup <- x$warmup
  if (warmup == 0L) {
    synthetic_index <- c(1L, seq_len(dim(samples)[1L]))
    samples <- samples[synthetic_index, , , drop = FALSE]
    if (!is.null(sampler_params)) {
      sampler_params <- lapply(
        sampler_params,
        function(z) z[synthetic_index, , drop = FALSE]
      )
    }
    if (!is.null(samples_unbounded) &&
        is.array(samples_unbounded) &&
        length(dim(samples_unbounded)) == 3L) {
      samples_unbounded <- samples_unbounded[
        synthetic_index,
        ,
        ,
        drop = FALSE
      ]
    }
    warmup <- 1L
  }

  out <- list(
    samples = samples,
    sampler_params = sampler_params,
    mle = x$mle,
    monitor = x$monitor,
    model = x$model,
    metric = x$metric,
    par_names = x$par_names,
    max_treedepth = x$max_treedepth,
    warmup = warmup,
    iter = x$iter,
    thin = x$thin,
    timing = x$timing,
    algorithm = x$algorithm,
    samples_unbounded = samples_unbounded,
    inits = x$inits
  )
  class(out) <- c("tmbfit", "list")
  out
}

#' Save and read portable opal fits
#'
#' `save_opal_fit()` writes atomically and refuses to replace an existing file
#' unless requested. `read_opal_fit()` verifies the portable payload and can
#' rebuild its transient RTMB objective.
#'
#' @param x An `opal_fit` object.
#' @param file Path to an RDS file.
#' @param compress Compression passed to [saveRDS()].
#' @param overwrite Replace an existing file.
#' @param strict Treat scientific-contract incompatibility as an error.
#' @param rebuild Rebuild and cache the RTMB objective after reading. Defaults
#'   to `strict`.
#' @return `save_opal_fit()` invisibly returns the normalized path;
#'   `read_opal_fit()` returns an `opal_fit`.
#' @name opal_fit_io
NULL

#' @rdname opal_fit_io
#' @export
#'
save_opal_fit <- function(x, file, compress = "gzip", overwrite = FALSE) {
  validate_opal_fit(x)
  x$runtime_id <- .opal_fit_runtime_id(x)
  file <- path.expand(file)
  directory <- dirname(file)
  if (!dir.exists(directory)) {
    stop("Output directory does not exist: ", directory, call. = FALSE)
  }
  if (file.exists(file) && !isTRUE(overwrite)) {
    stop(
      "File already exists; set `overwrite = TRUE` to replace it: ",
      file,
      call. = FALSE
    )
  }

  temporary_file <- tempfile(
    pattern = ".opal-fit-",
    tmpdir = directory,
    fileext = ".rds"
  )
  on.exit(unlink(temporary_file), add = TRUE)
  saveRDS(x, file = temporary_file, compress = compress, version = 3L)

  backup_file <- NULL
  if (file.exists(file)) {
    backup_file <- tempfile(
      pattern = ".opal-fit-backup-",
      tmpdir = directory,
      fileext = ".rds"
    )
    if (!file.rename(file, backup_file)) {
      stop("Could not stage the existing fit for replacement.", call. = FALSE)
    }
  }
  if (!file.rename(temporary_file, file)) {
    restored <- is.null(backup_file) || file.rename(backup_file, file)
    if (restored) {
      stop("Could not move the completed fit into place.", call. = FALSE)
    }
    stop(
      "Could not install the new fit or restore the previous file; backup: ",
      backup_file,
      call. = FALSE
    )
  }
  if (!is.null(backup_file) && unlink(backup_file) != 0L) {
    warning("The previous fit backup could not be removed: ", backup_file,
            call. = FALSE)
  }
  invisible(normalizePath(file, mustWork = TRUE))
}

#' @rdname opal_fit_io
#' @export
#'
read_opal_fit <- function(file, strict = FALSE, rebuild = strict) {
  if (!is.logical(strict) || length(strict) != 1L || is.na(strict) ||
      !is.logical(rebuild) || length(rebuild) != 1L || is.na(rebuild)) {
    stop("`strict` and `rebuild` must each be `TRUE` or `FALSE`.",
         call. = FALSE)
  }
  x <- readRDS(file)
  validate_opal_fit(x)
  if (!identical(x$runtime_id, .opal_fit_runtime_id(x))) {
    stop("The saved fit runtime identity does not match its payload.",
         call. = FALSE)
  }

  compatibility <- opal_fit_compatibility(x)
  messages <- c(compatibility$errors, compatibility$warnings)
  if (length(compatibility$errors) && isTRUE(strict)) {
    stop(paste(messages, collapse = " "), call. = FALSE)
  }
  if (length(messages)) warning(paste(messages, collapse = " "), call. = FALSE)
  if (isTRUE(rebuild)) {
    if (length(compatibility$errors)) {
      stop("Cannot rebuild an incompatible opal fit.", call. = FALSE)
    }
    rebuild_opal_object(x, strict = TRUE, silent = TRUE, cache = TRUE)
  }
  x
}

.opal_estimability_status <- function(x) {
  if (is.null(x)) return("not run")
  .opal_fit_or(x$message, x$status)
}

#' @export
print.opal_fit <- function(x, ...) {
  validate_opal_fit(x)
  opt <- x$fit$opt
  cat("<opal_fit>\n")
  cat("  Model:        ", x$model$name, " (schema ",
      x$model$schema_version, ")\n", sep = "")
  cat("  opal version: ", x$model$package_version, "\n", sep = "")
  cat("  Parameters:   ", length(opt$par), " active\n", sep = "")
  cat("  Objective:    ", format(signif(opt$objective, 10L)), "\n", sep = "")
  cat("  Convergence:  ", opt$convergence, "\n", sep = "")
  cat("  Estimability: ",
      .opal_estimability_status(x$fit$estimability), "\n", sep = "")
  if (is.null(x$mcmc)) {
    cat("  MCMC:         not stored\n")
  } else {
    cat("  MCMC:         stored [",
        paste(dim(x$mcmc$samples), collapse = " x "), "]\n", sep = "")
  }
  cat("  Derived sets: ", length(x$derived), "\n", sep = "")
  cat("  Created:      ", x$provenance$created_at, "\n", sep = "")
  invisible(x)
}

#' @export
summary.opal_fit <- function(object, ...) {
  validate_opal_fit(object)
  opt <- object$fit$opt
  out <- list(
    model = data.frame(
      model = object$model$name,
      model_schema = object$model$schema_version,
      scientific_version = object$model$scientific_version,
      opal_version = object$model$package_version,
      stringsAsFactors = FALSE
    ),
    optimization = data.frame(
      method = object$optimization$method,
      n_parameters = length(opt$par),
      objective = opt$objective,
      convergence = opt$convergence,
      message = .opal_fit_or(opt$message, NA_character_),
      stringsAsFactors = FALSE
    ),
    estimability = .opal_estimability_status(object$fit$estimability),
    mcmc_dimensions = if (is.null(object$mcmc)) {
      NULL
    } else {
      dim(object$mcmc$samples)
    },
    derived = names(object$derived),
    metadata = object$provenance$metadata
  )
  class(out) <- c("summary.opal_fit", "list")
  out
}

#' @export
print.summary.opal_fit <- function(x, ...) {
  cat("opal fitted-model summary\n\n")
  print(x$model, row.names = FALSE)
  cat("\nOptimization\n")
  print(x$optimization, row.names = FALSE)
  cat("\nEstimability: ", x$estimability, "\n", sep = "")
  if (!is.null(x$mcmc_dimensions)) {
    cat("MCMC dimensions: ",
        paste(x$mcmc_dimensions, collapse = " x "), "\n", sep = "")
  }
  if (length(x$derived)) {
    cat("Derived results: ", paste(x$derived, collapse = ", "), "\n", sep = "")
  }
  invisible(x)
}

#' @export
print.opal_mcmc <- function(x, ...) {
  cat("<opal_mcmc>\n")
  cat("  Samples:   ", paste(dim(x$samples), collapse = " x "), "\n",
      sep = "")
  cat("  Algorithm: ", x$algorithm, "\n", sep = "")
  cat("  Scope:     ", .opal_fit_or(x$parameter_scope, "unattached"), "\n",
      sep = "")
  invisible(x)
}
