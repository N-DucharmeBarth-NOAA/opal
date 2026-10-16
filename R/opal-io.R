#' Convert a legacy Opal fit to the staged object workflow
#'
#' Verifies the saved fitted objective without optimising or sampling. Existing
#' diagnostics, posterior draws, derived outputs, and provenance are preserved.
#' @param x A legacy `opal_fit` or an `opal_obj`.
#' @param integrity Legacy runtime verification mode; portable permits verified cross-R reads.
#' @return An `opal_obj`.
#' @family assessment workflow
#' @examples
#' legacy <- read_opal_fit(system.file("extdata", "opaka_quickstart_fit.rds",
#'                                    package = "opal"), integrity = "portable")
#' assessment <- opal_from_fit(legacy)
#' summary(assessment)
#' @export
opal_from_fit <- function(x, integrity = c("portable", "exact")) {
  if (inherits(x, "opal_obj")) {
    validate_opal_obj(x, results = TRUE)
    return(x)
  }
  integrity <- match.arg(integrity)
  validate_opal_fit(x)
  runtime <- opal_fit_object(x, fresh = TRUE, integrity = integrity)
  out <- opal_obj(x$data, x$parameters, x$map, x$random, x$bounds,
                   x$control, makeadfun_args = x$makeadfun_args,
                   metadata = x$provenance$metadata)
  out$model <- x$model
  out$build <- list(active_names = .opal_expand_parameter_names(names(runtime$par)),
                    objective = as.numeric(runtime$fn(runtime$par)))
  out$bounds <- .opal_normalize_bounds(
    .opal_fit_or(x$bounds, get_bounds(runtime, x$parameters)), runtime$par)
  out$fit <- x$fit
  out$fit$parameters <- x$parameters
  out$fit$last_par_best <- runtime$env$last.par.best
  out$provenance <- x$provenance
  out$provenance$migrated_from <- list(class = "opal_fit", schema_version = x$schema_version)
  out <- .opal_obj_seal(out)
  if (!is.null(x$mcmc)) out <- opal_attach_mcmc(out, x$mcmc, x$mcmc$settings, check = FALSE)
  out$derived <- x$derived
  out <- .opal_obj_seal(out)
  assign(.opal_obj_runtime_id(out), runtime, .opal_obj_cache)
  out
}

#' Save and read staged Opal objects
#'
#' Every stage is portable. Saved objects contain a checksum of the complete
#' payload; their transient RTMB objective is never serialised. Reads verify
#' the checksum and scientific contract. Legacy fits are converted with
#' objective verification. No read operation optimises or samples.
#' @param x An `opal_obj` or legacy `opal_fit` to save.
#' @param file RDS path.
#' @param compress RDS compression.
#' @param overwrite Allow replacement of an existing file.
#' @param rebuild Rebuild a configured object's runtime after reading.
#' @param strict Reject incompatible scientific contracts. FALSE permits
#'   inspection only; runtime access still rejects an incompatible model.
#' @param integrity Verification mode for legacy files: `"portable"` checks
#'   the rebuilt objective across R versions; `"exact"` also requires the
#'   original runtime checksum. New-format files always verify their payload
#'   checksum, independently of this legacy option.
#' @return Save returns the path invisibly; read returns an `opal_obj`.
#' @name opal_io
#' @details
#' Assign the return value of `opal_read()` to resume work in another session.
#' Use `rebuild = TRUE` to verify runtime reconstruction immediately. Existing
#' files are protected by default; set `overwrite = TRUE` deliberately when
#' saving an updated assessment to the same path.
#' @family assessment workflow
#' @examples
#' inputs <- opaka_quickstart_inputs()
#' assessment <- opal_obj(inputs$data, inputs$parameters, inputs$map)
#' path <- tempfile(fileext = ".rds")
#' opal_save(assessment, path)
#' restored <- opal_read(path, rebuild = TRUE)
#' summary(restored)
#' unlink(path)
#' @export
opal_save <- function(x, file, compress = "gzip", overwrite = FALSE) {
  .opal_obj_flag(overwrite, "overwrite")
  x <- opal_from_fit(x)
  validate_opal_obj(x, results = TRUE)
  envelope <- list(format = "opal_obj", version = 1L, object = x,
                    checksum = .opal_obj_hash(x))
  if (!is.character(file) || length(file) != 1L || is.na(file) || !nzchar(file)) stop("Supply one file path.")
  file <- path.expand(file)
  if (dir.exists(file)) stop("The output path is a directory.", call. = FALSE)
  if (!dir.exists(dirname(file))) stop("Output directory does not exist.", call. = FALSE)
  if (file.exists(file) && !overwrite) stop("File already exists; use overwrite = TRUE.", call. = FALSE)
  temporary <- tempfile(".opal-", tmpdir = dirname(file))
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(envelope, temporary, compress = compress, version = 3L)
  backup <- NULL
  if (file.exists(file)) {
    backup <- tempfile(".opal-backup-", tmpdir = dirname(file))
    if (!file.rename(file, backup)) stop("Could not back up the existing file.", call. = FALSE)
  }
  if (!file.rename(temporary, file)) {
    restored <- is.null(backup) || file.rename(backup, file)
    stop(if (restored) "Could not install the saved object; previous file restored." else
      paste("Could not restore the previous file; backup:", backup), call. = FALSE)
  }
  if (!is.null(backup)) unlink(backup)
  invisible(normalizePath(file))
}

#' @rdname opal_io
#' @export
opal_read <- function(file, rebuild = FALSE, strict = TRUE,
                      integrity = c("portable", "exact")) {
  .opal_obj_flag(rebuild, "rebuild")
  .opal_obj_flag(strict, "strict")
  integrity <- match.arg(integrity)
  saved <- readRDS(file)
  if (inherits(saved, "opal_fit")) return(opal_from_fit(saved, integrity))
  if (!is.list(saved) || !identical(saved$format, "opal_obj") ||
      !identical(saved$version, 1L) ||
      !identical(saved$checksum, .opal_obj_hash(saved$object))) {
    stop("Invalid or modified Opal file payload.", call. = FALSE)
  }
  x <- saved$object
  validate_opal_obj(x, results = TRUE)
  compatible <- tryCatch({ .opal_obj_compatible(x); TRUE }, error = identity)
  if (inherits(compatible, "error")) {
    if (strict || rebuild) stop(conditionMessage(compatible), call. = FALSE)
    warning(conditionMessage(compatible), call. = FALSE)
  }
  if (rebuild) x <- opal_build(x)
  x
}
