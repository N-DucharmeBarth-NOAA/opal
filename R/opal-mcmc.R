#' Attach posterior draws to an Opal object
#' @param x A configured or fitted `opal_obj`.
#' @param mcmc SparseNUTS output, sample matrix, or iteration-chain-variable array.
#' @param settings Named portable sampler settings.
#' @param check Run posterior diagnostics.
#' @param check_args Arguments passed to [opal_check()].
#' @return An updated `opal_obj`. Previous selected draws are retained in history.
#' @export
opal_attach_mcmc <- function(x, mcmc, settings = list(), check = TRUE,
                             check_args = list()) {
  validate_opal_obj(x, results = TRUE)
  x <- opal_build(x)
  .opal_obj_flag(check, "check")
  .opal_validate_named_list(settings, "settings")
  .opal_obj_portable(settings, "settings")
  candidate <- .opal_attach_mcmc(mcmc, opal_rtmb(x, fresh = TRUE), settings)
  if (is.null(candidate)) stop("Supply posterior draws.", call. = FALSE)
  candidate$target_id <- x$identity$target
  candidate$fit_id <- x$identity$fit
  candidate$payload_id <- .opal_obj_hash(unclass(candidate))
  old <- x
  x$mcmc <- candidate
  x$validation$mcmc <- NULL
  x$derived <- list()
  if (check) x <- .opal_run_check(x, "mcmc", check_args)
  # Never replace a checked, passing posterior with a failed or unchecked run.
  old_passed <- !is.null(old$mcmc) && isTRUE(old$validation$mcmc$passes) &&
    identical(old$validation$mcmc$identity, .opal_check_identity(old, "mcmc"))
  rejected <- old_passed && !isTRUE(x$validation$mcmc$passes)
  archived <- if (rejected) x else old
  selected <- if (rejected) old else x
  if (!is.null(archived$mcmc)) {
    selected$mcmc_history <- c(old$mcmc_history, list(list(
      status = if (rejected) "not_selected" else "previous",
      mcmc = archived$mcmc, validation = archived$validation$mcmc)))
  }
  .opal_obj_seal(selected)
}

#' Run MCMC for an Opal object
#'
#' Sampling uses a fresh RTMB objective. It never changes the stored optimum.
#' `init = "auto"` uses the fitted point when available and the configured
#' active parameters otherwise. Failed sampler attempts are recorded and
#' warned about, preserving any existing posterior. The default sampler skips
#' internal optimisation; an MLE is optional. The default stan metric supports
#' bounds. Random effects are sampled jointly by default, over their full
#' domain, while fixed effects retain the object's bounds. Change bounds
#' with opal_update(), not sampler arguments. Use laplace = TRUE and metric = 'unit'
#' for marginal sampling.
#' @param x A configured or fitted `opal_obj`.
#' @param sampler `"snuts"` or a function accepting an `obj` argument and sampler settings.
#' @param init Initial-value policy or explicit sampler initial values.
#' @param check Run posterior diagnostics.
#' @param check_args Arguments passed to [opal_check()].
#' @param ... Named sampler settings, including seed, chains, cores,
#'   num_samples, and num_warmup. The obj, globals, lower, and upper arguments are reserved.
#' @return An updated `opal_obj`, including portable samples and attempt history.
#' @export
opal_mcmc <- function(x, sampler = "snuts", init = "auto", check = TRUE,
                      check_args = list(), ...) {
  validate_opal_obj(x, results = TRUE)
  x <- opal_build(x)
  .opal_obj_flag(check, "check")
  .opal_validate_named_list(check_args, "check_args")
  supplied <- list(...)
  .opal_validate_named_list(supplied, "sampler arguments")
  if (length(intersect(names(supplied), c("obj", "globals", "lower", "upper")))) {
    stop("Sampler objective, globals, and bounds are managed by opal_mcmc(); change bounds with opal_update().", call. = FALSE)
  }
  fun <- if (is.function(sampler)) sampler else if (identical(sampler, "snuts")) {
    SparseNUTS::sample_snuts
  } else stop("Unknown sampler; supply 'snuts' or a function.", call. = FALSE)
  object <- opal_rtmb(x, fresh = TRUE)
  args <- utils::modifyList(list(
    init = init, metric = "stan", skip_optimization = TRUE, laplace = FALSE,
    num_samples = 500L, num_warmup = 1000L, chains = 4L, cores = 4L), supplied)
  if (length(x$random) && isTRUE(args$laplace) && identical(args$metric, "stan")) {
    stop("SparseNUTS's stan metric requires joint sampling; use metric = 'unit' for Laplace sampling.", call. = FALSE)
  }
  bounds <- x$bounds
  if (length(x$random) && !isTRUE(args$laplace)) {
    # Sampling the joint density requires a full parameter vector and bounds.
    parameters <- object$env$parList(par = object$env$last.par.best)
    object <- do.call(RTMB::MakeADFun, c(list(func = cmb(opal_model, x$data),
      parameters = parameters, map = x$map, random = NULL, silent = TRUE), x$makeadfun_args))
    bounds <- list(lower = stats::setNames(rep(-Inf, length(object$par)), names(object$par)),
                   upper = stats::setNames(rep(Inf, length(object$par)), names(object$par)))
    fixed <- match(.opal_expand_parameter_names(names(x$bounds$lower)),
                   .opal_expand_parameter_names(names(object$par)))
    bounds$lower[fixed] <- x$bounds$lower
    bounds$upper[fixed] <- x$bounds$upper
  }
  if (any(is.finite(c(bounds$lower, bounds$upper)))) {
    args$lower <- unname(bounds$lower)
    args$upper <- unname(bounds$upper)
  }
  if (identical(init, "auto")) args$init <- "last.par.best"
  # MakeADFun has already evaluated the initial state, so this also works
  # without an MLE. Explicit numeric starts are replicated per chain.
  if (is.numeric(args$init)) args$init <- rep(list(args$init), args$chains)
  if (is.null(args$seed)) args$seed <- sample.int(.Machine$integer.max, 1L)
  settings <- args
  settings$sampler <- if (is.function(sampler)) {
    list(type = "function", signature = .opal_text_checksum(.opal_function_text(sampler)))
  } else list(type = "package", name = "SparseNUTS::sample_snuts")
  .opal_obj_portable(settings, "sampler settings")
  result <- tryCatch(do.call(fun, c(list(obj = object, globals = opal_globals()), args)),
                     error = identity)
  if (inherits(result, "error")) {
    x$mcmc_history <- c(x$mcmc_history, list(list(status = "error",
      message = conditionMessage(result), settings = settings, mcmc = NULL,
      target_id = x$identity$target, fit_id = x$identity$fit)))
    warning("MCMC failed; previous results retained: ", conditionMessage(result), call. = FALSE)
    return(.opal_obj_seal(x))
  }
  opal_attach_mcmc(x, result, settings, check, check_args)
}
