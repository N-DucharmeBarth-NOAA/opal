# Diagnostic contract is independent of the numerical model contract.
.opal_validation_version <- "opal_validation_v2"

.opal_validation_current <- function(x, scope) {
  record <- x$validation[[scope]]
  !is.null(record) && identical(record$version, .opal_validation_version) &&
    identical(record$identity, .opal_check_identity(x, scope)) &&
    identical(record$payload_id, .opal_obj_hash(record[setdiff(names(record), "payload_id")]))
}

#' Diagnose biological feasibility of an Opal model
#'
#' Checks the initial-equilibrium and harvest penalties, population states,
#' biological inputs, catch reconstruction, and derived depletion. A numerical
#' optimum alone does not imply that these checks pass.
#' @param x A configured or fitted `opal_obj`.
#' @param penalty_tolerance Maximum allowed positive continuation penalty.
#' @param catch_tolerance Maximum catch error divided by `pmax(1, observed)`.
#' @return A list with `passes`, named logical `checks`, and numerical `metrics`.
#' @details Zero abundance is allowed, but unfished spawning output must be
#'   positive. Catch is conditioned on, so its difference is a reconstruction
#'   check, not an observation residual. Checks use the reported, draw-specific
#'   biology. They do not assess identifiability or scientific model adequacy.
#' @family assessment diagnostics
#' @export
opal_diagnose <- function(x, penalty_tolerance = 1e-10, catch_tolerance = 1e-6) {
  object <- opal_rtmb(x, fresh = TRUE)
  p <- object$env$last.par.best
  .opal_diagnose_report(x$data, object$env$parList(par = p), object$report(p),
                        penalty_tolerance, catch_tolerance)
}

.opal_diagnose_report <- function(data, parameters, report,
                                  penalty_tolerance = 1e-10, catch_tolerance = 1e-6) {
  tolerance <- c(penalty_tolerance, catch_tolerance)
  if (length(tolerance) != 2L || any(!is.finite(tolerance)) || any(tolerance < 0)) {
    stop("Biological tolerances must be non-negative finite scalars.", call. = FALSE)
  }
  valid <- function(z, positive = FALSE) is.numeric(z) && length(z) > 0L &&
    all(is.finite(z)) && all(if (positive) z > 0 else z >= 0)
  checks <- c(
    initial_equilibrium = valid(report$lp_init_penalty) && report$lp_init_penalty <= penalty_tolerance,
    harvest_penalty = valid(report$lp_penalty) && report$lp_penalty <= penalty_tolerance,
    population = valid(report$number_ysa) && valid(report$number0_ysa),
    mortality = valid(report$M_a),
    maturity = valid(report$maturity_a) && all(report$maturity_a <= 1),
    spawning_potential = valid(report$spawning_potential_a) && any(report$spawning_potential_a > 0),
    weight = valid(report$weight_fya_mod),
    selectivity = valid(report$sel_fya) && all(report$sel_fya <= 1 + 1e-10),
    steepness = valid(exp(parameters$log_h)) && all(exp(parameters$log_h) > 0.2 & exp(parameters$log_h) <= 1),
    recruitment = valid(report$R0, TRUE) && valid(report$sigma_r, TRUE),
    spawning_output = valid(report$B0, TRUE) && valid(report$spawning_biomass_y) &&
      valid(report$spawning_biomass0_y, TRUE),
    harvest = valid(report$hrate_ysa) && all(report$hrate_ysa <= 1),
    depletion = valid(report$static_depletion_y) && valid(report$dynamic_depletion_y))
  catch_error <- Inf
  if (identical(dim(report$catch_pred_ysf), dim(data$catch_obs_ysf)) &&
      valid(report$catch_pred_ysf) && valid(data$catch_obs_ysf)) {
    catch_error <- max(abs(report$catch_pred_ysf - data$catch_obs_ysf) /
                         pmax(1, data$catch_obs_ysf))
  }
  checks <- c(checks, catch_reconstruction = is.finite(catch_error) && catch_error <= catch_tolerance)
  list(passes = all(checks), checks = checks, metrics = list(
    initial_penalty = report$lp_init_penalty, total_penalty = report$lp_penalty,
    max_harvest = if (valid(report$hrate_ysa)) max(report$hrate_ysa) else NA_real_,
    max_relative_catch_error = catch_error))
}

# All posterior consumers share one draw order and require sampled latent states.
.opal_posterior_context <- function(x, draws = NULL) {
  validate_opal_obj(x, results = TRUE)
  if (is.null(x$mcmc)) stop("No posterior draws are stored.", call. = FALSE)
  m <- x$mcmc
  if (length(x$random) && m$parameter_scope != "complete") {
    stop("Biological posterior analysis requires complete joint draws, including random effects.", call. = FALSE)
  }
  keep <- seq.int(m$warmup + 1L, dim(m$samples)[1L])
  samples <- m$samples[keep, , m$par_names, drop = FALSE]
  values <- matrix(samples, ncol = length(m$par_names), dimnames = list(NULL, m$par_names))
  ids <- expand.grid(iteration = keep, chain = seq_len(m$chains))
  ids$draw <- seq_len(nrow(ids))
  if (is.null(draws)) draws <- ids$draw
  if (!is.numeric(draws) || !length(draws) || any(!is.finite(draws)) ||
      any(draws != floor(draws)) || any(draws < 1 | draws > nrow(ids)) || anyDuplicated(draws)) {
    stop("`draws` must select distinct retained draw indices.", call. = FALSE)
  }
  object <- opal_rtmb(x, fresh = TRUE)
  expected <- .opal_expand_parameter_names(names(object$env$last.par.best))
  if (!identical(expected, colnames(values))) stop("Posterior parameter layout is incompatible.", call. = FALSE)
  list(object = object, values = values[draws, , drop = FALSE], ids = ids[draws, , drop = FALSE])
}

.opal_check_draws <- function(x, penalty_tolerance, catch_tolerance) {
  context <- .opal_posterior_context(x)
  fixed <- match(.opal_expand_parameter_names(names(x$bounds$lower)), colnames(context$values))
  failures <- vector("list", nrow(context$values))
  for (i in seq_len(nrow(context$values))) {
    p <- context$values[i, ]
    reasons <- character()
    if (any(p[fixed] < x$bounds$lower | p[fixed] > x$bounds$upper)) reasons <- "bounds"
    biological <- tryCatch(.opal_diagnose_report(x$data, context$object$env$parList(par = p),
      context$object$report(p), penalty_tolerance, catch_tolerance), error = identity)
    if (inherits(biological, "error")) {
      reasons <- c(reasons, paste0("report: ", conditionMessage(biological)))
    } else reasons <- c(reasons, names(biological$checks)[!biological$checks])
    if (length(reasons)) failures[[i]] <- cbind(context$ids[i, ], reason = paste(reasons, collapse = "; "))
  }
  failed <- do.call(rbind, failures)
  if (is.null(failed)) failed <- cbind(context$ids[FALSE, ], reason = character())
  list(passes = nrow(failed) == 0L, checked = nrow(context$values), failures = failed,
       scope = "all retained joint draws")
}
