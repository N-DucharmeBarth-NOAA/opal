#' Calculate deterministic equilibrium reference points
#'
#' Calculates MSY under a fixed fleet allocation, selectivity, and biology.
#' This is an equilibrium calculation, independent of the projection engine.
#' @param x A fitted or sampled `opal_obj`.
#' @param fleet_weights Non-negative fleet allocation weights, normalised to
#'   sum to one. Required because the fleet mix is a scientific choice.
#' @param year Model-year index supplying selectivity and weight at age.
#' @param uncertainty Use the fitted point or complete joint MCMC draws.
#' @param draws Optional retained draw indices when `uncertainty = "mcmc"`.
#' @param u_max Maximum scalar seasonal harvest intensity, in `(0, 1]`.
#' @param grid_size Number of harvest intensities used to bracket maxima.
#' @param name Name of the stored result.
#' @return An updated object. [opal_derived()] returns per-draw MSY, spawning
#'   output at MSY, recruitment at MSY, seasonal intensity at MSY, depletion
#'   at MSY, terminal spawning output relative to MSY, and status probabilities.
#' @details Fleet harvest fractions are `u * fleet_weights[f] * selectivity[f,a]`.
#'   Natural mortality follows seasonal fishing, matching historical dynamics.
#'   Yield is total catch in the units of `weight_fya_mod`, including fleets
#'   whose input catches use numbers. Beverton-Holt equilibrium recruitment is
#'   solved analytically without recruitment deviations or bias corrections.
#'   No positive continuation is accepted as a biological equilibrium.
#'
#'   These are deterministic, fixed-biology reference points, not stochastic
#'   MSY or management advice. Spawning potential at recruitment age must be
#'   zero to match the dynamics' spawning-before-recruitment convention. An
#'   upper-bound optimum is flagged and should not be treated as a resolved
#'   MSY. Posterior calculations use each draw's own biology and do not discard
#'   invalid draws. Equal draw weights are used in status probabilities.
#' @family assessment reference points
#' @export
opal_msy <- function(x, fleet_weights, year = x$data$n_year,
                     uncertainty = c("fit", "mcmc"), draws = NULL,
                     u_max = 1, grid_size = 201L, name = "msy") {
  uncertainty <- match.arg(uncertainty)
  validate_opal_obj(x)
  d <- x$data
  if (!is.numeric(fleet_weights) || length(fleet_weights) != d$n_fishery ||
      any(!is.finite(fleet_weights)) || any(fleet_weights < 0) || sum(fleet_weights) <= 0) stop("Supply non-negative fleet weights with a positive sum.")
  fleet_weights <- fleet_weights / sum(fleet_weights)
  if (length(year) != 1L || !is.finite(year) || year != floor(year) || year < 1L || year > d$n_year) stop("Invalid reference year index.")
  if (length(u_max) != 1L || !is.finite(u_max) || u_max <= 0 || u_max > 1) stop("`u_max` must be in (0, 1].")
  if (length(grid_size) != 1L || !is.finite(grid_size) || grid_size < 11 || grid_size != floor(grid_size)) stop("Use at least 11 grid points.")
  if (uncertainty == "mcmc") {
    context <- .opal_posterior_context(x, draws)
  } else {
    if (is.null(x$fit$opt)) stop("Point reference points require a fitted model.")
    if (!is.null(draws)) stop("`draws` applies only to MCMC reference points.")
    object <- opal_rtmb(x, fresh = TRUE)
    context <- list(object = object, values = matrix(object$env$last.par.best, nrow = 1L),
      ids = data.frame(iteration = NA_integer_, chain = NA_integer_, draw = 1L))
  }
  results <- vector("list", nrow(context$values))
  for (i in seq_len(nrow(context$values))) {
    report <- context$object$report(context$values[i, ])
    if (!.opal_diagnose_report(d, context$object$env$parList(par = context$values[i, ]), report)$passes) {
      stop("Biologically invalid source at draw ", context$ids$draw[i], "; reference points were not stored.")
    }
    if (report$spawning_potential_a[1L] != 0) stop("MSY requires zero spawning potential at recruitment age.")
    sel <- matrix(report$sel_fya[, year, ], d$n_fishery, d$n_age)
    weight <- matrix(report$weight_fya_mod[, year, ], d$n_fishery, d$n_age)
    equilibrium <- function(u) .opal_equilibrium(u, fleet_weights, sel, weight,
      report$M_a, report$spawning_potential_a, report$alpha, report$beta, d$n_season)
    grid <- seq(0, u_max, length.out = grid_size)
    yield <- vapply(grid, function(u) equilibrium(u)$yield, numeric(1))
    peaks <- which(yield > 0 & yield >= c(-Inf, utils::head(yield, -1L)) &
                     yield >= c(utils::tail(yield, -1L), -Inf))
    candidates <- c(0, u_max, vapply(peaks, function(j) {
      stats::optimise(function(u) -equilibrium(u)$yield,
        c(grid[max(1L, j - 1L)], grid[min(grid_size, j + 1L)]), tol = 1e-9)$minimum
    }, numeric(1)))
    u <- candidates[which.max(vapply(candidates, function(u) equilibrium(u)$yield, numeric(1)))]
    state <- equilibrium(u)
    terminal <- utils::tail(report$spawning_biomass_y, 1L)
    results[[i]] <- cbind(context$ids[i, ], data.frame(msy = state$yield, b_msy = state$spawning,
      r_msy = state$recruitment, u_msy = u, depletion_msy = state$spawning / report$B0,
      terminal_b_bmsy = if (state$spawning > 0) terminal / state$spawning else NA_real_,
      at_upper_bound = abs(u - u_max) <= 1e-7, resolved = state$yield > 0 && state$spawning > 0 && abs(u-u_max) > 1e-7))
  }
  table <- do.call(rbind, results)
  result <- list(draws = table, probability_below_bmsy = if (all(table$resolved) &&
    all(is.finite(table$terminal_b_bmsy))) mean(table$terminal_b_bmsy < 1) else NA_real_,
    fleet_weights = fleet_weights, year = year, uncertainty = uncertainty,
    validation = summary(x)$checks[[uncertainty]],
    units = "Yield: model weight units; spawning output: model spawning units")
  .opal_store_derived(x, name, result, uncertainty,
    list(fleet_weights = fleet_weights, year = year, draws = draws, u_max = u_max, grid_size = grid_size))
}

.opal_equilibrium <- function(u, fleet_weights, sel, weight, mortality, spawning, alpha, beta, seasons) {
  harvest <- sweep(sel, 1L, u * fleet_weights, "*")
  seasonal_survival <- (1 - colSums(harvest)) * exp(-mortality / seasons)
  if (any(seasonal_survival < 0 | seasonal_survival > 1)) stop("Invalid equilibrium survival.")
  survival <- seasonal_survival^seasons
  n <- length(mortality)
  per_recruit <- c(1, if (n > 1) cumprod(survival[-n]) else numeric())
  per_recruit[n] <- per_recruit[n] / (1 - survival[n])
  spr <- sum(per_recruit * spawning)
  recruitment <- if (spr > 0) max(0, alpha - beta / spr) else 0
  abundance <- per_recruit * recruitment
  output <- sum(abundance * spawning)
  yield <- 0
  for (s in seq_len(seasons)) {
    yield <- yield + sum(sweep(harvest * weight, 2L, abundance, "*"))
    abundance <- abundance * seasonal_survival
  }
  list(yield = yield, spawning = output, recruitment = recruitment,
       spr = spr, numbers = per_recruit * recruitment)
}
