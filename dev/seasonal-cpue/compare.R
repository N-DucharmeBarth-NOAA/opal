# Quantify the CPUE indexing correction, without changing the package model.
# Run from the repository root: Rscript --vanilla dev/seasonal-cpue/compare.R
# Output files are read by vignettes/seasonal-cpue-impact.qmd.
pkgload::load_all('.', quiet = TRUE)
output <- 'inst/extdata/seasonal-cpue-impact'
dir.create(output, recursive = TRUE, showWarnings = FALSE)

legacy_environment <- new.env(parent = asNamespace('opal'))
sys.source('dev/seasonal-cpue/legacy-cpue-v4.R', legacy_environment)
legacy_model <- opal_model
environment(legacy_model) <- list2env(
  list(get_cpue_like = legacy_environment$legacy_cpue_like_v4),
  parent = environment(opal_model))
models <- list(v4_indexing = legacy_model, v5_indexing = opal_model)
make_object <- function(inputs, model) {
  RTMB::MakeADFun(cmb(model, inputs$data), inputs$parameters,
                 map = inputs$map, silent = TRUE)
}

make_inputs <- function(seasons, full_coverage = FALSE) {
  inputs <- opal_example_inputs(seed = 712L)
  d <- inputs$data
  d$n_season <- seasons
  # Preserve annual catches across scenarios; only their seasonal allocation changes.
  seasonal_share <- seq_len(seasons) / sum(seq_len(seasons))
  annual_catch <- 2 * seq(30, 70, length.out = d$n_year)
  d$catch_obs_ysf <- array(rep(outer(annual_catch, seasonal_share), d$n_fishery),
                           c(d$n_year, seasons, d$n_fishery))
  # The legacy code cannot evaluate ts > n_year. Retain exactly the same early
  # observations in both paired fits, and report that restricted coverage.
  times <- seq_len(if (full_coverage) d$n_year * seasons else d$n_year)
  d$cpue_data <- data.frame(ts = rep(times, d$n_index),
    fishery = rep(seq_len(d$n_index), each = length(times)),
    index = rep(seq_len(d$n_index), each = length(times)), units = 1L,
    value = 1, se = .1)
  inputs$data <- d
  object <- make_object(inputs, opal_model)
  report <- object$report(object$par)
  # Paired runs use identical simulated observations and the same starting point.
  # These are controlled examples, not real-stock estimates or a calibration study.
  set.seed(90210L)
  inputs$data$cpue_data$value <- exp(stats::rnorm(nrow(d$cpue_data),
    log(report$cpue_pred), report$cpue_sigma))
  for (type in c('lf', 'wf')) {
    obs <- unlist(lapply(report[[paste0(type, '_pred')]], function(pred) {
      unlist(lapply(seq_len(nrow(pred)), function(i) {
        as.numeric(stats::rmultinom(1, 120L, pred[i, ]))
      }))
    }), use.names = FALSE)
    inputs$data[[paste0(type, '_obs_flat')]] <- obs
    inputs$data[[paste0(type, '_obs_ints')]] <- as.integer(obs)
    inputs$data[[paste0(type, '_obs_prop')]] <- obs / 120
  }
  inputs
}

fit_object <- function(object, inputs) {
  bounds <- get_bounds(object, inputs$parameters)
  start <- object$par
  for (pass in 1:2) {
    opt <- stats::nlminb(start, object$fn, object$gr, lower = bounds$lower,
      upper = bounds$upper, control = list(iter.max = 10000, eval.max = 10000,
                                         rel.tol = 1e-10, x.tol = 1e-8))
    start <- opt$par
  }
  gradient <- max(abs(object$gr(opt$par)))
  hessian_min <- min(eigen(object$he(opt$par), symmetric = TRUE)$values)
  report <- object$report(opt$par)
  biology <- opal:::.opal_diagnose_report(inputs$data,
    object$env$parList(par = opt$par), report)
  if (opt$convergence != 0L || gradient >= 1e-3 || hessian_min <= 0 || !biology$passes) {
    print(list(seasons = inputs$data$n_season, opt = opt, gradient = gradient,
               hessian_min = hessian_min, biology = biology))
  }
  stopifnot(opt$convergence == 0L, gradient < 1e-3, hessian_min > 0, biology$passes)
  list(opt = opt, report = report, gradient = gradient,
       hessian_min = hessian_min, biology = biology$passes)
}

fixed <- fitted <- parameters <- predictions <- failures <- list()
for (seasons in c(1L, 2L, 4L)) {
  inputs <- make_inputs(seasons)
  objects <- lapply(models, function(model) make_object(inputs, model))
  reports <- lapply(objects, function(o) o$report(o$par))
  nll <- vapply(objects, function(o) as.numeric(o$fn(o$par)), numeric(1))
  # Verify that the paired implementation changes only the CPUE contribution.
  cpue_fields <- c('cpue_pred', 'cpue_res', 'cpue_sigma', 'lp_cpue')
  unchanged <- setdiff(names(reports[[1L]]), cpue_fields)
  stopifnot(identical(reports[[1L]][unchanged], reports[[2L]][unchanged]),
    isTRUE(all.equal(unname(nll[2L] - nll[1L]),
      sum(reports[[2L]]$lp_cpue) - sum(reports[[1L]]$lp_cpue), tolerance = 1e-10)))
  old <- reports$v4_indexing$cpue_pred
  new <- reports$v5_indexing$cpue_pred
  change <- 100 * (new / old - 1)
  fixed[[length(fixed) + 1L]] <- data.frame(seasons = seasons,
    years = inputs$data$n_year, observations = length(old),
    last_observed_year = ceiling(inputs$data$n_year / seasons),
    median_absolute_prediction_change_pct = median(abs(change)),
    max_absolute_prediction_change_pct = max(abs(change)),
    nll_v4 = nll[1], nll_v5 = nll[2], nll_change = nll[2] - nll[1],
    cpue_nll_change = sum(reports$v5_indexing$lp_cpue) - sum(reports$v4_indexing$lp_cpue))
  fits <- lapply(objects, fit_object, inputs = inputs)
  values <- lapply(fits, function(f) c(
    B0 = f$report$B0, q1 = unname(exp(f$opt$par[2L])), q2 = unname(exp(f$opt$par[3L])),
    terminal_spawning_output = tail(f$report$spawning_biomass_y, 1L),
    terminal_depletion = tail(f$report$static_depletion_y, 1L)))
  parameters[[length(parameters) + 1L]] <- data.frame(seasons = seasons,
    quantity = names(values[[1L]]), v4 = unname(values[[1L]]),
    v5 = unname(values[[2L]]), change_pct = 100 * (values[[2L]] / values[[1L]] - 1))
  for (method in names(fits)) {
    f <- fits[[method]]
    fitted[[length(fitted) + 1L]] <- data.frame(seasons = seasons, method = method,
      objective = f$opt$objective, convergence = f$opt$convergence,
      max_gradient = f$gradient, minimum_hessian_eigenvalue = f$hessian_min,
      biology_passes = f$biology)
  }
  predictions[[length(predictions) + 1L]] <- cbind(seasons = seasons,
    inputs$data$cpue_data, old_prediction_at_start = old, new_prediction_at_start = new,
    old_fitted_prediction = fits$v4_indexing$report$cpue_pred,
    new_fitted_prediction = fits$v5_indexing$report$cpue_pred)
  if (seasons == 1L) {
    stopifnot(identical(nll[[1L]], nll[[2L]]),
      identical(objects[[1L]]$gr(objects[[1L]]$par), objects[[2L]]$gr(objects[[2L]]$par)),
      identical(reports[[1L]], reports[[2L]]), identical(values[[1L]], values[[2L]]))
  } else {
    full <- make_inputs(seasons, full_coverage = TRUE)
    old_error <- tryCatch({ make_object(full, legacy_model); NA_character_ },
                          error = conditionMessage)
    current <- make_object(full, opal_model)
    stopifnot(!is.na(old_error), is.finite(current$fn(current$par)))
    failures[[length(failures) + 1L]] <- data.frame(seasons = seasons,
      observations = nrow(full$data$cpue_data), legacy_error = old_error,
      current_objective_finite = TRUE)
  }
}
save_table <- function(x, file) utils::write.csv(do.call(rbind, x),
  file.path(output, file), row.names = FALSE)
save_table(fixed, 'fixed-parameter-comparison.csv')
save_table(fitted, 'fit-diagnostics.csv')
save_table(parameters, 'fitted-quantities.csv')
save_table(predictions, 'predictions.csv')
save_table(failures, 'full-coverage.csv')

# Check real bundled input conventions independently of the simulated fits.
inventory <- lapply(c('opal_baseline', 'opakapaka', 'wcpo_bet'), function(model) {
  d <- get_data(model)
  data.frame(model = model, n_year = d$n_year, n_season = d$n_season,
             min_ts = min(d$cpue_data$ts), max_ts = max(d$cpue_data$ts))
})
save_table(inventory, 'bundled-inputs.csv')
stopifnot(all(vapply(inventory, function(d) d$n_season == 1L, logical(1))))
writeLines(c('Legacy implementation: 2dd8e3de2bbfb4f00a9489c1ae456d8dd2ce5e14',
             'Current scientific contract: opal_model_contract_v5',
             sub('[[:space:]]+$', '', capture.output(sessionInfo()))),
           file.path(output, 'session-info.txt'))
print(do.call(rbind, fixed))
print(do.call(rbind, parameters))
print(do.call(rbind, fitted))
