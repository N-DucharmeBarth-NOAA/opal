.opal_with_seed <- function(seed, code) {
  if (!is.null(seed)) {
    had <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
    old <- if (had) get(".Random.seed", .GlobalEnv) else NULL
    on.exit(if (had) assign(".Random.seed", old, .GlobalEnv) else
      if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv))
    set.seed(seed)
  }
  force(code)
}

# The sequential conditional distributions are binomial, beta, and
# beta-binomial. Their product is the corresponding composition likelihood.
.opal_composition_pit <- function(observed, predicted, family, concentration = NULL) {
  k <- length(predicted)
  out <- rep(NA_real_, k)
  if (k < 2L) return(out)
  if (family == 1L) observed <- floor(observed + 0.5)
  if (family == 2L && (any(observed <= 0) || abs(sum(observed) - 1) > 1e-6)) {
    stop("Dirichlet observations must be positive proportions summing to one.")
  }
  alpha <- predicted * concentration
  logsum <- function(z) {
    if (!length(z) || all(z == -Inf)) return(-Inf)
    m <- max(z)
    m + log(sum(exp(z - m)))
  }
  for (i in seq_len(k - 1L)) {
    rest <- seq.int(i, k)
    total <- sum(observed[rest])
    if (total <= 0) next
    if (family == 2L) {
      fraction <- observed[i] / total
      a <- alpha[i]; b <- sum(alpha[rest[-1L]])
      left <- stats::pbeta(fraction, a, b, log.p = TRUE)
      right <- stats::pbeta(fraction, a, b, lower.tail = FALSE, log.p = TRUE)
      out[i] <- if (left < log(0.5)) stats::qnorm(left, log.p = TRUE) else
        stats::qnorm(right, lower.tail = FALSE, log.p = TRUE)
    } else {
      q <- observed[i]
      if (family == 1L) {
        prob <- predicted[i] / sum(predicted[rest])
        lo <- stats::pbinom(q - 1, total, prob, log.p = TRUE)
        hi <- stats::pbinom(q, total, prob, lower.tail = FALSE, log.p = TRUE)
        mass <- stats::dbinom(q, total, prob, log = TRUE)
      } else {
        a <- alpha[i]; b <- sum(alpha[rest[-1L]])
        logmass <- function(z) lchoose(total, z) + lbeta(z + a, total - z + b) - lbeta(a, b)
        mass <- logmass(q)
        lo <- if (q == 0) -Inf else logsum(logmass(seq.int(0, q - 1)))
        hi <- if (q == total) -Inf else logsum(logmass(seq.int(q + 1, total)))
      }
      u <- stats::runif(1)
      left <- logsum(c(lo, log(u) + mass))
      right <- logsum(c(hi, log1p(-u) + mass))
      out[i] <- if (left < log(0.5)) stats::qnorm(left, log.p = TRUE) else
        stats::qnorm(right, lower.tail = FALSE, log.p = TRUE)
    }
  }
  out
}

.opal_composition_rows <- function(data, parameters, report, type) {
  field <- function(s) data[[paste0(type, "_", s)]]
  family <- field("switch")
  n_comp <- data[[paste0("n_", type)]]
  if (is.null(family) || family == 0L || is.null(n_comp) || n_comp == 0L) return(list())
  predictions <- report[[paste0(type, "_pred")]]
  if (!length(predictions)) return(list())
  years <- field("year_fi")
  if (is.null(years)) years <- split(field("year"), field("fishery"))
  sizes <- field("n_fi")
  if (is.null(sizes)) sizes <- split(field("n"), field("fishery"))
  ints <- field("n_int_fi")
  observed <- field(c("obs_flat", "obs_prop", "obs_ints")[family])
  offset <- row_id <- 0L
  out <- list()
  for (j in seq_along(predictions)) {
    f <- field("fishery_f")[j]
    bins <- seq.int(field("minbin")[f], field("maxbin")[f])
    for (i in seq_len(nrow(predictions[[j]]))) {
      row_id <- row_id + 1L
      n <- sizes[[j]][i]
      n_int <- if (is.null(ints)) floor(n + 0.5) else ints[[j]][i]
      include <- data$removal_switch_f[f] == 0 && if (family == 3L) n_int > 0 else n > 0
      obs <- observed[offset + seq_along(bins)]
      offset <- offset + length(bins)
      concentration <- exp(parameters[[paste0("log_", type, "_tau")]][f])
      if (family == 2L) concentration <- concentration * n
      year_index <- years[[j]][i]
      year <- if (length(data$years) >= year_index) data$years[year_index] else data$first_yr + year_index - 1L
      out[[row_id]] <- list(type = type, family = family, row = row_id, fishery = f,
        year = year, bin = bins, observed = obs, predicted = predictions[[j]][i, ],
        concentration = concentration, include = include,
        observation_index = offset - length(bins) + seq_along(bins))
    }
  }
  out
}

.opal_sdnr <- function(z, conf) {
  z <- z[is.finite(z)]
  n <- length(z)
  s <- if (n >= 2L) stats::sd(z) else NA_real_
  interval <- if (n >= 2L) sqrt((n - 1) * s^2 /
    stats::qchisq(c((1 + conf) / 2, (1 - conf) / 2), n - 1)) else c(NA_real_, NA_real_)
  data.frame(n = n, sdnr = s, lower = interval[1], upper = interval[2],
             mean = if (n) mean(z) else NA_real_, mar = if (n) stats::median(abs(z)) else NA_real_)
}

#' Calculate one-step-ahead observation residuals
#'
#' Adds CPUE, length-composition, and weight-composition OSA diagnostics to an
#' Opal object. Use [plot_osa_sdnr()] for the dataset-level SDNR comparison.
#' @param x A fitted `opal_obj`.
#' @param seed Random seed for discrete residuals; restored on exit. `NULL`
#'   uses and advances the current random-number stream.
#' @param conf Confidence level for approximate chi-squared SDNR intervals.
#' @param name Name under `x$derived` for the diagnostics.
#' @return An updated object. [opal_derived()] returns its residual table,
#'   dataset summaries, and calculation settings.
#' @details For fixed-effect fits, CPUE uses its lognormal CDF, and composition
#'   residuals use exact sequential binomial, beta, or beta-binomial CDFs for
#'   multinomial, Dirichlet, or Dirichlet-multinomial observations. Parameters
#'   are held at their fitted values. For random-effect models, RTMB's
#'   `oneStepPredict()` integrates latent states using its Laplace-based OSA
#'   approximation. Other observation streams remain conditioned on. There is
#'   no silent substitution of conditional residuals for marginal residuals.
#'
#'   Observation order is the stored input order, then ascending composition
#'   bin. The final composition bin is constrained and omitted. Removed
#'   fisheries, zero sample sizes, exhausted counts, and numerical failures
#'   remain in the table with explicit exclusion reasons. Multinomial counts
#'   are rounded as in the fitted likelihood; no new pseudocounts are added.
#'   SDNR should be near one, but its interval is approximate because fitted
#'   parameters and finite samples affect calibration. Catch is conditioned on
#'   and has no observation likelihood; process priors are not observations.
#' @family assessment diagnostics
#' @export
opal_osa <- function(x, seed = 123L, conf = 0.95, name = "osa") {
  validate_opal_obj(x)
  if (is.null(x$fit$opt)) stop("OSA diagnostics require a fitted model.", call. = FALSE)
  if (length(conf) != 1 || !is.finite(conf) || conf <= 0 || conf >= 1) stop("`conf` must be between zero and one.")
  result <- .opal_with_seed(seed, .opal_residual_result(x, conf))
  .opal_store_derived(x, name, result, "fit", list(seed = seed, conf = conf))
}

.opal_residual_result <- function(x, conf) {
  obj <- opal_rtmb(x, fresh = TRUE)
  r <- obj$report(obj$env$last.par.best)
  p <- obj$env$parList(par = obj$env$last.par.best)
  d <- x$data
  latent <- length(x$random) > 0L
  source <- if (latent) "RTMB marginal OSA" else "exact conditional CDF"
  pieces <- list()
  if (isTRUE(d$cpue_switch > 0L) && nrow(d$cpue_data)) {
    cpue <- d$cpue_data
    idx <- if (is.null(cpue$index)) rep(1L, nrow(cpue)) else cpue$index
    residual <- (log(cpue$value) - log(r$cpue_pred)) / r$cpue_sigma
    if (latent) residual <- RTMB::oneStepPredict(obj, observation.name = "cpue_log_obs",
      method = "oneStepGeneric", discrete = FALSE, seed = NULL, trace = FALSE)$residual
    pieces[[1L]] <- data.frame(dataset = paste("CPUE", idx), type = "cpue", family = "lognormal",
      fishery = cpue$fishery, row = seq_len(nrow(cpue)), bin = NA_integer_,
      year = d$first_yr + (cpue$ts - 1L) %/% d$n_season,
      observation_index = seq_len(nrow(cpue)), observed = cpue$value, predicted = r$cpue_pred,
      residual = residual, reason = ifelse(is.finite(residual), "", "non-finite residual"), source = source)
  }
  for (type in c("lf", "wf")) {
    rows <- .opal_composition_rows(d, p, r, type)
    if (!length(rows)) next
    flat <- lapply(rows, function(row) {
      k <- length(row$bin)
      z <- rep(NA_real_, k)
      reason <- rep(if (row$include) "" else "excluded by likelihood", k)
      if (row$include) {
        reason[k] <- "constrained final bin"
        if (!latent) z <- .opal_composition_pit(row$observed, row$predicted, row$family, row$concentration)
        if (row$family != 2L && k > 1L) {
          counts <- if (row$family == 1L) floor(row$observed + 0.5) else row$observed
          exhausted <- rev(cumsum(rev(counts))) == 0
          reason[exhausted & seq_len(k) < k] <- "exhausted counts"
        }
      }
      data.frame(dataset = paste(toupper(type), row$fishery), type = type,
        family = c("multinomial", "Dirichlet", "Dirichlet-multinomial")[row$family],
        fishery = row$fishery, row = row$row, bin = row$bin, year = row$year,
        observation_index = row$observation_index, observed = row$observed,
        predicted = row$predicted, residual = z, reason = reason, source = source)
    })
    df <- do.call(rbind, flat)
    if (latent && any(df$reason == "")) {
      family <- rows[[1L]]$family
      observation <- paste0(type, "_", c("obs_flat", "obs_prop", "obs_ints")[family])
      select <- which(df$reason == "")
      args <- list(obj = opal_rtmb(x, fresh = TRUE), observation.name = observation,
        method = "oneStepGeneric", discrete = family != 2L,
        subset = df$observation_index[select], seed = NULL, trace = FALSE)
      args$range <- if (family == 2L) c(0, 1) else c(0, Inf)
      residual <- do.call(RTMB::oneStepPredict, args)$residual
      if (length(residual) == nrow(df)) residual <- residual[select]
      if (length(residual) != length(select)) stop("RTMB returned an incompatible OSA layout.")
      df$residual[select] <- residual
    }
    df$reason[df$reason == "" & !is.finite(df$residual)] <- "non-finite residual"
    df$residual[df$reason != ""] <- NA_real_
    pieces[[length(pieces) + 1L]] <- df
  }
  if (!length(pieces)) stop("No active observation likelihoods were found.", call. = FALSE)
  residuals <- do.call(rbind, pieces)
  groups <- split(residuals, factor(residuals$dataset, levels = unique(residuals$dataset)))
  summary <- do.call(rbind, lapply(groups, function(z) cbind(
    data.frame(dataset = z$dataset[1], family = z$family[1], source = z$source[1],
      excluded = sum(z$reason != ""), failures = sum(z$reason == "non-finite residual")),
    .opal_sdnr(z$residual, conf))))
  rownames(summary) <- NULL
  if (any(summary$failures > 0)) warning("Some OSA residuals are non-finite; inspect the exclusion reasons.", call. = FALSE)
  list(residuals = residuals, summary = summary, conf = conf,
       conditioning = "Fitted parameters; other observation streams held observed")
}

#' Compare OSA residual dispersion across datasets
#' @param x An Opal object containing stored OSA diagnostics.
#' @param name Name of the stored diagnostic result.
#' @return A ggplot with one dataset per row, SDNR on the horizontal axis,
#'   approximate intervals, and a reference line at one. Plot data contain all
#'   dataset summaries, including counts of excluded and failed residuals.
#' @family assessment diagnostics
#' @export
plot_osa_sdnr <- function(x, name = "osa") {
  result <- opal_derived(x, name)
  df <- result$summary
  df$dataset <- factor(df$dataset, levels = rev(unique(df$dataset)))
  ggplot2::ggplot(df, ggplot2::aes(x = .data$sdnr, y = .data$dataset)) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", colour = "grey45") +
    ggplot2::geom_segment(ggplot2::aes(x = .data$lower, xend = .data$upper, yend = .data$dataset)) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::labs(x = "SDNR of OSA residuals", y = "Dataset",
      caption = paste0(round(100 * result$conf), "% approximate intervals; reference SDNR = 1")) +
    ggplot2::theme_bw()
}

#' Plot stored OSA residuals
#' @param x An Opal object containing stored OSA diagnostics.
#' @param name Name of the stored diagnostic result.
#' @param type Plot residuals against model year or as normal quantiles.
#' @return A ggplot, faceted by dataset.
#' @family assessment diagnostics
#' @export
plot_osa_residuals <- function(x, name = "osa", type = c("time", "qq")) {
  type <- match.arg(type)
  df <- opal_derived(x, name)$residuals
  df <- df[is.finite(df$residual), ]
  if (type == "qq") {
    p <- ggplot2::ggplot(df, ggplot2::aes(sample = .data$residual)) +
      ggplot2::stat_qq() + ggplot2::stat_qq_line() +
      ggplot2::labs(x = "Standard normal quantile", y = "OSA residual")
  } else {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$year, y = .data$residual)) +
      ggplot2::geom_hline(yintercept = 0, colour = "grey60") + ggplot2::geom_point() +
      ggplot2::labs(x = "Model year", y = "OSA residual")
  }
  p + ggplot2::facet_wrap(~dataset) + ggplot2::theme_bw()
}
