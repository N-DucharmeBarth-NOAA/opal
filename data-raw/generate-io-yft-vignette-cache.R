# Precompute what the Indian Ocean yellowfin article plots, so that it
# renders in seconds without rebuilding the RTMB model:
#
#   inst/extdata/io_yft_vignette_cache.rds
#
# Read by vignettes/io-yft-ss3-comparison.qmd. Run from the package root
# after generate-io-yft-fits.R, or after any change to opal's model code
# (about a minute: one model build per entry).
pkgload::load_all()
source("data-raw/io-yft-setup.R")

d <- io_yft_data
n_ts <- d$n_year
read_fit <- function(file) opal_read(file.path("inst", "extdata", file))

# Everything the articles take from one model: report quantities as plain
# data frames, and the fit statistics when the model has been fitted. Length
# compositions on the model's 5 cm grid go in `lf`; the index size
# compositions on SS3's sizefreq bins, when the model has them, in
# `sizefreq`.
extract <- function(x, par = NULL) {
  object <- opal_rtmb(x, fresh = TRUE)
  if (is.null(par)) par <- if (is.null(x$fit$opt)) object$par else x$fit$opt$par
  report <- object$report(par)
  sel <- report$sel_fya[, 1, ]
  xd <- x$data
  native <- vapply(report$lf_pred, ncol, integer(1)) == xd$n_len
  lf <- do.call(rbind, lapply(which(native), function(j) {
    n <- xd$lf_n_fi[[j]]
    data.frame(fishery = xd$lf_fishery_f[j], bin = xd$len_lower,
               predicted = colSums(report$lf_pred[[j]] * n) / sum(n))
  }))
  sf_bins <- io_yft_reference$ss3$sizefreq$bins
  sizefreq <- if (any(!native)) do.call(rbind, lapply(which(!native), function(j) {
    p <- report$lf_pred[[j]]
    p <- p / rowSums(p)
    data.frame(fishery = xd$lf_fishery_f[j],
               ts = rep(xd$lf_year_fi[[j]], each = ncol(p)),
               bin = rep(sf_bins, nrow(p)), predicted = as.vector(t(p)))
  })) else NULL
  lf_f <- rep(xd$lf_fishery_f, xd$lf_n_f)
  lf_native <- rep(native, xd$lf_n_f)
  traj <- data.frame(ts = seq_len(n_ts),
                     ssb = report$spawning_biomass_y[seq_len(n_ts)],
                     static_dep = report$static_depletion_y[seq_len(n_ts)])
  fit <- if (is.null(x$fit$opt) || !identical(unname(par), unname(x$fit$opt$par))) NULL else list(
    n_par = length(x$fit$opt$par),
    objective = x$fit$opt$objective,
    max_gradient = x$validation$fit$metrics$max_gradient,
    passes = isTRUE(x$validation$fit$passes),
    positive_hessian = isTRUE(x$validation$fit$metrics$positive_hessian),
    summary = utils::capture.output(summary(x)))
  list(
    B0 = report$B0,
    objective = as.numeric(object$fn(par)),
    nll = c(index = sum(report$lp_cpue), length = sum(report$lp_lf[lf_native]),
            recruitment = sum(report$lp_rec), prior = sum(report$lp_prior),
            penalty = sum(report$lp_penalty)),
    # lambda-weighted terms, as they enter the objective
    nll_w = c(index = sum(report$lp_cpue_w), length = sum(report$lp_lf_w[lf_native]),
              sizefreq = sum(report$lp_lf_w[!lf_native]),
              recruitment = sum(report$lp_rec_w), prior = sum(report$lp_prior_w),
              penalty = sum(report$lp_penalty_w)),
    sizefreq_nll = sum(report$lp_lf[!lf_native]),
    lf_nll_f = tapply(report$lp_lf, lf_f, sum),
    trajectory = traj,
    sel_age = data.frame(fishery = rep(seq_len(nrow(sel)), times = ncol(sel)),
                         age = rep(d$age_a, each = nrow(sel)), sel = as.vector(sel)),
    cpue = data.frame(ts = d$cpue_data$ts, fishery = d$cpue_data$fishery,
                      pred = report$cpue_pred),
    lf = lf,
    sizefreq = sizefreq,
    rdev = unname(par[grepl("^rdev_y", names(par))]),
    sel_free = sum(grepl("^par_sel", names(par))),
    fit = fit)
}

age_sf <- io_yft_inputs()$`AGE-SF`
models <- list(
  # The article's model
  age_sf = read_fit("io_yft_age_sf_fit.rds"),
  # AGE-SF's inputs at SS3's solution, not fitted: opal's predictions of the
  # index size compositions against SS3's expected values
  age_sf_start = opal_obj(age_sf$data, age_sf$parameters, age_sf$map,
                          bounds = io_yft_age_bounds(age_sf$parameters, age_sf$map)))
stopifnot(isTRUE(all.equal(models$age_sf$data, age_sf$data)))

cache <- list(
  created = format(Sys.time(), tz = "UTC", usetz = TRUE),
  opal = as.character(utils::packageVersion("opal")),
  commit = tryCatch(system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE),
                    error = function(e) NA_character_),
  models = lapply(models, extract),
  # The article's bridge table (generate-io-yft-fits.R)
  bridge = readRDS(file.path("inst", "extdata", "io_yft_bridge.rds")),
  # What the article checks its inline AGE-SF inputs against
  inputs_fingerprint = io_yft_fingerprint(age_sf$data, age_sf$parameters, age_sf$map))

saveRDS(cache, file.path("inst", "extdata", "io_yft_vignette_cache.rds"),
        compress = "xz", version = 3L)
for (nm in names(cache$models)) {
  m <- cache$models[[nm]]
  cat(sprintf("%-8s objective %.4f, B0 %.4g, terminal depletion %.3f\n", nm,
              m$objective, m$B0, m$trajectory$static_dep[n_ts]))
}
