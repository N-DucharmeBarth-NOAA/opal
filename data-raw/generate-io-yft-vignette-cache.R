# Precompute what the Indian Ocean yellowfin articles plot, so that they
# render in seconds without rebuilding the RTMB model:
#
#   inst/extdata/io_yft_vignette_cache.rds
#
# Read by vignettes/io-yft-ss3-comparison.qmd and
# vignettes/io-yft-age-selectivity.qmd. Run from the package root after
# generate-io-yft-c2-fit.R or generate-io-yft-age-fits.R, or after any change
# to opal's model code (about a minute: one model build per entry).
pkgload::load_all()
source("data-raw/io-yft-age-setup.R")

d <- io_yft_data
n_ts <- d$n_year
read_fit <- function(file) opal_read(file.path("inst", "extdata", file))

# Everything the articles take from one model: report quantities as plain
# data frames, and the fit statistics when the model has been fitted.
extract <- function(x) {
  object <- opal_rtmb(x, fresh = TRUE)
  par <- if (is.null(x$fit$opt)) object$par else x$fit$opt$par
  report <- object$report(par)
  sel <- report$sel_fya[, 1, ]
  lf <- do.call(rbind, lapply(seq_along(d$lf_fishery_f), function(j) {
    n <- d$lf_n_fi[[j]]
    data.frame(fishery = d$lf_fishery_f[j], bin = d$len_lower,
               predicted = colSums(report$lf_pred[[j]] * n) / sum(n))
  }))
  traj <- data.frame(ts = seq_len(n_ts),
                     ssb = report$spawning_biomass_y[seq_len(n_ts)],
                     static_dep = report$static_depletion_y[seq_len(n_ts)])
  fit <- if (is.null(x$fit$opt)) NULL else list(
    n_par = length(x$fit$opt$par),
    objective = x$fit$opt$objective,
    max_gradient = x$validation$fit$metrics$max_gradient,
    passes = isTRUE(x$validation$fit$passes),
    summary = utils::capture.output(summary(x)))
  list(
    B0 = report$B0,
    objective = as.numeric(object$fn(par)),
    nll = c(index = sum(report$lp_cpue), length = sum(report$lp_lf),
            recruitment = sum(report$lp_rec), prior = sum(report$lp_prior),
            penalty = sum(report$lp_penalty)),
    trajectory = traj,
    sel_age = data.frame(fishery = rep(seq_len(nrow(sel)), times = ncol(sel)),
                         age = rep(d$age_a, each = nrow(sel)), sel = as.vector(sel)),
    cpue = data.frame(ts = d$cpue_data$ts, fishery = d$cpue_data$fishery,
                      pred = report$cpue_pred),
    lf = lf,
    rdev = unname(par[grepl("^rdev_y", names(par))]),
    sel_free = sum(grepl("^par_sel", names(par))),
    fit = fit)
}

inputs <- io_yft_age_inputs()
sa <- io_yft_reference$ss3$sel_age
external_data <- inputs$data
external_data$sel_fa_external <- matrix(0, d$n_fishery, d$n_age)
external_data$sel_fa_external[cbind(sa$fishery, match(sa$age, d$age_a))] <- sa$sel
models <- list(
  C2 = read_fit("io_yft_c2_fit.rds"),
  A = read_fit("io_yft_age_fit.rds"),
  A_ss3 = read_fit("io_yft_age_ss3_fit.rds"),
  # Age-based selectivity at SS3's values, everything else at SS3's
  # solution, not fitted: the bridging check.
  A_start = opal_obj(inputs$data, inputs$parameters, inputs$map_c2,
                     bounds = io_yft_age_bounds(inputs$parameters, inputs$map_c2)),
  # The same, with SS3's selectivity-at-age supplied as a fixed matrix
  E_start = opal_obj(external_data, inputs$parameters, inputs$map_c2,
                     bounds = io_yft_age_bounds(inputs$parameters, inputs$map_c2)))

cache <- list(
  created = format(Sys.time(), tz = "UTC", usetz = TRUE),
  opal = as.character(utils::packageVersion("opal")),
  commit = tryCatch(system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE),
                    error = function(e) NA_character_),
  models = lapply(models, extract))

saveRDS(cache, file.path("inst", "extdata", "io_yft_vignette_cache.rds"),
        compress = "xz", version = 3L)
for (nm in names(cache$models)) {
  m <- cache$models[[nm]]
  cat(sprintf("%-8s objective %.4f, B0 %.4g, terminal depletion %.3f\n", nm,
              m$objective, m$B0, m$trajectory$static_dep[n_ts]))
}
