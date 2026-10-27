# Fit the Indian Ocean yellowfin models of the io-yft-ss3-comparison
# article (see io-yft-setup.R for what each one is). Run from the package
# root; takes about 45 minutes. Afterwards, rerun
# generate-io-yft-vignette-cache.R.
#
#   inst/extdata/io_yft_age_sf_fit.rds  AGE-SF, the article's model
#   inst/extdata/io_yft_bridge.rds      one summary row per fit in the
#                                       article's bridge table
#
# LEN and AGE start from SS3's solution. AGE-SF and AGE-SF-NP start from
# AGE's estimates: without priors, AGE-SF-NP started from SS3's solution
# stops at a different, worse optimum with a larger stock.
pkgload::load_all()
source("data-raw/io-yft-setup.R")

inputs <- io_yft_inputs()
control <- list(eval.max = 10000, iter.max = 10000)
metadata <- list(
  stock = "Indian Ocean yellowfin (simulated, replicate 4)",
  source = "SS3 model by M.N. Maunder, H. Xu and C. Minte-Vera (IATTC), converted to opal",
  converter = "https://github.com/N-DucharmeBarth/wcpfc-bet-single-region (code/rtmb)")
description <- c(
  LEN = "LEN: SS3 age-specific K growth, opal length-based selectivity",
  AGE = "AGE: SS3's age-based selectivity curves, every parameter SS3 estimated",
  `AGE-SF` = "AGE-SF: AGE plus SS3's index size compositions (lambda 0.51)",
  `AGE-SF-NP` = "AGE-SF-NP: AGE-SF with SS3's lambda of 0 on parameter priors")
build <- function(nm, parameters = inputs[[nm]]$parameters, map = inputs[[nm]]$map) {
  bounds <- if (nm == "LEN") NULL else io_yft_age_bounds(parameters, map)
  opal_obj(inputs[[nm]]$data, parameters, map, bounds = bounds, control = control,
           metadata = c(metadata, list(model = description[[nm]])))
}
par_list <- function(x, par) opal_rtmb(x, fresh = TRUE)$env$parList(par)
report_fit <- function(nm, x) {
  v <- x$validation$fit
  cat(sprintf("%s: objective %.4f, %d parameters, max gradient %.2e, check %s\n", nm,
              x$fit$opt$objective, length(x$fit$opt$par), v$metrics$max_gradient,
              isTRUE(v$passes)))
}

fits <- list()
for (nm in c("LEN", "AGE")) {
  fits[[nm]] <- io_yft_fit_polished(build(nm))
  report_fit(nm, fits[[nm]])
  stopifnot(isTRUE(fits[[nm]]$validation$fit$passes))
}
p_age <- par_list(fits$AGE, fits$AGE$fit$opt$par)

fits$`AGE-SF` <- io_yft_fit_polished(build("AGE-SF", p_age))
report_fit("AGE-SF", fits$`AGE-SF`)
stopifnot(isTRUE(fits$`AGE-SF`$validation$fit$passes))

# AGE-SF-NP ------------------------------------------------------------------
# Without the prior, several selectivity parameters are not identified: the
# likelihood is flat in them. The model is fitted freely first (AGE's map);
# parameters with a Hessian diagonal below 1e-3 at that estimate are then
# held at it, which leaves the objective unchanged, and any parameter still
# dominating a near-zero eigenvalue is added, before the final fit is
# polished and checked. The check is reported, not required: without priors
# the Hessian need not be positive definite.
hold <- function(map, cells) {
  lv <- matrix(as.integer(map$par_sel), nrow = inputs$LEN$data$n_fishery)
  lv[lv %in% lv[cells]] <- NA_integer_
  map$par_sel <- factor(as.vector(lv))
  map
}
first_cells <- function(map) {
  lv <- matrix(as.integer(map$par_sel), nrow = inputs$LEN$data$n_fishery)
  vapply(seq_len(max(lv, na.rm = TRUE)), function(k) which(lv == k)[1], integer(1))
}
x <- opal_fit(build("AGE-SF-NP", p_age), n_passes = 2L, check = FALSE)
hess <- opal_rtmb(x, fresh = TRUE)$he(x$fit$opt$par)
sel <- grep("^par_sel", names(x$fit$opt$par))
map <- hold(inputs$`AGE-SF-NP`$map, first_cells(inputs$`AGE-SF-NP`$map)[abs(diag(hess)[sel]) < 1e-3])
p_free <- par_list(x, x$fit$opt$par)
for (round in 1:3) {
  y <- opal_fit(build("AGE-SF-NP", p_free, map), n_passes = 2L, check = FALSE)
  e <- eigen(opal_rtmb(y, fresh = TRUE)$he(y$fit$opt$par), symmetric = TRUE)
  small <- which(abs(e$values) < 1e-6)
  if (!length(small)) break
  sel_y <- grep("^par_sel", names(y$fit$opt$par))
  add <- vapply(small, function(k) sel_y[which.max(abs(e$vectors[sel_y, k]))], integer(1))
  p_free <- par_list(y, y$fit$opt$par)
  map <- hold(map, first_cells(map)[match(add, sel_y)])
}
fits$`AGE-SF-NP` <- io_yft_fit_polished(build("AGE-SF-NP", p_free, map), guarded = TRUE)
report_fit("AGE-SF-NP", fits$`AGE-SF-NP`)

# Save -----------------------------------------------------------------------
change <- c(
  LEN = "Estimate opal's length-based selectivity",
  AGE = "SS3's age-based curves and mirrors; every parameter SS3 estimated",
  `AGE-SF` = "Add SS3's index size compositions (lambda 0.51)",
  `AGE-SF-NP` = "Switch off the priors (SS3's lambda of 0)")
bridge <- do.call(rbind, lapply(names(change), function(nm)
  io_yft_bridge_row(fits[[nm]], nm, change[[nm]])))
print(bridge)
saveRDS(bridge, file.path("inst", "extdata", "io_yft_bridge.rds"), compress = "xz", version = 3L)
opal_save(fits$`AGE-SF`, file.path("inst", "extdata", "io_yft_age_sf_fit.rds"),
          compress = "xz", overwrite = TRUE)
