# Fit the Indian Ocean yellowfin example with age-based selectivity
# (sel_basis_f = 2), starting from SS3's own age_selex parameters, for the
# io-yft-age-selectivity article. Run from the package root; takes several
# minutes. Afterwards, rerun generate-io-yft-vignette-cache.R.
#
#   inst/extdata/io_yft_age_fit.rds      A: C2's map (49 selectivity
#                                        parameters), age-based curves
#   inst/extdata/io_yft_age_ss3_fit.rds  A-SS3: every selectivity parameter
#                                        SS3 estimated (58)
pkgload::load_all()
source("data-raw/io-yft-age-setup.R")

inputs <- io_yft_age_inputs()
control <- list(eval.max = 10000, iter.max = 10000)
metadata <- list(
  stock = "Indian Ocean yellowfin (simulated, replicate 4)",
  source = "SS3 model by M.N. Maunder, H. Xu and C. Minte-Vera (IATTC), converted to opal",
  converter = "https://github.com/N-DucharmeBarth/wcpfc-bet-single-region (code/rtmb)")

fits <- list(
  io_yft_age_fit = list(
    map = inputs$map_c2,
    model = "A: SS3 age-specific K growth, opal age-based selectivity, C2's map"),
  io_yft_age_ss3_fit = list(
    map = inputs$map_ss3,
    model = "A-SS3: as A, estimating every selectivity parameter SS3 estimated"))

for (nm in names(fits)) {
  x <- opal_obj(inputs$data, inputs$parameters, fits[[nm]]$map,
                bounds = io_yft_age_bounds(inputs$parameters, fits[[nm]]$map),
                control = control,
                metadata = c(metadata, list(model = fits[[nm]]$model)))
  t0 <- proc.time()[["elapsed"]]
  x <- io_yft_fit_polished(x)
  cat(sprintf("%s: objective %.4f, %d parameters, check %s, %.0f s\n", nm,
              x$fit$opt$objective, length(x$fit$opt$par),
              isTRUE(x$validation$fit$passes), proc.time()[["elapsed"]] - t0))
  stopifnot(isTRUE(x$validation$fit$passes))
  opal_save(x, file.path("inst", "extdata", paste0(nm, ".rds")),
            compress = "xz", overwrite = TRUE)
}
