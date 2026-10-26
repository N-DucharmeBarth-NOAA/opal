# Recreate the stored Indian Ocean yellowfin fit used by the
# io-yft-ss3-comparison article. Run from the package root; takes a few
# minutes. The inputs (data/io_yft_*.rda) were exported from an SS3 model by
# code/rtmb/scripts/04-export-io-yft-for-opal.r in
# https://github.com/N-DucharmeBarth/wcpfc-bet-single-region.
pkgload::load_all()

x <- opal_obj(io_yft_data, io_yft_parameters, io_yft_map,
              control = list(eval.max = 10000, iter.max = 10000))
x <- opal_fit(x, n_passes = 2L, check = FALSE)

# nlminb stops on relative function change, which for an objective near 1.4e4
# leaves gradients of order 1e-3. A few Newton steps take them below 1e-8.
object <- opal_rtmb(x, fresh = TRUE)
opt <- x$fit$opt
for (i in seq_len(5L)) {
  g <- as.numeric(object$gr(opt$par))
  if (max(abs(g)) < 1e-8) break
  opt$par <- stats::setNames(opt$par - as.numeric(solve(object$he(opt$par), g)),
                             names(opt$par))
}
opt$objective <- object$fn(opt$par)

x <- opal_attach_fit(x, opt, check = TRUE,
                     check_args = list(gradient_tolerance = 1e-4))
stopifnot(isTRUE(x$validation$fit$passes))
x <- opal_update(x, metadata = list(
  stock = "Indian Ocean yellowfin (simulated, replicate 4)",
  model = "C2: SS3 age-specific K growth, opal length-based selectivity",
  source = "SS3 model by M.N. Maunder, H. Xu and C. Minte-Vera (IATTC), converted to opal",
  converter = "https://github.com/N-DucharmeBarth/wcpfc-bet-single-region (code/rtmb)"))

opal_save(x, "inst/extdata/io_yft_c2_fit.rds", compress = "xz", overwrite = TRUE)
print(summary(x))
