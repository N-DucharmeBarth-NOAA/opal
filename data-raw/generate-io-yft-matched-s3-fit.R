# Fit S3 for the io-yft-ss3-matched article: S2 with opal's default priors
# kept, i.e. A-SS3 plus SS3's index size compositions (lambda 0.51). With
# the priors, every selectivity parameter SS3 estimated is identified, so it
# uses A-SS3's map and starts from A-SS3's estimates, like S1 and S2. Run
# from the package root; takes about 10 minutes. Afterwards, rerun
# generate-io-yft-vignette-cache.R.
#
#   inst/extdata/io_yft_matched_s3_fit.rds
pkgload::load_all()
source("data-raw/io-yft-ss3-matched-setup.R")

inputs <- io_yft_matched_inputs()
a_ss3 <- opal_read(file.path("inst", "extdata", "io_yft_age_ss3_fit.rds"))
p_a <- opal_rtmb(a_ss3, fresh = TRUE)$env$parList(a_ss3$fit$opt$par)

x <- opal_obj(inputs$data_s3, p_a, inputs$map_ss3,
              bounds = io_yft_age_bounds(p_a, inputs$map_ss3),
              control = list(eval.max = 10000, iter.max = 10000),
              metadata = list(
                stock = "Indian Ocean yellowfin (simulated, replicate 4)",
                model = "S3: as S2, with opal's default priors (A-SS3 plus the index size compositions)",
                source = "SS3 model by M.N. Maunder, H. Xu and C. Minte-Vera (IATTC), converted to opal",
                converter = "https://github.com/N-DucharmeBarth/wcpfc-bet-single-region (code/rtmb)"))
t0 <- proc.time()[["elapsed"]]
x <- io_yft_fit_polished(x)
cat(sprintf("S3: objective %.4f, %d parameters, check %s, %.0f s\n", x$fit$opt$objective,
            length(x$fit$opt$par), isTRUE(x$validation$fit$passes), proc.time()[["elapsed"]] - t0))
print(x$validation$fit$checks)
opal_save(x, file.path("inst", "extdata", "io_yft_matched_s3_fit.rds"),
          compress = "xz", overwrite = TRUE)
