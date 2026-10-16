# Regenerate only the object-format fixture. Scientific goldens are separate.
# Run deliberately from the repository root under the R version being tested.
stopifnot(identical(Sys.getenv("OPAL_REFRESH_OBJECT_FIXTURE"), "true"))
pkgload::load_all()
source("tests/testthat/helper-opal-fit.R")
source("tests/testthat/helper-opal-object.R")
x <- opal_fit(small_opal_object(), check = FALSE)
draws <- matrix(c(x$fit$opt$par, x$fit$opt$par + 0.001), ncol = 1,
                 dimnames = list(NULL, "log_B0"))
x <- opal_attach_mcmc(x, draws, check = FALSE)
dir.create("tests/testthat/_fixtures", showWarnings = FALSE)
opal_save(x, "tests/testthat/_fixtures/opal-obj-v1.rds", overwrite = TRUE)
