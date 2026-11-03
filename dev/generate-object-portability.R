# Regenerate only the object-format fixture. Scientific goldens are separate.
# Run deliberately from the repository root under the R version being tested.
stopifnot(identical(Sys.getenv("OPAL_REFRESH_OBJECT_FIXTURE"), "true"))
pkgload::load_all()
source("tests/testthat/helper-opal-fit.R")
source("tests/testthat/helper-opal-object.R")
previous <- Sys.getenv("OPAL_REBUILD_FROM")
if (nzchar(previous)) {
  old <- readRDS(previous)$object
  x <- opal_obj(old$data, old$parameters, old$map, random = old$random,
    bounds = old$bounds, control = old$control,
    makeadfun_args = old$makeadfun_args, metadata = old$provenance$metadata)
  object <- opal_rtmb(x, fresh = TRUE)
  stopifnot(identical(names(object$par), names(old$fit$opt$par)),
            abs(object$fn(old$fit$opt$par) - old$fit$opt$objective) < 1e-8,
            all(is.finite(object$gr(old$fit$opt$par))))
  x <- opal_attach_fit(x, old$fit$opt, check = FALSE)
} else {
  x <- opal_fit(small_opal_object(), check = FALSE)
}
draws <- matrix(c(x$fit$opt$par, x$fit$opt$par + 0.001), ncol = 1,
                 dimnames = list(NULL, "log_B0"))
x <- opal_attach_mcmc(x, draws, check = FALSE)
dir.create("tests/testthat/_fixtures", showWarnings = FALSE)
opal_save(x, "tests/testthat/_fixtures/opal-obj-v1.rds", overwrite = TRUE)
