# Recreate the simulated posterior used by the assessment-tools article.
# Run from the package root. Sampling is performed here, not during website builds.
devtools::load_all()
inputs <- opal_example_inputs()
x <- opal_obj(inputs$data,inputs$parameters,inputs$map,
  metadata=list(stock="Simulated tutorial",purpose="Software example, not stock advice"))
# Rebuild reports and diagnostics from verified saved estimates/draws when
# only the reporting contract changes, avoiding cross-toolchain fit drift.
previous <- Sys.getenv("OPAL_REBUILD_FROM")
if (nzchar(previous)) {
  old <- readRDS(previous)$object
  stopifnot(isTRUE(all.equal(x$data, old$data)),
            isTRUE(all.equal(x$parameters, old$parameters)),
            identical(x$map, old$map))
  object <- opal_rtmb(x, fresh = TRUE)
  stopifnot(identical(names(object$par), names(old$fit$opt$par)),
            abs(object$fn(old$fit$opt$par) - old$fit$opt$objective) < 1e-8,
            all(is.finite(object$gr(old$fit$opt$par))))
  x <- opal_attach_fit(x, old$fit$opt)
} else {
  x <- opal_fit(x)
}
stopifnot(isTRUE(x$validation$fit$passes))
if (nzchar(previous)) {
  x <- opal_attach_mcmc(x, old$mcmc, settings = old$mcmc$settings)
  stopifnot(identical(x$mcmc$samples, old$mcmc$samples),
            identical(x$mcmc$sampler_params, old$mcmc$sampler_params))
} else {
  x <- opal_mcmc(x,seed=927,chains=4,cores=1,num_warmup=1000,num_samples=1000,
                 adapt_delta=.95)
}
stopifnot(isTRUE(x$validation$mcmc$passes))
# Store both the draw identities and verified analysis outputs for reproducibility.
x <- opal_posterior(x,quantities=c("B0","spawning_biomass_y","static_depletion_y"))
x <- opal_osa(x)
opal_save(x,'inst/extdata/assessment_example.rds',overwrite=TRUE)
print(summary(x))
print(x$validation$mcmc$metrics$parameters)
