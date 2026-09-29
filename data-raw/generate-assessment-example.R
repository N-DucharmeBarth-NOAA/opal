# Recreate the simulated posterior used by the assessment-tools article.
# Run from the package root. Sampling is performed here, not during website builds.
devtools::load_all()
inputs <- opal_example_inputs()
x <- opal_fit(opal_obj(inputs$data,inputs$parameters,inputs$map,
  metadata=list(stock="Simulated tutorial",purpose="Software example, not stock advice")))
stopifnot(isTRUE(x$validation$fit$passes))
x <- opal_mcmc(x,seed=927,chains=4,cores=1,num_warmup=1000,num_samples=1000,
               adapt_delta=.95)
stopifnot(isTRUE(x$validation$mcmc$passes))
# Store both the draw identities and verified analysis outputs for reproducibility.
x <- opal_posterior(x,quantities=c("B0","spawning_biomass_y","static_depletion_y"))
x <- opal_osa(x)
opal_save(x,'inst/extdata/assessment_example.rds',overwrite=TRUE)
print(summary(x))
print(x$validation$mcmc$metrics$parameters)
