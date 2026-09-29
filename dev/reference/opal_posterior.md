# Summarise posterior parameters and derived model quantities

Evaluates each selected joint draw using its own biological inputs and
reports. Summaries are stored in the portable assessment object.

## Usage

``` r
opal_posterior(
  x,
  quantities = c("spawning_biomass_y", "static_depletion_y", "dynamic_depletion_y"),
  probs = c(0.025, 0.5, 0.975),
  draws = NULL,
  name = "posterior"
)
```

## Arguments

- x:

  An `opal_obj` with posterior draws.

- quantities:

  Numeric report names to summarise.

- probs:

  Distinct, increasing quantile probabilities between zero and one.

- draws:

  Optional retained draw indices, in iteration-within-chain order.

- name:

  Name of the stored result.

## Value

An updated object.
[`opal_derived()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_derived.md)
returns long-form `parameters` and `reports` summaries, report
dimensions, draw identifiers, and the source MCMC validation status.
Element indices follow R's column order.

## Details

Summaries never establish posterior acceptance. Non-finite report values
cause an error, rather than silently discarding draws. Marginal
posterior samples lacking random effects cannot be used. Use
[`opal_check()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_check.md)
to assess mixing, sampler diagnostics, bounds, and biology.

## See also

Other assessment workflow:
[`opal_attach_fit()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_attach_fit.md),
[`opal_attach_mcmc()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_attach_mcmc.md),
[`opal_build()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_build.md),
[`opal_check()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_check.md),
[`opal_derived()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_derived.md),
[`opal_example_inputs()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_example_inputs.md),
[`opal_fit()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_fit.md),
[`opal_from_fit()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_from_fit.md),
[`opal_io`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_io.md),
[`opal_mcmc()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_mcmc.md),
[`opal_obj()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_obj.md),
[`opal_project()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_project.md),
[`opal_report()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_report.md),
[`opal_update()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_update.md)
