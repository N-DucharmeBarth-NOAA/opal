# Update an Opal model configuration

Changes to the scientific configuration clear dependent fits, samples,
and reports. Metadata and optimiser controls do not alter the scientific
target. Defaults are resolved again when their upstream configuration
changes.

## Usage

``` r
opal_update(
  x,
  data,
  parameters,
  map,
  random,
  bounds,
  control,
  makeadfun_args,
  priors,
  metadata
)
```

## Arguments

- x:

  An `opal_obj`.

- data, parameters, map, random, bounds, control, makeadfun_args:

  As in
  [`opal_obj()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_obj.md).

- priors:

  Optional replacement for `data$priors`.

- metadata:

  Named metadata to merge.

## Value

An updated `opal_obj`; the input object is not modified.

## Details

Changes to data, priors, parameters, maps, random effects, bounds, or
RTMB construction settings clear dependent fits, posteriors,
diagnostics, and projections. Metadata is merged, and optimiser controls
can change without invalidating results. A supplied map is retained
unless explicitly replaced; check its dimensions when changing parameter
structure.

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
[`opal_posterior()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_posterior.md),
[`opal_project()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_project.md),
[`opal_report()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_report.md)

## Examples

``` r
inputs <- opaka_quickstart_inputs()
assessment <- opal_obj(inputs$data, inputs$parameters, inputs$map)
labelled <- opal_update(assessment, metadata = list(label = "Baseline"))
summary(labelled)
#> <opal_obj> configured
#> Active parameters: NA 
#> Fit check: not run  | MCMC check: not run 
```
