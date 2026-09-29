# Create a portable Opal assessment object

Creates an S3 object before fitting. Configuration, fitted results, and
posterior samples remain ordinary R data. RTMB objectives are cached
only within the current session. Explicit inputs are recommended for
custom stocks; `model` identifies a bundled parameter configuration when
needed.

## Usage

``` r
opal_obj(
  data,
  parameters = NULL,
  map = NULL,
  random = character(),
  bounds = NULL,
  control = NULL,
  model = NULL,
  makeadfun_args = list(),
  metadata = list()
)
```

## Arguments

- data:

  Named model data list.

- parameters:

  Named initial parameter list, or NULL for bundled defaults.

- map:

  RTMB map. NULL resolves defaults;
  [`list()`](https://rdrr.io/r/base/list.html) explicitly frees all
  parameters.

- random:

  Names of random-effect parameters.

- bounds:

  A data frame with `lower` and `upper` columns, or a list of numeric
  `lower` and `upper` vectors in active RTMB parameter order. Use `NULL`
  to derive bounds with
  [`get_bounds()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/get_bounds.md).

- control:

  Optimiser controls.

- model:

  Optional bundled configuration name for
  [`get_parameters()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/get_parameters.md).

- makeadfun_args:

  Additional portable arguments for `RTMB::MakeADFun()`.

- metadata:

  Named user metadata.

## Value

An `opal_obj`. Construction never optimises or samples.

## Details

Use one object throughout the workflow: construct, optionally build,
fit, sample, report, and save.
[`summary()`](https://rdrr.io/r/base/summary.html) reports the lifecycle
stage separately from fit and MCMC checks. Results live in `$fit`,
`$mcmc`, `$validation`, and `$derived`; labels live in
`$provenance$metadata`. Use
[`opal_update()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_update.md)
to change configuration rather than assigning directly to scientific
fields.

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
[`opal_posterior()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_posterior.md),
[`opal_project()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_project.md),
[`opal_report()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_report.md),
[`opal_update()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_update.md)

## Examples

``` r
inputs <- opaka_quickstart_inputs()
x <- opal_obj(inputs$data, inputs$parameters, inputs$map)
summary(x)
#> <opal_obj> configured
#> Active parameters: NA 
#> Fit check: not run  | MCMC check: not run 
```
