# Attach posterior draws to an Opal object

Validate and store posterior draws for the object's current
configuration. Parameter names must identify the active model layout.
Matrices contain draws in rows and named variables in columns; arrays
use iteration, chain, and variable dimensions. Use
[`opal_mcmc()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_mcmc.md)
to run a new sampler.

## Usage

``` r
opal_attach_mcmc(x, mcmc, settings = list(), check = TRUE, check_args = list())
```

## Arguments

- x:

  A configured or fitted `opal_obj`.

- mcmc:

  SparseNUTS output, sample matrix, or iteration-chain-variable array.

- settings:

  Named portable sampler settings.

- check:

  Run posterior diagnostics.

- check_args:

  Arguments passed to
  [`opal_check()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_check.md).

## Value

An updated `opal_obj`. Previous selected draws are retained in history.

## Details

A failed or unchecked import does not replace a previously checked,
passing posterior; inspect `x$mcmc_history` for unselected attempts.
Parameter-only imports lack sampler diagnostics, so a passing MCMC check
requires richer sampler output. Full joint draws are required for
random-effect projections.

## See also

Other assessment workflow:
[`opal_attach_fit()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_attach_fit.md),
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
[`opal_report()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_report.md),
[`opal_update()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_update.md)
