# Attach an externally optimised Opal fit

Attach a result from an external optimiser to the configuration that
produced it. The objective, parameter layout, and bounds are verified
before retaining the point. For ordinary optimisation, use
[`opal_fit()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_fit.md).

## Usage

``` r
opal_attach_fit(x, opt, check = TRUE, check_args = list())
```

## Arguments

- x:

  A configured `opal_obj`.

- opt:

  Optimiser result with named `par`, `objective`, and `convergence`.

- check:

  Run fitting diagnostics after attachment.

- check_args:

  Named arguments passed to
  [`opal_check()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_check.md).

## Value

A fitted `opal_obj`.

## See also

Other assessment workflow:
[`opal_attach_mcmc()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_attach_mcmc.md),
[`opal_build()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_build.md),
[`opal_check()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_check.md),
[`opal_fit()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_fit.md),
[`opal_from_fit()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_from_fit.md),
[`opal_io`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_io.md),
[`opal_mcmc()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_mcmc.md),
[`opal_obj()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_obj.md),
[`opal_project()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_project.md),
[`opal_report()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_report.md),
[`opal_update()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_update.md)

## Examples

``` r
assessment <- opal_read(system.file("extdata", "opaka_quickstart_fit.rds",
                                   package = "opal"))
#> Warning: The saved fit runtime identity differs because it was serialized under a different R version; verifying by rebuild.
#> Warning: The saved fit runtime identity differs because it was serialized under a different R version; verifying by rebuild.
# Reattach the already verified optimum; no optimisation is performed.
assessment <- opal_attach_fit(assessment, assessment$fit$opt, check = FALSE)
```
