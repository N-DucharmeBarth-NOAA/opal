# Convert a legacy Opal fit to the staged object workflow

Verifies the saved fitted objective without optimising or sampling.
Existing diagnostics, posterior draws, derived outputs, and provenance
are preserved.

## Usage

``` r
opal_from_fit(x, integrity = c("portable", "exact"))
```

## Arguments

- x:

  A legacy `opal_fit` or an `opal_obj`.

- integrity:

  Legacy runtime verification mode; portable permits verified cross-R
  reads.

## Value

An `opal_obj`.

## See also

Other assessment workflow:
[`opal_attach_fit()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_attach_fit.md),
[`opal_attach_mcmc()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_attach_mcmc.md),
[`opal_build()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_build.md),
[`opal_check()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_check.md),
[`opal_derived()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_derived.md),
[`opal_example_inputs()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_example_inputs.md),
[`opal_fit()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_fit.md),
[`opal_io`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_io.md),
[`opal_mcmc()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_mcmc.md),
[`opal_obj()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_obj.md),
[`opal_posterior()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_posterior.md),
[`opal_project()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_project.md),
[`opal_report()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_report.md),
[`opal_update()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_update.md)

## Examples

``` r
legacy <- read_opal_fit(system.file("extdata", "opaka_quickstart_fit.rds",
                                   package = "opal"), integrity = "portable")
#> Warning: The saved fit runtime identity differs because it was serialized under a different R version; verifying by rebuild.
#> Warning: Portable runtime identity was not independently verified because `rebuild = FALSE`.
assessment <- opal_from_fit(legacy)
#> Warning: The saved fit runtime identity differs because it was serialized under a different R version; verifying by rebuild.
#> Warning: The saved fit runtime identity differs because it was serialized under a different R version; verifying by rebuild.
summary(assessment)
#> <opal_obj> fitted
#> Active parameters: 127 
#> Objective: 1832.412 
#> Fit check: not run  | MCMC check: not run 
```
