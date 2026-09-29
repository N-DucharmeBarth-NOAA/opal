# Report a configured or fitted Opal model

Rebuilds an isolated RTMB objective when needed and evaluates the report
at the stored fitted point, or the configured point if the object is
unfitted. This returns point predictions, not posterior summaries.

## Usage

``` r
opal_report(x)
```

## Arguments

- x:

  An `opal_obj`.

## Value

The named model report at the stored parameters.

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
[`opal_update()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_update.md)

## Examples

``` r
assessment <- opal_read(system.file("extdata", "opaka_quickstart_fit.rds",
                                   package = "opal"))
#> Warning: The saved fit runtime identity differs because it was serialized under a different R version; verifying by rebuild.
#> Warning: The saved fit runtime identity differs because it was serialized under a different R version; verifying by rebuild.
head(opal_report(assessment)$spawning_biomass_y)
#> [1] 4029.267 3939.677 3803.894 3690.690 3623.467 3581.016
```
