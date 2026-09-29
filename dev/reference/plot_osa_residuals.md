# Plot stored OSA residuals

Plot stored OSA residuals

## Usage

``` r
plot_osa_residuals(x, name = "osa", type = c("time", "qq"))
```

## Arguments

- x:

  An Opal object containing stored OSA diagnostics.

- name:

  Name of the stored diagnostic result.

- type:

  Plot residuals against model year or as normal quantiles.

## Value

A ggplot, faceted by dataset.

## See also

Other assessment diagnostics:
[`opal_diagnose()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_diagnose.md),
[`opal_osa()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_osa.md),
[`plot_composition()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/plot_composition.md),
[`plot_osa_sdnr()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/plot_osa_sdnr.md),
[`plot_prior_posterior()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/plot_prior_posterior.md)
