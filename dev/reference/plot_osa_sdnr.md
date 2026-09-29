# Compare OSA residual dispersion across datasets

Compare OSA residual dispersion across datasets

## Usage

``` r
plot_osa_sdnr(x, name = "osa")
```

## Arguments

- x:

  An Opal object containing stored OSA diagnostics.

- name:

  Name of the stored diagnostic result.

## Value

A ggplot with one dataset per row, SDNR on the horizontal axis,
approximate intervals, and a reference line at one. Plot data contain
all dataset summaries, including counts of excluded and failed
residuals.

## See also

Other assessment diagnostics:
[`opal_diagnose()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_diagnose.md),
[`opal_osa()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_osa.md),
[`plot_composition()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/plot_composition.md),
[`plot_osa_residuals()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/plot_osa_residuals.md),
[`plot_prior_posterior()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/plot_prior_posterior.md)
