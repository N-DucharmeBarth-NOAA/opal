# Compare parameter priors and posterior distributions

Compare parameter priors and posterior distributions

## Usage

``` r
plot_prior_posterior(x, parameters = NULL)
```

## Arguments

- x:

  An `opal_obj` with joint posterior samples and `data$priors`.

- parameters:

  Optional expanded active parameter names to plot.

## Value

A ggplot with marginal posterior histograms and prior densities on each
parameter's stored scale. Fixed and shared mapped elements are omitted
when there is no unambiguous scalar prior correspondence.

## Details

The displayed priors are the specified, untruncated densities.
Optimiser/sampler bounds may truncate the realised target. For shared
parameter maps, multiple prior contributions need not equal one density;
such parameters and blocks with multiple priors are excluded with a
warning.

## See also

Other assessment diagnostics:
[`opal_diagnose()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_diagnose.md),
[`opal_osa()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_osa.md),
[`plot_composition()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/plot_composition.md),
[`plot_osa_residuals()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/plot_osa_residuals.md),
[`plot_osa_sdnr()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/plot_osa_sdnr.md)
