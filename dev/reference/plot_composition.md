# Plot observed and fitted compositions

Plot observed and fitted compositions

## Usage

``` r
plot_composition(x, type = c("lf", "wf"), fishery = NULL)
```

## Arguments

- x:

  A fitted `opal_obj`.

- type:

  Length (`"lf"`) or weight (`"wf"`) compositions.

- fishery:

  Optional fishery indices to include.

## Value

A ggplot with observation proportions and fitted probabilities. Removed
and zero-size rows are excluded. Bin numbers refer to the prepared data,
including any tail aggregation.

## See also

Other assessment diagnostics:
[`opal_diagnose()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_diagnose.md),
[`opal_osa()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_osa.md),
[`plot_osa_residuals()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/plot_osa_residuals.md),
[`plot_osa_sdnr()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/plot_osa_sdnr.md),
[`plot_prior_posterior()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/plot_prior_posterior.md)
