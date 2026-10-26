# Plot observed and fitted compositions

Plot observed and fitted compositions

## Usage

``` r
plot_composition(
  x,
  type = c("lf", "wf"),
  fishery = NULL,
  scale = c("auto", "proportion", "density")
)
```

## Arguments

- x:

  A fitted `opal_obj`.

- type:

  Length (`"lf"`) or weight (`"wf"`) compositions.

- fishery:

  Optional fishery indices to include.

- scale:

  Display proportions or probability densities. `"auto"` uses density
  when bin widths differ or a fishery-year has different bin schemes.

## Value

A ggplot with observed points and a fitted line at actual bin midpoints.
Horizontal observed and fitted segments show the bin extents. Removed
and zero-size rows are excluded. Density divides both observed
proportions and fitted probabilities by bin width.

## See also

Other assessment diagnostics:
[`opal_diagnose()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_diagnose.md),
[`opal_osa()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_osa.md),
[`plot_osa_residuals()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/plot_osa_residuals.md),
[`plot_osa_sdnr()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/plot_osa_sdnr.md),
[`plot_prior_posterior()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/plot_prior_posterior.md)
