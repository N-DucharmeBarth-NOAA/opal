# Sample accepted members of an assessment grid

Sample accepted members of an assessment grid

## Usage

``` r
opal_grid_mcmc(grid, ...)
```

## Arguments

- grid:

  An `opal_grid` returned by
  [`opal_grid()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_grid.md).

- ...:

  Arguments passed to
  [`opal_mcmc()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_mcmc.md).
  Use an explicit seed.

## Value

The grid with updated models and a `mcmc_passes` summary column. Failed
fits are skipped, and failed sampling attempts remain in each model.

## See also

Other assessment sensitivity:
[`opal_grid()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_grid.md),
[`opal_grid_draws()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_grid_draws.md),
[`opal_profile()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_profile.md),
[`plot_opal_profile()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/plot_opal_profile.md)
