# Select balanced posterior draws from a model grid

Select balanced posterior draws from a model grid

## Usage

``` r
opal_grid_draws(grid, n_per_model, seed = 123L)
```

## Arguments

- grid:

  An `opal_grid` with accepted posteriors in every member.

- n_per_model:

  Number of draws selected without replacement per model.

- seed:

  Reproducible selection seed, restored on exit.

## Value

A table of model, chain, iteration, retained draw index, and equal model
weights. Selection is balanced across chains within each model, with any
remainder allocated to the first chains. Models are not silently
dropped, and equal model weighting is an explicit assumption.

## See also

Other assessment sensitivity:
[`opal_grid()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_grid.md),
[`opal_grid_mcmc()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_grid_mcmc.md),
[`opal_profile()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_profile.md),
[`plot_opal_profile()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/plot_opal_profile.md)
