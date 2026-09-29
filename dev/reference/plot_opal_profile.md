# Plot an objective profile

Plot an objective profile

## Usage

``` r
plot_opal_profile(x, name = "profile")
```

## Arguments

- x:

  An Opal object with a stored profile.

- name:

  Name of the stored profile.

## Value

A ggplot. Failed conditional fits are shown as crosses when their
objectives are available; errors remain in the stored table.

## See also

Other assessment sensitivity:
[`opal_grid()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_grid.md),
[`opal_grid_draws()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_grid_draws.md),
[`opal_grid_mcmc()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_grid_mcmc.md),
[`opal_profile()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_profile.md)
