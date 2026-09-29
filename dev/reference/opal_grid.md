# Fit a reproducible grid of assessment scenarios

Fit a reproducible grid of assessment scenarios

## Usage

``` r
opal_grid(x, scenarios, directory = NULL, resume = TRUE, fit_args = list())
```

## Arguments

- x:

  A configured or fitted `opal_obj`.

- scenarios:

  Named list of argument lists passed to
  [`opal_update()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_update.md).

- directory:

  Optional directory for per-scenario portable checkpoints.

- resume:

  Reuse checkpoints only when target and fitting settings match.

- fit_args:

  Named arguments passed to
  [`opal_fit()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_fit.md).

## Value

An `opal_grid` list with `models`, `summary`, and `settings`. Every
scenario is retained, including failures. `passes` uses the full fit
check, not just the optimiser code. Checkpoint identities include the
scientific configuration and fitting controls.

## See also

Other assessment sensitivity:
[`opal_grid_draws()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_grid_draws.md),
[`opal_grid_mcmc()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_grid_mcmc.md),
[`opal_profile()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_profile.md),
[`plot_opal_profile()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/plot_opal_profile.md)
