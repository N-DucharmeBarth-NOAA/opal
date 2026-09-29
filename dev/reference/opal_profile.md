# Profile an assessment parameter

Fixes a scalar parameter or a mapped parameter group at each supplied
value, then re-optimises the remaining parameters using the usual
two-pass fit.

## Usage

``` r
opal_profile(
  x,
  parameter,
  values,
  element = 1L,
  fit_args = list(),
  name = "profile"
)
```

## Arguments

- x:

  A fitted `opal_obj`.

- parameter:

  Name of a parameter block, on its stored scale.

- values:

  Distinct finite profile values on that scale.

- element:

  One-based element within the block, in R column order.

- fit_args:

  Named arguments passed to
  [`opal_fit()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_fit.md).

- name:

  Name of the stored result.

## Value

An updated object containing a table, component contributions,
conditional fits, and failed-point messages, accessible with
[`opal_derived()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_derived.md).

## Details

This profiles the complete penalised objective, including priors and
process penalties; it is a likelihood profile only when those terms are
absent. Fixed effects can be profiled in a Laplace model, but random
effect blocks cannot. Shared map elements are fixed together. The
original fit is the reference; a negative objective difference indicates
that a conditional fit improved on it. Failed points are never
interpolated.

## See also

Other assessment sensitivity:
[`opal_grid()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_grid.md),
[`opal_grid_draws()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_grid_draws.md),
[`opal_grid_mcmc()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_grid_mcmc.md),
[`plot_opal_profile()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/plot_opal_profile.md)
