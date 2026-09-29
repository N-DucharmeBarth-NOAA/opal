# Diagnose biological feasibility of an Opal model

Checks the initial-equilibrium and harvest penalties, population states,
biological inputs, catch reconstruction, and derived depletion. A
numerical optimum alone does not imply that these checks pass.

## Usage

``` r
opal_diagnose(x, penalty_tolerance = 1e-10, catch_tolerance = 1e-06)
```

## Arguments

- x:

  A configured or fitted `opal_obj`.

- penalty_tolerance:

  Maximum allowed positive continuation penalty.

- catch_tolerance:

  Maximum catch error divided by `pmax(1, observed)`.

## Value

A list with `passes`, named logical `checks`, and numerical `metrics`.

## Details

Zero abundance is allowed, but unfished spawning output must be
positive. Catch is conditioned on, so its difference is a reconstruction
check, not an observation residual. Checks use the reported,
draw-specific biology. They do not assess identifiability or scientific
model adequacy.

## See also

Other assessment diagnostics:
[`opal_osa()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_osa.md),
[`plot_composition()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/plot_composition.md),
[`plot_osa_residuals()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/plot_osa_residuals.md),
[`plot_osa_sdnr()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/plot_osa_sdnr.md),
[`plot_prior_posterior()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/plot_prior_posterior.md)
