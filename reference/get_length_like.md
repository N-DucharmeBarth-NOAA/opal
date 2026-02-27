# Length Composition Likelihood

Computes likelihood for observed length compositions using a
probability-of-length-at-age (PLA) matrix. Gradients propagate through
to growth parameters when they are estimated.

## Usage

``` r
get_length_like(
  lf_obs_flat,
  lf_obs_ints,
  lf_obs_prop,
  catch_pred_fya,
  pla,
  lf_n_f,
  lf_fishery_f,
  lf_year_fi,
  lf_n_fi,
  lf_minbin,
  lf_maxbin,
  removal_switch_f,
  lf_switch,
  n_len,
  n_lf,
  log_lf_tau
)
```

## Arguments

- lf_obs_flat:

  numeric vector of unrounded length composition counts (used when
  lf_switch=1).

- lf_obs_ints:

  integer vector of integer length composition counts (used when
  lf_switch=3).

- lf_obs_prop:

  numeric vector of length composition proportions (used when
  lf_switch=2).

- catch_pred_fya:

  3D `array` `[n_fishery, n_year, n_age]` of predicted catch-at-age from
  [`do_dynamics()`](https://n-ducharmebarth-noaa.github.io/opal/reference/do_dynamics.md).

- pla:

  matrix `[n_len, n_age]` probability-of-length-at-age from
  [`get_pla()`](https://n-ducharmebarth-noaa.github.io/opal/reference/get_pla.md).
  On the AD tape when growth parameters are estimated.

- lf_n_f:

  integer vector `[n_fishery]` of number of observations per fishery.

- lf_fishery_f:

  integer vector of fishery indices for each observation group.

- lf_year_fi:

  list of integer vectors of year indices for each observation.

- lf_n_fi:

  list of integer vectors of sample sizes for each observation.

- lf_minbin:

  integer vector `[n_fishery]` of minimum bin index per fishery.

- lf_maxbin:

  integer vector `[n_fishery]` of maximum bin index per fishery.

- removal_switch_f:

  integer vector `[n_fishery]` indicating if fishery is removed
  (0=included, 1=removed).

- lf_switch:

  integer likelihood type selector (1=multinomial, 2=dirichlet,
  3=dirichlet-multinomial).

- n_len:

  integer number of length bins.

- n_lf:

  integer total number of length composition observations.

- log_lf_tau:

  numeric vector `[n_fishery]` of log-scale variance adjustment
  parameters.

## Value

a `numeric` vector of negative log-likelihood contributions, one per
observation.
