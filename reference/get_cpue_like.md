# CPUE index likelihood

Computes the likelihood for a standardized CPUE index using a log-linear
model.

## Usage

``` r
get_cpue_like(data, parameters, number_ysa, sel_fya, creep_init = 1)
```

## Arguments

- data:

  a `list` of data inputs (cpue_data, cpue_switch, etc.).

- parameters:

  a `list` of parameter values (log_cpue_tau, log_cpue_omega,
  cpue_creep, log_cpue_q, weight_fya, etc.).

- number_ysa:

  a 3D `array` `[n_year, n_season, n_age]` of numbers-at-age.

- sel_fya:

  a 3D `array` `[n_fishery, n_year, n_age]` of selectivity by fishery,
  year, and age.

- creep_init:

  scalar initialization value for creeping adjustment (default 1).

## Value

a `numeric` vector of negative log-likelihood contributions.
