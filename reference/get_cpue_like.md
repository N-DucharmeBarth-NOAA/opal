# CPUE index likelihood (multi-index)

Computes the likelihood for one or more standardised CPUE/survey indices
using a log-linear model. Each index has its own catchability (q), extra
variance (tau), power parameter (omega), and effort creep.
Mean-centering of predicted CPUE is performed within each index.

## Usage

``` r
get_cpue_like(
  cpue_data,
  parameters,
  number_ysa,
  sel_fya,
  weight_fya,
  cpue_switch = 1L
)
```

## Arguments

- cpue_data:

  a `list` of data inputs. Must contain:

  cpue_data

  :   data.frame with columns `ts`, `fishery`, `value`, `se`, `units`,
      and `index`.

  cpue_switch

  :   integer switch (0 = skip likelihood).

  n_index

  :   integer number of distinct indices.

- parameters:

  a `list` of parameter values. Must contain:

  log_cpue_q

  :   numeric vector `[n_index]`.

  log_cpue_tau

  :   numeric vector `[n_index]`.

  log_cpue_omega

  :   numeric vector `[n_index]`.

- number_ysa:

  a 3D `array` `[n_year, n_season, n_age]` of numbers-at-age.

- sel_fya:

  a 3D `array` `[n_fishery, n_year, n_age]` of selectivity by fishery,
  year, and age.

- weight_fya:

  a 3D `array` `[n_fishery, n_year, n_age]` of weight-at-age by fishery
  and year.

- cpue_switch:

  boolean flag to calculate the cpue likelihood.

## Value

numeric vector of length `nrow(cpue_data)` with per-observation negative
log-likelihood contributions.
