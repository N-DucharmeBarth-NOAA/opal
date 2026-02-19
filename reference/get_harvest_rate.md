# Harvest rate calculation

Computes age-specific harvest rates by fishery for a single year-season
combination, using the Baranov catch equation.

## Usage

``` r
get_harvest_rate(data, y, s, number_ysa, sel_fya, weight_fya)
```

## Arguments

- data:

  A `list` of model data. Must contain: `n_fishery`, `n_age`,
  `catch_obs_ysf`, `catch_units_f`.

- y:

  Integer. Year index (1-based).

- s:

  Integer. Season index (1-based).

- number_ysa:

  Numeric array `[n_year+1, n_season, n_age]`. Current numbers-at-age.

- sel_fya:

  Numeric array `[n_fishery, n_year, n_age]`. Selectivity at age by
  fishery and year.

- weight_fya:

  Numeric array `[n_fishery, n_year, n_age]`. Mean weight at age by
  fishery and year. Passed explicitly so AD gradients propagate if
  growth is ever estimated.

## Value

A named list with:

- h_rate_fa:

  Harvest rate array `[n_fishery, n_age]`.

- h_rate_a:

  Total harvest rate vector (length `n_age`).

- penalty:

  Penalty from
  [`posfun`](https://n-ducharmebarth-noaa.github.io/opal/reference/posfun.md)
  for harvest-rate constraint.

## Details

`weight_fya` is passed as an explicit argument (not read from `data`) so
that AD gradients propagate correctly if growth parameters are estimated
in the future.
