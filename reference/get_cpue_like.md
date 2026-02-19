# CPUE index likelihood

Computes the likelihood for a standardized CPUE index using a log-linear
model.

## Usage

``` r
get_cpue_like(data, parameters, number_ysa, sel_fya, creep_init = 1)
```

## Arguments

- data:

  Integer switch to activate the likelihood.

- parameters:

  a `vector` of year indices for CPUE observations.

- number_ysa:

  a 3D `array` year, season, age of numbers-at-age.

- sel_fya:

  a 3D `array` of selectivity by fishery, year, and age.

## Value

a `list` with predicted CPUE, residuals, and likelihood vector.
