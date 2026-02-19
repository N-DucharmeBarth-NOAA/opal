# Length Composition Likelihood

Computes likelihood for observed length compositions using ALKs.

## Usage

``` r
get_length_like(data, parameters, log_lf_alpha, catch_pred_fya, alk_ysal)
```

## Arguments

- catch_pred_fya:

  3D array of predicted catch.

- alk_ysal:

  4D array year, season, age, length_bin of ALKs.

- lf_switch:

  Vectors of indices for observations.

- removal_switch_f:

  Vectors of indices for observations.

- lf_year, lf_season, lf_fishery:

  Vectors of indices for observations.

- lf_minbin:

  Minimum size bin to be aggregated.

- lf_obs:

  Matrix of observed length proportions.

- lf_n:

  Vector of effective sample sizes.

- par_log_lf_alpha:

  Vector of effective sample sizes.

## Value

List with predicted compositions and likelihood contributions.
