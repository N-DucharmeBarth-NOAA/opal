# Project selectivity

Project selectivity

## Usage

``` r
project_selectivity(
  data,
  obj,
  first_yr = NULL,
  last_yr = NULL,
  n_proj = 5,
  n_iter = 1
)
```

## Arguments

- data:

  a `list` of parameter values.

- obj:

  a `list` of parameter values.

- first_yr:

  the first year sampled. Defaults to the first model year.

- last_yr:

  the last year.

- n_proj:

  the number of projection years.

- n_iter:

  the number of simulated trajectories.

## Value

An array of projected selectivity by iteration, fishery, year, and age.
