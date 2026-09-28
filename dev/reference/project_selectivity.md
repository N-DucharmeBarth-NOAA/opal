# Project selectivity

Generate future selectivity from a selected historical period in the
object's point report. For each age and fishery, variable positive
selectivities are sampled independently on the log scale; constant
values are repeated. This does not propagate posterior uncertainty in
selectivity parameters.

## Usage

``` r
project_selectivity(
  data,
  obj = NULL,
  first_yr = NULL,
  last_yr = NULL,
  n_proj = 5,
  n_iter = 1
)
```

## Arguments

- data:

  An Opal object or model data list.

- obj:

  For legacy data-list calls, the RTMB objective. Omit for an Opal
  object.

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

## See also

[`opal_project()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_project.md),
[`project_dynamics()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/project_dynamics.md)
