# Project dynamics

Forward-projects population dynamics for `n_proj` years using either
stored MCMC posterior draws (`uncertainty = "mcmc"`) or
multivariate-normal draws from a fitted fixed-effects model
(`uncertainty = "mvn"`). Pass the assessment as `data`; its runtime and
posterior are used internally. Use
[`opal_project()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_project.md)
to retain the result and its provenance in the object.

## Usage

``` r
project_dynamics(
  data,
  object = NULL,
  mcmc = NULL,
  n_proj = 5,
  n_iter = 1,
  rdev_y,
  sel_fya,
  catch_ysf,
  return_hist = FALSE,
  uncertainty = NULL
)
```

## Arguments

- data:

  An `opal_obj`, legacy `opal_fit`, or model data list.

- object:

  For legacy data-list calls, the fitted RTMB objective. Omit when
  `data` is an Opal object.

- mcmc:

  Legacy data-list calls only: a SparseNUTS fit. Omit for Opal objects,
  which use their stored posterior. In the legacy interface, when
  supplied, posterior draws are used for the projection. When `NULL`
  (default), MVN draws are generated from the Hessian-derived
  variance-covariance matrix.

- n_proj:

  Integer. Number of projection years.

- n_iter:

  Integer. Number of iterations (posterior draws or MVN samples).

- rdev_y:

  Numeric matrix `[n_iter, n_proj]`. Projected recruitment deviates
  (e.g., from
  [`project_rec_devs`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/project_rec_devs.md)).

- sel_fya:

  Numeric array `[n_iter, n_fishery, n_proj, n_age]`. Projected
  selectivity (e.g., from
  [`project_selectivity`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/project_selectivity.md)).

- catch_ysf:

  Numeric array `[n_proj, n_season, n_fishery]`. Projected observed
  catch by year, season, and fishery.

- return_hist:

  Logical (default `FALSE`). When `TRUE` the function returns a named
  list with elements `dyn` (the projection results) and `hist_sbio`, a
  numeric matrix `[n_iter, n_year + 1]` containing the full
  spawning-biomass trajectory from each parameter draw over the
  historical period. The last column of `hist_sbio` corresponds to the
  same terminal state as the `proj_year = 0` bridge point in the
  projection, so the two can be plotted seamlessly without any join
  discontinuity.

- uncertainty:

  For Opal objects, choose `mcmc` or `mvn` explicitly when a posterior
  exists. Random effects require complete joint posterior draws.

## Value

When `return_hist = FALSE` (default): a `list` of length `n_iter`, each
element being the named list returned by
[`do_dynamics`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/do_dynamics.md)
for that iteration. When `return_hist = TRUE`: a named list with
elements `dyn` and `hist_sbio`.
