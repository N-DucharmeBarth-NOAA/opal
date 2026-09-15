# Project dynamics

Forward-projects population dynamics for `n_proj` years using either
MCMC posterior draws (when `mcmc` is supplied) or multivariate-normal
(MVN) draws derived from the Hessian-based variance-covariance matrix at
the MLE (when `mcmc = NULL`).

## Usage

``` r
project_dynamics(
  data,
  object,
  mcmc = NULL,
  n_proj = 5,
  n_iter = 1,
  rdev_y,
  sel_fya,
  catch_ysf,
  return_hist = FALSE
)
```

## Arguments

- data:

  A `list` of model data (as passed to `opal_model`).

- object:

  The RTMB AD object returned by `RTMB::MakeADFun`, after optimisation.

- mcmc:

  Optional. MCMC fit object returned by `SparseNUTS`. When supplied,
  posterior draws are used for the projection. When `NULL` (default),
  MVN draws are generated from the Hessian-derived variance-covariance
  matrix.

- n_proj:

  Integer. Number of projection years.

- n_iter:

  Integer. Number of iterations (posterior draws or MVN samples).

- rdev_y:

  Numeric matrix `[n_iter, n_proj]`. Projected recruitment deviates
  (e.g., from
  [`project_rec_devs`](https://n-ducharmebarth-noaa.github.io/opal/reference/project_rec_devs.md)).

- sel_fya:

  Numeric array `[n_iter, n_fishery, n_proj, n_age]`. Projected
  selectivity (e.g., from
  [`project_selectivity`](https://n-ducharmebarth-noaa.github.io/opal/reference/project_selectivity.md)).

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

## Value

When `return_hist = FALSE` (default): a `list` of length `n_iter`, each
element being the named list returned by
[`do_dynamics`](https://n-ducharmebarth-noaa.github.io/opal/reference/do_dynamics.md)
for that iteration. When `return_hist = TRUE`: a named list with
elements `dyn` and `hist_sbio`.
