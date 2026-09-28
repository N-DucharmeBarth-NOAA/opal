# Project recruitment deviates

Opal-object inputs reconstruct the full historical recruitment vector,
including fixed and mapped years, for each draw. Set first_yr and
last_yr to choose the period used to estimate future variability.
Constant histories produce constant projections. Legacy list inputs
retain their original path.

## Usage

``` r
project_rec_devs(
  data,
  obj = NULL,
  mcmc = NULL,
  first_yr = NULL,
  last_yr = NULL,
  n_proj = 5,
  n_iter = NULL,
  max.p = 5,
  max.d = 5,
  max.q = 5,
  arima = TRUE,
  uncertainty = NULL
)
```

## Arguments

- data:

  An Opal object or model data list.

- obj:

  For legacy data-list calls, the RTMB objective. Omit for an Opal
  object.

- mcmc:

  Legacy data-list calls only: a SparseNUTS fit. Opal objects use their
  stored posterior when `uncertainty = "mcmc"`.

- first_yr:

  the first year sampled. Defaults to the first model year.

- last_yr:

  Last historical year used to estimate future variability. Defaults to
  the final model year.

- n_proj:

  Number of future model years.

- n_iter:

  Number of trajectories. Defaults to one for fitted values, or all
  retained draws for MCMC. MCMC uses the first `n_iter` draws in the
  same order as
  [`project_dynamics()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/project_dynamics.md).

- max.p:

  Maximum value of p, or the maximum value of p (the AR order) to
  consider.

- max.d:

  Maximum differencing order considered by
  [`forecast::auto.arima()`](https://pkg.robjhyndman.com/forecast/reference/auto.arima.html).

- max.q:

  Maximum moving-average order considered by
  [`forecast::auto.arima()`](https://pkg.robjhyndman.com/forecast/reference/auto.arima.html).

- arima:

  If `TRUE`, select an ARIMA model and bootstrap future innovations. If
  `FALSE`, fit an autoregressive model with
  [`stats::ar()`](https://rdrr.io/r/stats/ar.html) and simulate with
  [`stats::arima.sim()`](https://rdrr.io/r/stats/arima.sim.html).

- uncertainty:

  For Opal objects, choose `fit` or `mcmc`.

## Value

a `list` of projected recruitment deviates and ARIMA specifications.

## See also

[`opal_project()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_project.md),
[`project_dynamics()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/project_dynamics.md)
