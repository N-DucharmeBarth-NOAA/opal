# Project recruitment deviates

Project recruitment deviates

## Usage

``` r
project_rec_devs(
  data,
  obj,
  mcmc = NULL,
  first_yr = 1931,
  last_yr = NULL,
  n_proj = 5,
  n_iter = NULL,
  max.p = 5,
  max.d = 5,
  max.q = 5,
  arima = TRUE
)
```

## Arguments

- data:

  a `list` of parameter values.

- obj:

  a `list` of parameter values.

- mcmc:

  a `list` of parameter values.

- first_yr:

  a `list` of inputs.

- last_yr:

  a `list` of inputs.

- n_proj:

  a `list` of inputs.

- n_iter:

  a `list` of inputs.

- max.p:

  Maximum value of p, or the maximum value of p (the AR order) to
  consider.

- max.d:

  Maximum value of d, or the maximum value of q (the MA order) to
  consider.

- max.q:

  Maximum value of q

- arima:

  default = TRUE, FALSE = "lognormal"

## Value

a `list` of projected recruitment deviates and ARIMA specifications.
