# Project recruitment deviates

Project recruitment deviates

## Usage

``` r
project_selectivity(
  data,
  obj,
  mcmc = NULL,
  first_yr = 2000,
  last_yr = NULL,
  n_proj = 5,
  n_iter = NULL,
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

  the first year.

- last_yr:

  the last year.

- n_proj:

  the number of projection years.

- n_iter:

  a `list` of inputs.

- arima:

  default = TRUE, FALSE = "lognormal"

## Value

the negative log-likelihood (NLL) value.
