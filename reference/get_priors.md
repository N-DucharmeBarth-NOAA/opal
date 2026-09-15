# Get priors

Get priors

## Usage

``` r
get_priors(parameters, data = NULL)
```

## Arguments

- parameters:

  A `list` specifying the parameters to be passed to `MakeADFun`. Can be
  generated using the
  [`get_parameters()`](https://n-ducharmebarth-noaa.github.io/opal/reference/get_parameters.md)
  function. Vector parameters (e.g., `log_cpue_q`) are supported.

- data:

  A `list` of data inputs (optional). Used to retrieve prior center
  values for growth/variability parameters (e.g.,
  `data$prior_log_L1_mean`).

## Value

A `list` of priors.
