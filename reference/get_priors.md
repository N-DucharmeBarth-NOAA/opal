# Get priors

Don't include priors for recruitment deviates (par_rdev_y) or
selectivity (e.g., par_log_sel_1) here because they are dealt with in
the
[`get_recruitment_prior()`](https://n-ducharmebarth-noaa.github.io/opal/reference/get_recruitment_prior.md)
and `get_selectivity_prior()` functions.

## Usage

``` r
get_priors(parameters, data = NULL)
```

## Arguments

- parameters:

  A `list` specifying the parameters to be passed to `MakeADFun`. Can be
  generated using the
  [`get_parameters()`](https://n-ducharmebarth-noaa.github.io/opal/reference/get_parameters.md)
  function.

- data:

  A `list` of data inputs (optional). Used to retrieve prior center
  values for growth/variability parameters (e.g.,
  `data$prior_log_L1_mean`).

## Value

A `list` of priors.
