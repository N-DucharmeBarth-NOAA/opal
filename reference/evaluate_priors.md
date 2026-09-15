# Evaluate priors

Evaluate priors

## Usage

``` r
evaluate_priors(parameters, priors)
```

## Arguments

- parameters:

  A `list` specifying the parameters to be passed to `MakeADFun`. Can be
  generated using the
  [`get_parameters()`](https://n-ducharmebarth-noaa.github.io/opal/reference/get_parameters.md)
  function.

- priors:

  A `list` of named `list`s specifying priors for the parameters. Can be
  generated using the
  [`get_priors()`](https://n-ducharmebarth-noaa.github.io/opal/reference/get_priors.md)
  function.

## Value

A `numeric` value.

## Examples

``` r
{
  parameters <- list(log_B0 = log(1e6))
  priors <- list(
    log_B0 = list(
      type = "normal", par1 = log(1e6), par2 = 1,
      index = which("log_B0" == names(parameters))
    )
  )
  evaluate_priors(parameters, priors)
}
#> [1] 0.9189385
```
