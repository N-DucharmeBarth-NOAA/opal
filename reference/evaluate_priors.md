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
if (FALSE) { # \dontrun{
  parameters <- list(par_log_m4 = log(0.167))
  priors <- list(
    par_log_m4 = list(type = "normal", par1 = log(0.12), par2 = 0.4, 
                      index = which("par_log_m4" == names(parameters)))
  )
  evaluate_priors(parameters, priors)
} # }
```
