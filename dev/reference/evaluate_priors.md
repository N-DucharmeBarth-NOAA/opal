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
  [`get_parameters()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/get_parameters.md)
  function.

- priors:

  A `list` of named `list`s specifying priors for the parameters. Can be
  generated using the
  [`get_priors()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/get_priors.md)
  function.

## Value

A `numeric` value.

## Details

Supported distributions are `normal`, `student` (three degrees of
freedom), `lognormal`, and `beta` (mean and precision). Priors act on
the stored parameter scale; for example, a normal prior on `log_B0` is a
normal density on log spawning output. Locations and scales must be
scalars or match the parameter block length. Indices identify blocks in
`parameters`.

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
