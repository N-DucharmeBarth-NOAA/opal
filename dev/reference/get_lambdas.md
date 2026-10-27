# Resolve likelihood-component lambdas

Expands the optional `data$lambdas` list to one multiplier per CPUE
index or composition fleet. Scalar values are recycled; unspecified
components use 1.

## Usage

``` r
get_lambdas(data)
```

## Arguments

- data:

  A model data list.

## Value

A fully expanded named list of lambda values.

## Details

The supported entries are `cpue` (one or `n_index` values), `lf` and
`wf` (one or `n_fishery` values each), and scalar `rec`, `init_rec`,
`prior`, and `penalty` multipliers. `prior_mapped_only` is a logical
flag that removes prior contributions for elements fixed by the
parameter map.

[`opal_model()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_model.md)
cannot see the parameter map, so `prior_mapped_only` is applied only
when the objective is built from an `opal_obj` (through
[`opal_build()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_build.md),
[`opal_rtmb()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_build.md),
[`opal_fit()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_fit.md)
or
[`opal_mcmc()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_mcmc.md)).
An objective built directly with
`RTMB::MakeADFun(cmb(opal_model, data), ...)`, or rebuilt from a legacy
`opal_fit`, still includes priors on mapped elements, and
[`opal_model()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_model.md)
warns that the flag was not applied.
