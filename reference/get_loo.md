# Obtain the leave-one-out information criterion (LOO IC)

This function extracts point-wise log-likelihood values and then
calculates LOO statistics.

## Usage

``` r
get_loo(data, object, posterior)
```

## Arguments

- data:

  a `list` containing the data that was passed to `MakeADFun`.

- object:

  A `list` specifying the AD object created using `MakeADFun`.

- posterior:

  An `snutsfit` object created using the `sample_snuts` function.

## Value

a `psis_loo` object.
