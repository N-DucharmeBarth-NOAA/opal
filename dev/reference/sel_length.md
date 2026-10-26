# Selectivity-at-length for a single fishery, dispatched on type code

Selectivity-at-length for a single fishery, dispatched on type code

## Usage

``` r
sel_length(len, par, sel_type)
```

## Arguments

- len:

  Numeric vector of length-bin midpoints, or of ages for an age-based
  fishery (`sel_basis_f = 2`); see
  [`get_selectivity()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/get_selectivity.md).

- par:

  Numeric vector of length 6 containing selectivity parameters.

- sel_type:

  Selectivity type: 1 = logistic, 2 = double-normal, 3 = double
  Richards.

## Value

Numeric vector of selectivity values.
