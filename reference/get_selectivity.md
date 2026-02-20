# Compute selectivity-at-age from length-based selectivity curves

Defines selectivity as a parametric function of length (logistic or
double-normal), then converts to selectivity-at-age by
matrix-multiplying with the probability-of-length-at-age (PLA/ALK).

## Usage

``` r
get_selectivity(data, par_sel, pla, len_mid)
```

## Arguments

- data:

  A list containing model data. Required elements: n_fishery, n_year,
  n_age, sel_type_f.

- par_sel:

  Numeric matrix of dimensions `[n_fishery, 6]`. Each row is a real-line
  parameter vector. For logistic (sel_type_f == 1), only columns 1:2 are
  used. For double-normal (sel_type_f == 2), all 6 are used.

- pla:

  Numeric matrix of dimensions `[n_len, n_age]` containing
  probability-of-length-at-age (age-length key). Computed via
  [`get_pla()`](https://n-ducharmebarth-noaa.github.io/opal/reference/get_pla.md).

- len_mid:

  Numeric vector (length n_len) of length-bin midpoints.

## Value

3D array `sel_fya` of dimensions `[n_fishery, n_year, n_age]`.

## Details

Selectivity is currently time-invariant within each fishery (constant
across years).
