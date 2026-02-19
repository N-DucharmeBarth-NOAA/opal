# Compute selectivity-at-age from length-based selectivity curves

Defines selectivity as a parametric function of length (logistic or
double-normal), then converts to selectivity-at-age by
matrix-multiplying with the probability-of-length-at-age (PLA/ALK).

## Usage

``` r
get_selectivity(data, par_sel, mu_a, sd_a, len_lower, len_upper, len_mid)
```

## Arguments

- data:

  A list containing model data. Required elements: n_fishery, n_year,
  n_age, sel_type_f.

- par_sel:

  Numeric matrix of dimensions `[n_fishery, 6]`. Each row is a real-line
  parameter vector. For logistic (sel_type_f == 1), only columns 1:2 are
  used. For double-normal (sel_type_f == 2), all 6 are used.

- mu_a:

  Numeric vector (length n_age) of mean length at age. May be AD if
  growth parameters are estimated.

- sd_a:

  Numeric vector (length n_age) of SD of length at age. May be AD if
  growth parameters are estimated.

- len_lower:

  Numeric vector (length n_len) of lower bounds of length bins. Derived
  internally from `len_bin_start`, `len_bin_width`, and `n_len` inside
  [`bet_model()`](https://n-ducharmebarth-noaa.github.io/opal/reference/bet_model.md).

- len_upper:

  Numeric vector (length n_len) of upper bounds of length bins.

- len_mid:

  Numeric vector (length n_len) of length-bin midpoints.

## Value

3D array `sel_fya` of dimensions `[n_fishery, n_year, n_age]`.

## Details

The PLA is computed internally from `mu_a` and `sd_a` so that AD
gradients propagate correctly if growth parameters are estimated in the
future.

Selectivity is currently time-invariant within each fishery (constant
across years).
