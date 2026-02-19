# Convert RTMB real-line selectivity parameters back to SS3 natural scale

Inverse of `convert_ss3_selex_to_rtmb`. Useful for verifying round-trip
conversion and reporting parameter values in natural units.

## Usage

``` r
convert_rtmb_selex_to_ss3(par_sel, sel_type_f, sel_lengths)
```

## Arguments

- par_sel:

  Numeric matrix `[n_fishery, 6]` of RTMB real-line parameters.

- sel_type_f:

  Integer vector (1 = logistic, 2 = double-normal).

- sel_lengths:

  Numeric vector of selectivity length-bin midpoints.

## Value

Data.frame with SS3-scale parameter values per fishery.
