# Calculate Recruitment Bias Adjustment Ramp

Calculate Recruitment Bias Adjustment Ramp

## Usage

``` r
get_bias_adj_vector(years, do_rec_bias_ramp, bias_years, max_bias_adj)
```

## Arguments

- years:

  Numeric vector of model years.

- do_rec_bias_ramp:

  Integer flag (0 = off, 1 = on).

- bias_years:

  Numeric vector of length 4 defining the ramp (ascend start, plateau
  start, plateau end, descend end).

- max_bias_adj:

  Numeric scalar for the maximum bias adjustment fraction.

## Value

A numeric vector of length `length(years)` containing the bias
adjustment scalars.
