# Convert maturity-at-length to maturity-at-age

Uses the probability-of-length-at-age matrix (PLA) to convert a
maturity-at-length vector to maturity-at-age: `mat_a = t(pla) %*% mat_l`

## Usage

``` r
get_maturity_at_age(pla, maturity_at_length)
```

## Arguments

- pla:

  Matrix (n_len x n_age). Probability of length at age (from
  [`get_pla`](https://n-ducharmebarth-noaa.github.io/opal/reference/get_pla.md)).

- maturity_at_length:

  Numeric vector (length n_len). Maturity at each length bin (data).

## Value

Numeric vector (length n_age). Maturity at each age.
