# Compute mean length-at-age using the Schnute parameterization of VB growth

Uses the Schnute parameterization of the von Bertalanffy growth curve,
matching the SS3 formulation with `CV_Growth_Pattern = 2`.

## Usage

``` r
get_growth(n_age, A1, A2, L1, L2, k)
```

## Arguments

- n_age:

  Integer. Number of age classes.

- A1:

  Integer. Reference age for L1 (data).

- A2:

  Integer. Reference age for L2 (data).

- L1:

  Numeric. Length at age A1 (may be AD).

- L2:

  Numeric. Length at age A2 (may be AD).

- k:

  Numeric. VB growth coefficient (may be AD).

## Value

Numeric vector of length `n_age`: mean length at each age
`a = 1, ..., n_age`.
