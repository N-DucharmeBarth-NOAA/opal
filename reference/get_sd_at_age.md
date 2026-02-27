# Compute SD of length-at-age from CV1 and CV2

Linearly interpolates CV as a function of mean length, matching SS3's
`CV_Growth_Pattern = 2`.

## Usage

``` r
get_sd_at_age(mu_a, L1, L2, log_CV1, log_CV2)
```

## Arguments

- mu_a:

  Numeric vector. Mean length at age (from
  [`get_growth`](https://n-ducharmebarth-noaa.github.io/opal/reference/get_growth.md)).
  May be AD.

- L1:

  Numeric. Length at age A1 (may be AD).

- L2:

  Numeric. Length at age A2 (may be AD).

- log_CV1:

  Numeric. CV at age A1 (may be AD).

- log_CV2:

  Numeric. CV at age A2 (may be AD).

## Value

Numeric vector of SD at each age.
