# Natural mortality at length

Constructs a vector of M-at-age values using a declining early-age curve
and late-age increase.

## Usage

``` r
get_M_length(min_age, max_age, age_increase_M, m0, m30, length_mu_ysa)
```

## Arguments

- min_age:

  minimum model age.

- max_age:

  maximum model age.

- age_increase_M:

  age at which M begins to increase again.

- m0:

  M at age-1 and 2.

- m30:

  M at age-30 (terminal age M).

- length_mu_ysa:

  M at age-30 (terminal age M).

## Value

vector of M at age.
