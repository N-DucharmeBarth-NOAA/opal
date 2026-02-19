# Natural mortality at age

Constructs a vector of M-at-age values using a declining early-age curve
and late-age increase.

## Usage

``` r
get_M(min_age, max_age, age_increase_M, m0, m4, m10, m30)
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

- m4:

  M at age-4 (controls slope of decline).

- m10:

  M at age-10 (base for flat zone).

- m30:

  M at age-30 (terminal age M).

## Value

vector of M at age.
