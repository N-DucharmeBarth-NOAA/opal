# Probability of length at age matrix (age-length key)

Computes the probability distribution of fish lengths across different
age classes. AD-compatible: uses RTMB::pnorm and ADoverload so that
gradients propagate if growth parameters (mu_a, sd_a) are estimated.

## Usage

``` r
get_pla(len_lower, len_upper, mu_a, sd_a)
```

## Arguments

- len_lower:

  Numeric vector of lower bounds of length bins (length L). Data only.

- len_upper:

  Numeric vector of upper bounds of length bins (length L). Data only.

- mu_a:

  Numeric vector of mean length at age (length A). May be AD.

- sd_a:

  Numeric vector of SD of length at age (length A). May be AD.

## Value

Matrix of dimensions L x A where columns sum to 1.
