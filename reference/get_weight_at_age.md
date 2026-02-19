# Obtain the mean weight at age

Obtain the mean weight (kg) at age for each fishery each year. This
function is used within `get_data` so is generally not needed directly.

## Usage

``` r
get_weight_at_age(length_mu_ysa, length_sd_a)
```

## Arguments

- length_mu_ysa:

  an `array` containing the mean length at age (a) for each year (y) and
  season (s). This `array` can be generated using the function
  `get_length_at_age`.

- length_sd_a:

  a `vector` containing the standard deviation of mean length at age.

## Value

an `array`.

## Examples

``` r
length_mu_ysa <- get_length_at_age(length_mean = opal::length_mean)
#> Error: 'length_mean' is not an exported object from 'namespace:opal'
weight_fya <- get_weight_at_age(length_mu_ysa = length_mu_ysa, 
                                length_sd_a = opal::length_sd$SD)
#> Error: object 'length_mu_ysa' not found
```
