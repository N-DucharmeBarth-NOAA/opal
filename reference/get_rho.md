# Estimate temporal autocorrelation in recruitment deviations

Calculates the AR1 autocorrelation coefficient (phi) for recruitment
deviations.

## Usage

``` r
get_rho(first_yr = 1931, last_yr = 2022, rdev)
```

## Arguments

- first_yr:

  First model year.

- last_yr:

  Last model year.

- rdev:

  Vector of recruitment deviations.

## Value

Estimated autocorrelation.

## Examples

``` r
first_yr <- 1931
last_yr <- 2022
N <- length(first_yr:last_yr)
rdev <- arima.sim(list(order = c(1, 0, 0), ar = 0.5), n = N)
get_rho(first_yr, last_yr, rdev)
#> Error in `[.default`(rdev, i1:(i2 - 1)): only 0's may be mixed with negative subscripts
```
