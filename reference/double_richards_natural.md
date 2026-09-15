# Convert double Richards parameters to natural scale

This numeric reporting helper returns the transformed parameter values
and Maunder's `alpha1`, `beta1`, `alpha2`, and `beta2`. Under the
50%-point parameterization, `alpha2` can be non-positive when the shape
parameters differ.

## Usage

``` r
double_richards_natural(len, par)
```

## Arguments

- len:

  Numeric vector of length-bin midpoints.

- par:

  Numeric vector of length 6 containing selectivity parameters.

## Value

Named numeric vector of natural-scale parameters.
