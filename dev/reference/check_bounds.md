# Check if parameters are up against the bounds

A `data.frame` containing the parameters that are against their lower or
upper bound.

## Usage

``` r
check_bounds(opt, lower, upper)
```

## Arguments

- opt:

  an optimized TMB object.

- lower:

  a vector of lower bounds.

- upper:

  a vector of upper bounds.

## Value

a `data.frame`.
