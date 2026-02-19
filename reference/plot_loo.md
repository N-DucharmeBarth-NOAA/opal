# Plot PSIS LOO

This is the same as `plot(loo1)` except it colours/groups the different
types of data.

## Usage

``` r
plot_loo(x, exclude = NULL)
```

## Arguments

- x:

  A `psis_loo` object.

- exclude:

  Any variables to exclude from the plot (e.g., aerial, af, lf, cpue).

## Value

a `ggplot2`.
