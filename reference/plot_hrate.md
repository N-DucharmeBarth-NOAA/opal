# Plot harvest rate

Plot harvest rate by year, season, and age.

## Usage

``` r
plot_hrate(
  data,
  object,
  posterior = NULL,
  probs = c(0.025, 0.975),
  years = 2013:2022,
  ...
)
```

## Arguments

- data:

  a `list` containing the data that was passed to `MakeADFun`.

- object:

  a `list` specifying the AD object created using `MakeADFun`.

- posterior:

  an `rstan` objected created using the `tmbstan` function.

- probs:

  a numeric vector of probabilities with values in `[0,1]` for plotting
  quantiles of the posterior distribution.

- years:

  the years to show on the plot.

- ...:

  options passed on to `geom_density_ridges`.

## Value

a `ggplot2` object.
