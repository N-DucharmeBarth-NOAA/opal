# Plot LF

Plot LF fit.

## Usage

``` r
plot_lf(
  data,
  object,
  posterior = NULL,
  probs = c(0.025, 0.975),
  years = NULL,
  fishery = "LL1",
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

  the years to plot.

- fishery:

  the fishery to plot (LL1, LL2, LL3, or LL4).

- ...:

  additional parameters passed on to `facet_wrap`.

## Value

a `ggplot2` object.
