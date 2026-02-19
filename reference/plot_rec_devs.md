# Plot recruitment deviates

Plot recruitment deviates by year.

## Usage

``` r
plot_rec_devs(
  data,
  object,
  proj = NULL,
  posterior = NULL,
  probs = c(0.025, 0.975)
)
```

## Arguments

- data:

  a `list` containing the data that was passed to `MakeADFun`.

- object:

  a `list` specifying the AD object created using `MakeADFun`.

- proj:

  a `list` containing the data that was passed to `MakeADFun`.

- posterior:

  an `rstan` objected created using the `tmbstan` function.

- probs:

  a numeric vector of probabilities with values in `[0,1]` for plotting
  quantiles of the posterior distribution.

## Value

a `ggplot2` object.
