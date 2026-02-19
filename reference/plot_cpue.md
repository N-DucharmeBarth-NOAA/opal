# Plot CPUE

Plot CPUE fit.

## Usage

``` r
plot_cpue(data, object, posterior = NULL, probs = c(0.025, 0.975), nsim = 10)
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

- nsim:

  number of simulations to plot.

## Value

a `ggplot2` object.
