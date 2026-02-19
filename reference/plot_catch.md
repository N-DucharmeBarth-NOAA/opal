# Plot catch

Plot catch (in thousands of tonnes) by year, season, and fishery.

## Usage

``` r
plot_catch(
  data,
  obj,
  posterior = NULL,
  proj = NULL,
  probs = c(0.05, 0.95),
  plot_resid = FALSE
)
```

## Arguments

- data:

  a `list` containing the data that was passed to `MakeADFun`.

- obj:

  a `list` specifying the AD object created using `MakeADFun`.

- posterior:

  an `rstan` objected created using the `tmbstan` function.

- proj:

  an `rstan` objected created using the `tmbstan` function.

- probs:

  a numeric vector of probabilities with values in `[0,1]` for plotting
  quantiles of the posterior distribution.

- plot_resid:

  plot the residual or not.

## Value

a `ggplot2` object.
