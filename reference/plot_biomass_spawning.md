# Plot spawning biomass

Plot the spawning biomass (tonnes) or relative spawning biomass by year
for a single model run, a grid of model runs, or an MCMC.

## Usage

``` r
plot_biomass_spawning(
  data_list,
  object_list,
  posterior = NULL,
  probs = c(0.025, 0.975),
  relative = TRUE,
  labels = NULL
)
```

## Arguments

- data_list:

  a `list` containing the data that was passed to `MakeADFun`.

- object_list:

  a `list` specifying the AD object created using the `MakeADFun`
  function.

- posterior:

  an `rstan` objected created using the `tmbstan` function.

- probs:

  a numeric vector of probabilities with values in `[0,1]` for plotting
  quantiles of the posterior distribution. Defaults to the 90% credible
  interval.

- relative:

  if the plot should be relative spawning biomass.

- labels:

  a `vector` of labels for the model runs.

## Value

a `ggplot2` object.
