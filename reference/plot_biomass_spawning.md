# Plot spawning biomass

Plot spawning biomass or relative spawning biomass by year for one or
more model runs.

## Usage

``` r
plot_biomass_spawning(data_list, object_list, relative = TRUE, labels = NULL)
```

## Arguments

- data_list:

  A list of model data lists passed to `MakeADFun`.

- object_list:

  A list of AD objects created using `MakeADFun`.

- relative:

  Logical; plot spawning biomass relative to unfished biomass.

- labels:

  Optional labels for the model runs.

## Value

A `ggplot2` object.
