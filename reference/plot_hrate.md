# Plot harvest rate

Plot harvest rate by year, season, and age.

## Usage

``` r
plot_hrate(data, object, years = NULL, ...)
```

## Arguments

- data:

  A model data list passed to `MakeADFun`.

- object:

  The AD object created using `MakeADFun`.

- years:

  Optional years to show. The default plots every model year from the
  first catch year onwards.

- ...:

  Options passed to `geom_density_ridges`.

## Value

A `ggplot2` object.
