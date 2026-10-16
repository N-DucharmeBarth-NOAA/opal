# Plot harvest rate

Plot harvest rate by year, season, and age.

## Usage

``` r
plot_hrate(data, object = NULL, years = NULL, ...)
```

## Arguments

- data:

  An `opal_obj`, legacy `opal_fit`, or model data list.

- object:

  Optional for Opal objects. The AD object created using `MakeADFun`.

- years:

  Optional years to show. The default plots every model year from the
  first catch year onwards.

- ...:

  Options passed to `geom_density_ridges`.

## Value

A `ggplot2` object.
