# Plot catch

Plot catch by year and fishery.

## Usage

``` r
plot_catch(data, obj, plot_resid = FALSE)
```

## Arguments

- data:

  A model data list passed to `MakeADFun`.

- obj:

  The AD object created using `MakeADFun`.

- plot_resid:

  Logical; plot catch residuals instead of observed and predicted catch.

## Value

A `ggplot2` object.
