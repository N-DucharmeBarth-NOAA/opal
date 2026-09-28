# Plot catch

Plot catch by year and fishery.

## Usage

``` r
plot_catch(data, obj = NULL, plot_resid = FALSE)
```

## Arguments

- data:

  An `opal_obj`, legacy `opal_fit`, or model data list.

- obj:

  Optional for Opal objects. The AD object created using `MakeADFun`.

- plot_resid:

  Logical; plot catch residuals instead of observed and predicted catch.

## Value

A `ggplot2` object.
