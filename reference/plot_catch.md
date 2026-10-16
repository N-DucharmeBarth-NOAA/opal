# Plot catch

Plot catch by year and fishery.

## Usage

``` r
plot_catch(data, obj = NULL, plot_resid = FALSE, weight_units = NULL)
```

## Arguments

- data:

  An `opal_obj`, legacy `opal_fit`, or model data list.

- obj:

  Optional for Opal objects. The AD object created using `MakeADFun`.

- plot_resid:

  Logical; plot catch residuals instead of observed and predicted catch.

- weight_units:

  Label for weight catches, in the input data's units. Defaults to
  `data$catch_weight_units`, or `"weight units"` when unspecified.

## Value

A `ggplot2` object.
