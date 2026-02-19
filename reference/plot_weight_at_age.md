# Plot weight at age

Plot the mean weight (kg) at age by fishery for a subset of years.

## Usage

``` r
plot_weight_at_age(data, years = c(1931, 2000, 2010, 2022))
```

## Arguments

- data:

  a `list` containing the data that was passed to `MakeADFun`.

- years:

  the years to show on the plot.

## Value

a `ggplot2` object.
