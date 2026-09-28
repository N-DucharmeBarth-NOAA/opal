# Plot CPUE

Plot observed and predicted CPUE by season and fishery.

## Usage

``` r
plot_cpue(data, object = NULL)
```

## Arguments

- data:

  An `opal_obj`, legacy `opal_fit`, or model data list.

- object:

  Optional for Opal objects. The AD object created using `MakeADFun`.

## Value

A `ggplot2` object.

## Examples

``` r
assessment <- opal_read(system.file("extdata", "opaka_quickstart_fit.rds",
                                   package = "opal"))
#> Warning: The saved fit runtime identity differs because it was serialized under a different R version; verifying by rebuild.
#> Warning: The saved fit runtime identity differs because it was serialized under a different R version; verifying by rebuild.
plot_cpue(assessment)
#> Warning: Unknown or uninitialised column: `season`.
```
