# Plot spawning biomass

Plot spawning biomass or relative spawning biomass by year for one or
more model runs.

## Usage

``` r
plot_biomass_spawning(
  data_list,
  object_list = NULL,
  relative = TRUE,
  labels = NULL
)
```

## Arguments

- data_list:

  An Opal object, list of Opal objects, or list of model data lists.

- object_list:

  Omit for Opal objects. For legacy calls, a list of RTMB objectives.

- relative:

  Logical; plot spawning biomass relative to unfished biomass.

- labels:

  Optional labels for the model runs.

## Value

A `ggplot2` object.

## Examples

``` r
assessment <- opal_read(system.file("extdata", "opaka_quickstart_fit.rds",
                                   package = "opal"))
#> Warning: The saved fit runtime identity differs because it was serialized under a different R version; verifying by rebuild.
#> Warning: The saved fit runtime identity differs because it was serialized under a different R version; verifying by rebuild.
plot_biomass_spawning(assessment, relative = FALSE, labels = "Baseline")
```
