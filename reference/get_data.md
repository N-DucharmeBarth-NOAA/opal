# Get bundled model data

Loads one of the packaged opal model data objects. This replaces the
legacy data-construction helper, which depended on historical raw inputs
that are no longer bundled with the package.

## Usage

``` r
get_data(
  model = c("opal_baseline", "opakapaka", "wcpo_bet"),
  include_parameters = FALSE
)
```

## Arguments

- model:

  Character model identifier. Supported values are `"opal_baseline"`,
  `"opakapaka"`, and `"wcpo_bet"`. Aliases `"baseline"`, `"opaka"`, and
  `"bet"` are also accepted.

- include_parameters:

  Logical; if `TRUE`, return a list with both `data` and matching
  initial `parameters`.

## Value

A data list ready for
[`opal_model`](https://n-ducharmebarth-noaa.github.io/opal/reference/opal_model.md),
or a list with elements `data` and `parameters` when
`include_parameters = TRUE`.
