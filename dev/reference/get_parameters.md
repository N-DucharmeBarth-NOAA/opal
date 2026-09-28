# Get bundled initial parameter values

Loads the packaged initial parameter list that matches one of the
bundled opal model data sets.

## Usage

``` r
get_parameters(data = NULL, model = NULL)
```

## Arguments

- data:

  Optional model data list used to infer the bundled parameter set.

- model:

  Optional character model identifier. Supported values are
  `"opal_baseline"`, `"opakapaka"`, and `"wcpo_bet"`. Aliases
  `"baseline"`, `"opaka"`, and `"bet"` are also accepted. Supplying
  `model` is preferred when `data` has been modified after loading.

## Value

A `list` of initial parameter values.
