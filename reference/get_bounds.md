# Get default parameter bounds

Get `data.frame` of default parameter bounds.

## Usage

``` r
get_bounds(obj, parameters)
```

## Arguments

- obj:

  a `list` specifying the AD object created using the `MakeADFun`
  function.

- parameters:

  The parameter list used to construct `obj`. Retained for API
  compatibility.

## Value

a `data.frame` of parameter bounds.
