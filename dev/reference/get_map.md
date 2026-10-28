# Get default parameter mapping

Get a default parameter mapping. Parameter mapping is used by
`MakeADFun` to turn parameters on/off.

## Usage

``` r
get_map(parameters)
```

## Arguments

- parameters:

  a `list` containing the initial parameter values to be passed to
  `MakeADFun`.

## Value

a named `list` of parameter mapping.

## Details

Growth parameters `log_L1`, `log_L2`, `log_k`, `log_k_mult` and
`richards_b` are fixed. The growth cessation parameters (`log_L0`,
`log_rmax`, `log_k_gc`, `t50_gc`) are estimated when present. Under
`growth_model = "cessation"`, `log_L1`, `log_L2` and `log_k` are unused,
and fixing them keeps the Hessian non-singular.
