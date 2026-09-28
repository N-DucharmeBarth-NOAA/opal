# Summarise model parameters in a table

Builds a data.frame combining initial values, estimated values, bounds,
gradient, gradient check, and bounds check for all parameters. Each row
is one scalar element of the unlisted parameter vector.

## Usage

``` r
get_par_table(
  obj,
  parameters,
  map,
  lower = NULL,
  upper = NULL,
  grad_tol = 1e-04,
  bnds_tol = 0.025,
  include = c("core", "all_est", "all"),
  digits = 3,
  show_map = FALSE
)
```

## Arguments

- obj:

  An RTMB AD object from `MakeADFun`.

- parameters:

  Named list of initial parameter values passed to `MakeADFun`.

- map:

  Named list of map factors passed to `MakeADFun`.

- lower:

  Numeric vector of lower bounds (length `length(obj$par)`), or `NULL`
  for `-Inf` everywhere.

- upper:

  Numeric vector of upper bounds (length `length(obj$par)`), or `NULL`
  for `Inf` everywhere.

- grad_tol:

  Threshold for gradient check: absolute gradient below this is `"OK"`,
  otherwise `"BAD"`. Default `1e-4`.

- bnds_tol:

  Fraction of finite bound range defining the proximity margin for the
  bounds check. Default `0.025` (2.5%).

- include:

  Which rows to return: `"core"` (default) returns estimated parameters
  excluding `rdev_y`; `"all_est"` includes `rdev_y`; `"all"` includes
  fixed parameters too.

- digits:

  Number of significant digits for numeric columns. Default `3`. Use
  `NULL` to suppress rounding.

- show_map:

  Logical. Include `map` and `fixed` columns. Default `FALSE`.

## Value

A data.frame with columns `par`, `init`, `est`, `lwr`, `upr`, `grd`,
`gr_chk`, `bd_chk`, and optionally `group`, `map` and `fixed` (when
`show_map = TRUE`).
