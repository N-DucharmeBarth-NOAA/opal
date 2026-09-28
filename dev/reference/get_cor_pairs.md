# Return strongly-correlated parameter pairs

Computes the parameter correlation matrix from the (inverse) Hessian and
returns a data.frame of pairs whose absolute correlation exceeds
`threshold`. Pairs are sorted by descending absolute correlation. If no
pairs meet the criterion the function returns an informative string
instead of a data.frame.

## Usage

``` r
get_cor_pairs(
  obj,
  h = NULL,
  threshold = 0.9,
  include = c("core", "all"),
  digits = 3
)
```

## Arguments

- obj:

  An RTMB AD object from `MakeADFun` after optimisation.

- h:

  Optional pre-computed Hessian matrix (square, same dimension as
  `length(obj$par)`). When `NULL` (default) the Hessian is obtained via
  `obj$he()` if available, otherwise
  [`optimHess`](https://rdrr.io/r/stats/optim.html).

- threshold:

  Absolute correlation threshold. Only pairs with `|cor| > threshold`
  are returned. Default `0.9`.

- include:

  Which parameters to consider: `"core"` (default) excludes `rdev_y`
  parameters; `"all"` includes everything in `obj$par`.

- digits:

  Number of significant digits applied to the correlation column.
  Default `3`. Use `NULL` to suppress rounding.

## Value

A data.frame with columns `par1`, `par2`, `correlation` (sorted by
descending `|correlation|`), or a character string when no pairs exceed
`threshold`.
