# Double Richards selectivity as a function of length

A smooth, six-parameter curve formed by multiplying an ascending
generalized-logistic (Richards) limb by one minus a descending Richards
limb: \$\$(1 + \exp(-z_1))^{-1/\gamma_1} \[1 - (1 +
\exp(-z_2))^{-1/\gamma_2}\].\$\$ Parameters are on the real line and
transformed as follows:

- `par[1]`:

  Ascending 50% point: `mean(len) + par[1] * sd(len)`.

- `par[2]`:

  Ascending width: `exp(par[2]) * sd(len)`.

- `par[3]`:

  Ascending shape: `gamma1 = exp(par[3])`.

- `par[4]`:

  Descending 50% point is ascending 50% point plus
  `exp(par[4]) * sd(len)`.

- `par[5]`:

  Descending width: `exp(par[5]) * sd(len)`.

- `par[6]`:

  Descending shape: `gamma2 = exp(par[6])`.

The Maunder parameterization is recovered with
`alpha_k = L50_k + log(2^gamma_k - 1) / beta_k` and
`beta_k = log(19) / width_k`. The peak can be below one when the limbs
overlap; it is intentionally not rescaled. As either shape parameter
tends to negative infinity, its limb tends to a Gompertz curve with the
same 50% point. Large shape values flatten the limb toward 0.5 (at
`par[3] = par[6] = 7`, the curve is approximately 0.25). Mapping
`par[3]` and `par[6]` to zero gives a four-parameter double-logistic
curve. Suggested bounds for columns 3 and 6 are `c(-5, 5)`; use the
existing `par_sel` bounds for the remaining columns.

## Usage

``` r
sel_double_richards(len, par)
```

## Arguments

- len:

  Numeric vector of length-bin midpoints.

- par:

  Numeric vector of length 6 containing selectivity parameters.

## Value

Numeric vector of selectivity values in \[0, 1).

## References

Maunder, M. (2025). Double Richards selectivity. Unpublished technical
note, 19 November 2025.
