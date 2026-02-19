# Convert SS3 selectivity parameters to RTMB real-line parameterization

Takes SS3 natural-scale selectivity parameters and converts them to the
centered real-line parameterization used by the RTMB `sel_logistic` and
`sel_double_normal` functions.

## Usage

``` r
convert_ss3_selex_to_rtmb(ss3_pars, sel_type_f, sel_lengths)
```

## Arguments

- ss3_pars:

  A data.frame or matrix with one row per fishery. For double-normal
  (pattern 24): columns are peak, top_logit, ascend_se, descend_se,
  start_logit, end_logit. For logistic (pattern 1): columns are
  inflection, width.

- sel_type_f:

  Integer vector (length n_fishery). 1 = logistic, 2 = double-normal.

- sel_lengths:

  Numeric vector of selectivity length-bin midpoints (same vector that
  will be passed to sel_logistic/sel_double_normal).

## Value

Numeric matrix `[n_fishery, 6]` of RTMB real-line parameters.
