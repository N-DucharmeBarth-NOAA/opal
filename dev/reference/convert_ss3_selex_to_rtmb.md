# Convert SS3 selectivity parameters to RTMB real-line parameterization

Takes SS3 natural-scale selectivity parameters and converts them to the
centered real-line parameterization used by the RTMB `sel_logistic` and
`sel_double_normal` functions.

## Usage

``` r
convert_ss3_selex_to_rtmb(
  ss3_pars,
  sel_type_f,
  sel_lengths,
  basis = "length",
  sel_ages = NULL
)
```

## Arguments

- ss3_pars:

  A data.frame or matrix with one row per fishery. For double-normal
  (pattern 24 or 20): columns are peak, top_logit, ascend_se,
  descend_se, start_logit, end_logit. For logistic (pattern 1 or 12):
  columns are inflection, width.

- sel_type_f:

  Integer vector (length n_fishery). 1 = logistic, 2 = double-normal, 3
  = double Richards.

- sel_lengths:

  Numeric vector of selectivity length-bin midpoints (same vector that
  will be passed to sel_logistic/sel_double_normal). Used for
  length-based fisheries.

- basis:

  Character vector of length 1 or n_fishery: `"length"` (default) or
  `"age"`, giving the selectivity basis of each fishery.

- sel_ages:

  Numeric vector of ages at which age-based curves are evaluated, i.e.
  `seq(min_age, by = 1, length.out = n_age)`. Required when any fishery
  has `basis = "age"`.

## Value

Numeric matrix `[n_fishery, 6]` of RTMB real-line parameters.

## Details

Length-based fisheries (`basis = "length"`) take SS3 `size_selex`
patterns 1 (logistic) and 24 (double normal), scaled by `sel_lengths`.
Age-based fisheries (`basis = "age"`, used with `sel_basis_f = 2`) take
SS3 `age_selex` patterns 12 (logistic) and 20 (double normal), which
share the parameter definitions of patterns 1 and 24, scaled by
`sel_ages`.

SS3's `-999` codes for the initial and final double-normal selectivity
("no adjustment") have no exact equivalent in opal. They are converted
to `-9`, which sets the initial or final selectivity close to zero. With
a final value of `-999`, SS3 leaves the descending limb unadjusted at
the oldest age or length bin, so the two curves can differ there.
