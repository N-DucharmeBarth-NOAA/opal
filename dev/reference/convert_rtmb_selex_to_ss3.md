# Convert RTMB real-line selectivity parameters back to SS3 natural scale

Inverse of `convert_ss3_selex_to_rtmb`. Useful for verifying round-trip
conversion and reporting parameter values in natural units.

## Usage

``` r
convert_rtmb_selex_to_ss3(
  par_sel,
  sel_type_f,
  sel_lengths,
  basis = "length",
  sel_ages = NULL
)
```

## Arguments

- par_sel:

  Numeric matrix `[n_fishery, 6]` of RTMB real-line parameters.

- sel_type_f:

  Integer vector (1 = logistic, 2 = double-normal, 3 = double Richards).

- sel_lengths:

  Numeric vector of selectivity length-bin midpoints. Used for
  length-based fisheries.

- basis:

  Character vector of length 1 or n_fishery: `"length"` (default) or
  `"age"`, giving the selectivity basis of each fishery.

- sel_ages:

  Numeric vector of ages at which age-based curves are evaluated, i.e.
  `seq(min_age, by = 1, length.out = n_age)`. Required when any fishery
  has `basis = "age"`.

## Value

Data.frame with SS3-scale parameter values per fishery. Double Richards
rows retain `NA` parameter values because there is no SS3 equivalent.
Age-based rows correspond to SS3 `age_selex` patterns 12 and 20.
