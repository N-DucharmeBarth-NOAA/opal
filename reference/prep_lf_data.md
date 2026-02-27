# Prepare length composition data for model input

Transforms a pivoted wide-format length-frequency data frame (one row
per fishery x timestep, bins as columns) into the set of model-ready
arrays and vectors that `opal_model` and `get_length_like` expect. The
function validates bin alignment against the length structure in `data`,
optionally filters to a subset of fisheries, and appends every derived
object directly to `data`.

## Usage

``` r
prep_lf_data(
  data,
  lf_wide,
  lf_keep_fisheries = NULL,
  lf_switch = 1L,
  lf_minbin = NULL,
  lf_maxbin = NULL,
  lf_var_adjust = NULL,
  lf_cap = NULL
)
```

## Arguments

- data:

  a `list` containing at minimum the length structure scalars
  `len_bin_start`, `len_bin_width`, `n_len`, and `n_fishery`.

- lf_wide:

  a `data.frame` in wide format: columns `fishery`, `year`, `month`,
  `ts` plus one numeric column per length bin (named by the bin
  lower-bound value). Typically produced by
  [`tidyr::pivot_wider()`](https://tidyr.tidyverse.org/reference/pivot_wider.html)
  on the long-format `wcpo_bet_lf` object.

- lf_keep_fisheries:

  integer vector of fishery indices to retain. Pass `NULL` (the default)
  to keep all fisheries present in `lf_wide`.

- lf_switch:

  integer likelihood type selector passed straight through to
  `data$lf_switch` (1 = multinomial, 2 = Dirichlet, 3 =
  Dirichlet-multinomial). Default `1L`.

- lf_minbin:

  integer vector length `n_fishery` giving the minimum bin index
  (1-based) used for each fishery. Defaults to `rep(1L, n_fishery)`.

- lf_maxbin:

  integer vector length `n_fishery` giving the maximum bin index
  (1-based) used for each fishery. Defaults to `rep(n_len, n_fishery)`.

- lf_var_adjust:

  numeric vector of length `n_fishery` used to scale the effective
  sample size for each fishery. For observation row \\i\\ belonging to
  fishery \\f\\, the sample size is adjusted as \\n_i /
  \text{lf\\var\\adjust}\_f\\. Values greater than 1 downweight (reduce
  effective \\N\\) and values less than 1 upweight the fishery. Defaults
  to `rep(1, n_fishery)` (no adjustment).

- lf_cap:

  positive integer (or `NULL`). When supplied, each effective sample
  size is capped: \\n_i = \min(n_i, \text{lf\\cap})\\. Applied after
  proportions are computed so that observed compositions are unaffected.
  Default `NULL` (no cap). This is a Multifan-CL legacy feature.

## Value

The input `data` list with the following elements appended or updated:

- `len_lower`, `len_upper`, `len_mid`:

  Length bin boundary and midpoint vectors derived from the scalar
  inputs.

- `lf_switch`:

  Passed through from the argument.

- `lf_obs_in`:

  Matrix of observed proportions, `[n_obs x n_len]`.

- `lf_n`:

  Numeric vector of sample sizes per observation row.

- `lf_fishery`:

  Integer vector of fishery index per observation row.

- `lf_fishery_f`:

  Integer vector of unique fishery indices with LF data.

- `lf_n_f`:

  Integer vector of observation counts per fishery.

- `lf_year`:

  Integer vector of model timestep index (1-based) per observation row.

- `lf_season`:

  Integer vector of season index (all 1).

- `lf_minbin`, `lf_maxbin`:

  Passed through from arguments.

- `lf_var_adjust`:

  Passed through from argument; numeric vector `[n_fishery]` of
  variance-adjustment divisors applied to `lf_n`.

- `removal_switch_f`:

  Integer vector `[n_fishery]` of removal flags (all 0).

- `n_lf`:

  Total number of observation rows.

- `lf_obs_data`:

  Named list `list(obs = lf_obs_list)` for plain-R access in
  `get_length_like`.

- `lf_obs_flat`:

  Flattened numeric vector of counts (for multinomial, `lf_switch = 1`).

- `lf_obs_ints`:

  Flattened integer vector of rounded counts (for Dirichlet-multinomial,
  `lf_switch = 3`).

- `lf_obs_prop`:

  Flattened numeric vector of normalised proportions (for Dirichlet,
  `lf_switch = 2`).

- `lf_nbins`:

  Number of bins used in each observation (scalar, derived from the
  first fishery's min/max bin setting).

## Details

Rows in `lf_wide` with a total sample size of zero are silently removed
before any other processing. Bin alignment between the data frame
columns and the model's length structure is checked with
[`stopifnot()`](https://rdrr.io/r/base/stopifnot.html).

## See also

[`get_length_like`](https://n-ducharmebarth-noaa.github.io/opal/reference/get_length_like.md),
[`opal_model`](https://n-ducharmebarth-noaa.github.io/opal/reference/opal_model.md)
