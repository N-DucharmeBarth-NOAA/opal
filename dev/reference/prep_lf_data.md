# Prepare length composition data for model input

Transforms wide or long length-frequency observations (one composition
per fishery x timestep and optional source) into model-ready arrays and
vectors that `opal_model` and `get_length_like` expect. The function
validates bin alignment against the length structure in `data`,
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
  lf_cap = NULL,
  lf_addtocomp = 1e-08,
  lf_upper_edge = NULL,
  lf_tails = c("accumulate", "truncate"),
  lf_edge_tol = 1e-08 * data$len_bin_width,
  lf_allow_misaligned = FALSE
)
```

## Arguments

- data:

  a `list` containing at minimum the length structure scalars
  `len_bin_start`, `len_bin_width`, `n_len`, and `n_fishery`.

- lf_wide:

  A wide data frame, list of wide data frames, or long data frame. Wide
  frames contain `fishery`, `year`, `month`, `ts`, and numeric count
  columns named by bin lower edges. Long input contains the same
  metadata plus `bin_lower`, `bin_upper`, and `count`. Long observations
  are grouped by `fishery`, `ts`, and optionally `source`; calendar
  metadata is taken from the first row.

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

- lf_addtocomp:

  small non-negative numeric constant added to observed composition
  proportions after tail compression and before renormalisation. This
  robustifies zero bins. Default `1e-08`.

- lf_upper_edge:

  Final upper edge, or one per element of a list of wide data frames.
  Defaults to the model length-grid upper edge,
  `data$len_bin_start + data$n_len * data$len_bin_width`. Ignored for
  long input, which supplies both boundaries explicitly.

- lf_tails:

  Prediction tail handling: `"accumulate"` (default) extends the outer
  matrix edges, or `"truncate"` excludes outside mass.

- lf_edge_tol:

  Absolute edge tolerance; defaults to `1e-8 * data$len_bin_width`.
  Nearby model edges snap silently.

- lf_allow_misaligned:

  Allow observed edges to split model bins, with a diagnostic warning.
  The default rejects all edges not matching a model boundary, including
  edges beyond the model range.

## Value

The input `data` list with the following elements appended or updated:

- `len_lower`, `len_upper`, `len_mid`:

  Length bin boundary and midpoint vectors derived from the scalar
  inputs.

- `lf_switch`:

  Passed through from the argument.

- `lf_obs_in`:

  Observed proportions in retained input order: a matrix with one row
  per observation when bin counts agree, otherwise a list of row
  vectors. Native input has `n_len` columns.

- `lf_n`:

  Numeric vector of sample sizes per observation row.

- `lf_n_int`:

  Integer Dirichlet-multinomial sample sizes obtained by half-up
  rounding of `lf_n`.

- `lf_fishery`:

  Integer vector of fishery index per observation row.

- `lf_fishery_f`:

  Integer vector of unique fishery indices with LF data.

- `lf_n_f`:

  Integer vector of observation counts per fishery.

- `lf_year`:

  Integer vector of model timestep index (1-based) per observation row.

- `lf_year_fi`, `lf_n_fi`, `lf_n_int_fi`, `lf_row_fi`:

  Lists split by fishery containing model timesteps, effective sample
  sizes, and row indices. Precomputed so
  [`opal_model()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_model.md)
  does not rebuild them on every objective evaluation.

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

  Flattened integer counts that sum to `lf_n_int` per observation (for
  Dirichlet-multinomial, `lf_switch = 3`).

- `lf_obs_prop`:

  Flattened numeric vector of normalised proportions (for Dirichlet,
  `lf_switch = 2`).

- `lf_addtocomp`:

  Stored value of the add-to-composition constant used to robustify
  observed composition bins.

- `lf_nbins`:

  Legacy scalar bin count for the first observation in fishery-grouped
  flat order, or zero for empty input. Use `lf_nbins_i` for
  heterogeneous observations.

## Details

Inputs may also be lists of wide data frames or long data frames.
Observed counts remain on their original grid. `lf_scheme_edges` is a
deduplicated list of observed edges, including native min/max pooling.
`lf_rebin_s` contains prediction matrices (observed bins by model length
bins). `lf_scheme` and `lf_nbins_i` follow fishery-grouped flat
observation order, while `lf_row_fi` indexes retained input order. The
logical `lf_native_i`, in flat observation order, identifies
observations originally on the full native grid before min/max pooling.
Heterogeneous observations use lists instead of rectangular matrices.

Splitting a model length bin at an observed boundary assumes uniform
prediction mass within that model bin. This approximation is not needed
when observed boundaries match model boundaries, so nonmatching
boundaries are rejected by default. Setting `lf_allow_misaligned = TRUE`
opts into the approximation with a diagnostic warning. Check bin
labeling and units before changing the grid. Diagnostics suggest a finer
width dividing the observed offsets; edges beyond the model range also
require extending that range. Changing the grid changes the mean and SD
of `len_mid`, which scale length-based selectivity parameters: existing
`par_sel` values are not directly transferable.

Observations with a total sample size of zero are silently removed after
validating counts and edges. Bin alignment against the model's length
structure is checked with the edge tolerance. Dirichlet-multinomial
counts use half-up rounding of effective sample sizes and
largest-remainder allocation. Observations that round to zero are
skipped by the Dirichlet-multinomial likelihood.

## See also

[`get_length_like`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/get_length_like.md),
[`opal_model`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_model.md)
