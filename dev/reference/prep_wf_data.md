# Prepare weight composition data for model input

Transforms wide or long weight-frequency observations into model-ready
arrays and vectors that get_weight_like() expects. Also precomputes the
rebinning matrix for converting predicted length compositions to weight
compositions.

## Usage

``` r
prep_wf_data(
  data,
  wf_wide,
  wf_keep_fisheries = NULL,
  wf_switch = 1L,
  wf_minbin = NULL,
  wf_maxbin = NULL,
  wf_var_adjust = NULL,
  wf_cap = NULL,
  wf_addtocomp = 1e-08,
  wf_upper_edge = NULL,
  wf_tails = c("accumulate", "truncate"),
  wf_edge_tol = 1e-08 * data$wt_bin_width
)
```

## Arguments

- data:

  list containing at minimum: lw_a, lw_b, len_bin_start, len_bin_width,
  n_len, n_fishery, and the weight bin scalars wt_bin_start,
  wt_bin_width, n_wt.

- wf_wide:

  A wide data frame, list of wide data frames, or long data frame. Wide
  frames contain `fishery`, `year`, `month`, `ts`, and numeric count
  columns named by bin lower edges. Long input contains the same
  metadata plus `bin_lower`, `bin_upper`, and `count`. Long observations
  are grouped by `fishery`, `ts`, and optionally `source`; calendar
  metadata is taken from the first row.

- wf_keep_fisheries:

  integer vector of fishery indices to retain (NULL = keep all).

- wf_switch:

  integer likelihood type (1=multinomial, 2=Dirichlet,
  3=Dirichlet-multinomial). Default 1L.

- wf_minbin:

  integer vector (length n_fishery) min weight bin index. Defaults to
  rep(1L, n_fishery).

- wf_maxbin:

  integer vector (length n_fishery) max weight bin index. Defaults to
  rep(n_wt, n_fishery).

- wf_var_adjust:

  numeric vector (length n_fishery) variance adjustment divisors.
  Defaults to rep(1, n_fishery).

- wf_cap:

  positive integer (or `NULL`). When supplied, each effective sample
  size is capped: \\n_i = \min(n_i, \text{wf\\cap})\\. Applied after
  proportions are computed so that observed compositions are unaffected.
  Default `NULL` (no cap). This is a Multifan-CL legacy feature.

- wf_addtocomp:

  small non-negative numeric constant added to observed composition
  proportions after tail compression and before renormalisation. This
  robustifies zero bins. Default `1e-08`.

- wf_upper_edge:

  Final upper edge, or one per element of a list of wide data frames.
  Defaults to the largest observed bin lower edge plus
  `data$wt_bin_width`. Ignored for long input, which supplies both
  boundaries explicitly.

- wf_tails:

  Prediction tail handling: `"accumulate"` (default) extends custom
  scheme matrix edges to include the converted length range;
  `"truncate"` excludes outside prediction mass. Native schemes retain
  the legacy finite intermediate weight range; see Details.

- wf_edge_tol:

  Absolute tolerance for deduplicating scheme edges; defaults to
  `1e-8 * data$wt_bin_width`.

## Value

data list with the following weight composition elements appended:

- `wf_switch`:

  Passed through from the argument.

- `wt_lower`, `wt_upper`, `wt_mid`:

  Weight bin boundary and midpoint vectors derived from the scalar
  inputs.

- `wt_bin_edges`:

  Weight bin boundary vector (length n_wt + 1).

- `wf_rebin_matrix`:

  Precomputed rebinning matrix (n_wt x n_len) for converting predicted
  length compositions to weight compositions.

- `n_wf`:

  Total number of WF observation rows.

- `wf_obs_in`:

  Observed proportions in retained input order: a matrix with one row
  per observation when bin counts agree, otherwise a list of row
  vectors. Native input has `n_wt` columns.

- `wf_obs_flat`:

  Flattened numeric vector of counts (for multinomial, `wf_switch = 1`).

- `wf_obs_ints`:

  Flattened integer counts that sum to `wf_n_int` per observation (for
  Dirichlet-multinomial, `wf_switch = 3`).

- `wf_obs_prop`:

  Flattened numeric vector of normalised proportions (for Dirichlet,
  `wf_switch = 2`).

- `wf_addtocomp`:

  Stored value of the add-to-composition constant used to robustify
  observed composition bins.

- `wf_n`:

  Numeric vector of sample sizes per observation row.

- `wf_n_int`:

  Integer Dirichlet-multinomial sample sizes obtained by half-up
  rounding of `wf_n`.

- `wf_fishery`:

  Integer vector of fishery index per observation row.

- `wf_fishery_f`:

  Integer vector of unique fishery indices with WF data.

- `wf_n_f`:

  Integer vector of observation counts per fishery.

- `wf_year`:

  Integer vector of model timestep per observation row.

- `wf_year_fi`, `wf_n_fi`, `wf_n_int_fi`, `wf_row_fi`:

  Lists split by fishery containing model timesteps, effective sample
  sizes, and row indices. Precomputed so
  [`opal_model()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_model.md)
  does not rebuild them on every objective evaluation.

- `wf_minbin`, `wf_maxbin`:

  Passed through from arguments.

- `wf_var_adjust`:

  Passed through from argument; numeric vector `[n_fishery]` of
  variance-adjustment divisors applied to `wf_n`.

## Details

Input may be one wide data frame, a list of wide frames, or long data.
`wf_scheme_edges` records observed edges, including native min/max
pooling; `wf_rebin_s` contains observed-bin by model-length-bin
prediction matrices. Length edges are converted using `lw_a * L^lw_b`.
Native schemes preserve the legacy two-stage weight mapping and, by
default, tail pooling. Even under `"accumulate"`, native predictions
outside the full native weight grid remain excluded, preserving the old
numerical mapping rather than accumulating the entire converted length
range. With `wf_tails = "truncate"`, their prediction matrices select
native min/max rows without accumulating excluded prediction mass.
Prediction tail policy does not change observed native min/max count
pooling. Native and custom weight operators remain distinct even when
their displayed edges coincide, because their prediction tail handling
can differ. Within each mapping type, schemes deduplicate observed edges
within `wf_edge_tol`. Registry edges are never expanded for prediction
tails. `wf_scheme` and `wf_nbins_i` follow fishery-grouped flat order;
`wf_row_fi` indexes retained input order. Heterogeneous observations use
lists where a rectangular matrix is not possible. The logical
`wf_native_i`, in flat observation order, identifies observations
originally on the full native grid before min/max pooling. Observed bins
narrower than overlapping converted-length source bins trigger a
warning, as do weight schemes entirely outside the source range.

Dirichlet-multinomial counts use half-up rounding of effective sample
sizes and largest-remainder allocation. Observations that round to zero
are skipped by the Dirichlet-multinomial likelihood.

## See also

[`prep_lf_data`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/prep_lf_data.md),
[`opal_model`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_model.md)
