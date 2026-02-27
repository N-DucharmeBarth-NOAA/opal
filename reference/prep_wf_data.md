# Prepare weight composition data for model input

Transforms a pivoted wide-format weight-frequency data frame into
model-ready arrays and vectors that get_weight_like() expects. Also
precomputes the rebinning matrix for converting predicted length
compositions to weight compositions.

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
  wf_cap = NULL
)
```

## Arguments

- data:

  list containing at minimum: lw_a, lw_b, len_bin_start, len_bin_width,
  n_len, n_fishery, and the weight bin scalars wt_bin_start,
  wt_bin_width, n_wt.

- wf_wide:

  data.frame in wide format: columns fishery, year, month, ts plus one
  numeric column per weight bin (named by bin value).

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

  Matrix of observed proportions (n_wf x n_wt).

- `wf_obs_flat`:

  Flattened numeric vector of counts (for multinomial, `wf_switch = 1`).

- `wf_obs_ints`:

  Flattened integer vector of rounded counts (for Dirichlet-multinomial,
  `wf_switch = 3`).

- `wf_obs_prop`:

  Flattened numeric vector of normalised proportions (for Dirichlet,
  `wf_switch = 2`).

- `wf_n`:

  Numeric vector of sample sizes per observation row.

- `wf_fishery`:

  Integer vector of fishery index per observation row.

- `wf_fishery_f`:

  Integer vector of unique fishery indices with WF data.

- `wf_n_f`:

  Integer vector of observation counts per fishery.

- `wf_year`:

  Integer vector of model timestep per observation row.

- `wf_minbin`, `wf_maxbin`:

  Passed through from arguments.

- `wf_var_adjust`:

  Passed through from argument; numeric vector `[n_fishery]` of
  variance-adjustment divisors applied to `wf_n`.

## See also

[`prep_lf_data`](https://n-ducharmebarth-noaa.github.io/opal/reference/prep_lf_data.md),
[`opal_model`](https://n-ducharmebarth-noaa.github.io/opal/reference/opal_model.md)
