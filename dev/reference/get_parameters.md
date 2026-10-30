# Get bundled initial parameter values

Loads the packaged initial parameter list that matches one of the
bundled opal model data sets.

## Usage

``` r
get_parameters(data = NULL, model = NULL)
```

## Arguments

- data:

  Optional model data list used to infer the bundled parameter set.

- model:

  Optional character model identifier. Supported values are
  `"opal_baseline"`, `"opakapaka"`, and `"wcpo_bet"`. Aliases
  `"baseline"`, `"opaka"`, and `"bet"` are also accepted. Supplying
  `model` is preferred when `data` has been modified after loading.

## Value

A `list` of initial parameter values.

## Details

If `data$k_mult_ages` is supplied, adds `log_k_mult` initialized to zero
(unit K multipliers).

Parameters for the alternative growth curves (see
[`get_growth`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/get_growth.md))
are added only when `data$growth_model` selects them:

- `growth_model = "richards"` adds `richards_b = 1`, which gives the VB
  curve.

- `growth_model = "cessation"` adds `log_L0` (initialized at `log_L1`),
  `log_rmax` (initialized at `log((L2 - L1) / (A2 - A1))`, length per
  model age step), `log_k_gc = 0` and `t50_gc` (initialized at the
  midpoint of the modeled ages, `(1 + n_age) / 2 + growth_age_offset`).

If `data$growth_sd_model = "sd"` (see
[`get_sd_at_age`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/get_sd_at_age.md)),
adds `log_SD1` and `log_SD2`, the log SDs of length-at-age at the
reference lengths, initialized at `log_CV1 + log_L1` and
`log_CV2 + log_L2` (the SDs implied by the CVs). `log_CV1` and `log_CV2`
are then unused.
