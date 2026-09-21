# Weight Composition Likelihood

Computes likelihood for observed weight compositions using the PLA
matrix and a precomputed rebinning matrix. The prediction pipeline is:
catch_at_age -\> PLA -\> pred_at_length -\> rebin_matrix -\>
pred_at_weight. Gradients propagate through PLA to growth parameters
when estimated.

## Usage

``` r
get_weight_like(
  wf_obs_flat,
  wf_obs_ints,
  wf_obs_prop,
  catch_pred_fya,
  pla,
  wf_rebin_matrix,
  wf_n_f,
  wf_fishery_f,
  wf_year_fi,
  wf_n_fi,
  wf_minbin,
  wf_maxbin,
  removal_switch_f,
  wf_switch,
  n_wt,
  n_wf,
  log_wf_tau,
  wf_addtocomp = 1e-08,
  wf_n_int_fi = NULL
)
```

## Arguments

- wf_obs_flat:

  numeric vector of unrounded weight comp counts (wf_switch=1).

- wf_obs_ints:

  integer vector of integer weight comp counts (wf_switch=3).

- wf_obs_prop:

  numeric vector of weight comp proportions (wf_switch=2).

- catch_pred_fya:

  3D array `[n_fishery, n_year, n_age]` of predicted catch-at-age from
  do_dynamics().

- pla:

  matrix `[n_len, n_age]` probability-of-length-at-age from get_pla().
  On the AD tape when growth parameters are estimated.

- wf_rebin_matrix:

  matrix `[n_wt, n_len]` precomputed rebinning weights from
  prep_wf_data().

- wf_n_f:

  integer vector `[n_fishery]` of observation counts per fishery.

- wf_fishery_f:

  integer vector of fishery indices with WF data.

- wf_year_fi:

  list of integer vectors of year indices per fishery.

- wf_n_fi:

  list of integer vectors of sample sizes per fishery.

- wf_minbin:

  integer vector `[n_fishery]` minimum weight bin index.

- wf_maxbin:

  integer vector `[n_fishery]` maximum weight bin index.

- removal_switch_f:

  integer vector `[n_fishery]` removal flags.

- wf_switch:

  integer likelihood type (1=multinomial, 2=Dirichlet, 3=DM).

- n_wt:

  integer number of weight bins.

- n_wf:

  integer total number of WF observations.

- log_wf_tau:

  numeric vector `[n_fishery]` log-scale variance adjustment.

- wf_addtocomp:

  small non-negative numeric constant added to predicted proportions
  before normalisation to robustify zero bins. Default `1e-08`.

- wf_n_int_fi:

  optional list of integer Dirichlet-multinomial sample sizes.

## Value

numeric vector of negative log-likelihood contributions, one per
observation.

## Details

For Dirichlet, concentration is `pred * n_i * exp(log_wf_tau[f])`; for
Dirichlet-multinomial, it is `pred * exp(log_wf_tau[f])`. The
multinomial branch uses RTMB's integer-rounded count behavior. The
Dirichlet-multinomial is intended for real counts with variance
adjustment equal to one; use Dirichlet for fractional effective sample
sizes.
