# `opal` model regression baseline

A saved reference run for the WCPO bigeye tuna assessment that documents
the objective, gradient, and RTMB report output used by the
`vignettes/baseline.Rmd` workflow to detect numerical regressions and
timing changes before and after refactors.

## Usage

``` r
opal_baseline
```

## Format

A named list accessible via `data(opal_baseline)` and packaged as
`data/opal_baseline.rda`.

## Source

`vignettes/baseline.Rmd` (the "Save or compare" chunk) which reproduces
the reference run saved to `data/opal_baseline.rda`.

## Details

### Baseline contents

- timestamp:

  Time when the baseline was created.

- description:

  Free-form label ("Pre-refactoring baseline" by default).

- opal_version:

  Version string of the installed `opal` package.

- dimensions:

  Named list with `n_year`, `n_age`, `n_fishery`, `n_season`, `n_len`,
  `n_lf`, `n_wt`, and `n_wf` from the model inputs.

- nll:

  Negative log-likelihood at the reference parameter set.

- gradient:

  Gradient vector evaluated at `obj$par` for the baseline.

- max_gr:

  Maximum absolute gradient entry.

- timing:

  List with `fn_elapsed`, `gr_elapsed`, `fn_median`, and `gr_median`
  values that summarize runtime measurements used in the comparison
  vignette.

- report:

  Sublist containing the baseline's `number_ysa`, `spawning_biomass_y`,
  catch and harvest rate predictions (`catch_pred_fya`,
  `catch_pred_ysf`, `hrate_ysa`, `hrate_ysfa`), selectivity curves
  (`sel_fya`), and log-likelihood components (`lp_prior`, `lp_penalty`,
  `lp_rec`, `lp_cpue`, `lp_lf`, `lp_wf`).

- par_values:

  Vector of `obj$par` values that produced the baseline.

## Examples

``` r
if (FALSE) { # \dontrun{
  data(opal_baseline)
  str(opal_baseline$dimensions)
  opal_baseline$timing$fn_median
} # }
```
