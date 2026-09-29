# Calculate one-step-ahead observation residuals

Adds CPUE, length-composition, and weight-composition OSA diagnostics to
an Opal object. Use
[`plot_osa_sdnr()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/plot_osa_sdnr.md)
for the dataset-level SDNR comparison.

## Usage

``` r
opal_osa(x, seed = 123L, conf = 0.95, name = "osa")
```

## Arguments

- x:

  A fitted `opal_obj`.

- seed:

  Random seed for discrete residuals; restored on exit. `NULL` uses and
  advances the current random-number stream.

- conf:

  Confidence level for approximate chi-squared SDNR intervals.

- name:

  Name under `x$derived` for the diagnostics.

## Value

An updated object.
[`opal_derived()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_derived.md)
returns its residual table, dataset summaries, and calculation settings.

## Details

For fixed-effect fits, CPUE uses its lognormal CDF, and composition
residuals use exact sequential binomial, beta, or beta-binomial CDFs for
multinomial, Dirichlet, or Dirichlet-multinomial observations.
Parameters are held at their fitted values. For random-effect models,
RTMB's `oneStepPredict()` integrates latent states using its
Laplace-based OSA approximation. Other observation streams remain
conditioned on. There is no silent substitution of conditional residuals
for marginal residuals.

Observation order is the stored input order, then ascending composition
bin. The final composition bin is constrained and omitted. Removed
fisheries, zero sample sizes, exhausted counts, and numerical failures
remain in the table with explicit exclusion reasons. Multinomial counts
are rounded as in the fitted likelihood; no new pseudocounts are added.
SDNR should be near one, but its interval is approximate because fitted
parameters and finite samples affect calibration. Catch is conditioned
on and has no observation likelihood; process priors are not
observations.

## See also

Other assessment diagnostics:
[`opal_diagnose()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_diagnose.md),
[`plot_composition()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/plot_composition.md),
[`plot_osa_residuals()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/plot_osa_residuals.md),
[`plot_osa_sdnr()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/plot_osa_sdnr.md),
[`plot_prior_posterior()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/plot_prior_posterior.md)
