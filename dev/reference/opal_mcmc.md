# Run MCMC for an Opal object

Sampling uses a fresh RTMB objective. It never changes the stored
optimum. `init = "auto"` uses the fitted point when available and the
configured active parameters otherwise. Failed sampler attempts are
recorded and warned about, preserving any existing posterior. The
default sampler skips internal optimisation; an MLE is optional. The
default stan metric supports bounds. Random effects are sampled jointly
by default, over their full domain, while fixed effects retain the
object's bounds. Change bounds with
[`opal_update()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_update.md).
Use `laplace = TRUE` and `metric = "unit"` for marginal sampling.

## Usage

``` r
opal_mcmc(
  x,
  sampler = "snuts",
  init = "auto",
  check = TRUE,
  check_args = list(),
  ...
)
```

## Arguments

- x:

  A configured or fitted `opal_obj`.

- sampler:

  `"snuts"` or a function accepting an `obj` argument and sampler
  settings.

- init:

  Initial-value policy or explicit sampler initial values.

- check:

  Run posterior diagnostics.

- check_args:

  Arguments passed to
  [`opal_check()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_check.md).

- ...:

  Named sampler settings, including seed, chains, cores, num_samples,
  and num_warmup. The obj, globals, lower, and upper arguments are
  reserved.

## Value

An updated `opal_obj`, including portable samples and attempt history.

## Details

Defaults are four chains, four cores, 1,000 warm-up iterations, and 500
retained iterations per chain. Set `seed` and choose effort appropriate
to the assessment. After sampling, inspect `x$validation$mcmc` and
`x$mcmc_history`; a failed or unchecked attempt cannot replace a
previously checked, passing posterior.
[`opal_as_tmbfit()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_as_tmbfit.md)
exposes the selected samples to SparseNUTS plotting and diagnostic
tools.

## See also

Other assessment workflow:
[`opal_attach_fit()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_attach_fit.md),
[`opal_attach_mcmc()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_attach_mcmc.md),
[`opal_build()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_build.md),
[`opal_check()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_check.md),
[`opal_derived()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_derived.md),
[`opal_example_inputs()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_example_inputs.md),
[`opal_fit()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_fit.md),
[`opal_from_fit()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_from_fit.md),
[`opal_io`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_io.md),
[`opal_obj()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_obj.md),
[`opal_posterior()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_posterior.md),
[`opal_project()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_project.md),
[`opal_report()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_report.md),
[`opal_update()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_update.md)

## Examples

``` r
if (FALSE) { # \dontrun{
assessment <- opal_read("assessment.rds")
assessment <- opal_mcmc(assessment, seed = 42, chains = 4, cores = 4,
                        num_warmup = 1000, num_samples = 500)
summary(assessment)
opal_save(assessment, "assessment-sampled.rds")
} # }
```
