# Check fitted or sampled Opal results

Checks retain the object even when diagnostics fail. Validation records
are tied to the exact fitted or sampled state and are separate from
lifecycle stage. Short MCMC smoke tests are not evidence of convergence.

## Usage

``` r
opal_check(
  x,
  scope = c("fit", "mcmc"),
  gradient_tolerance = 0.001,
  max_rhat = 1.01,
  min_ess = 100,
  stop_on_failure = FALSE
)
```

## Arguments

- x:

  An `opal_obj`.

- scope:

  Check the fitted optimum or MCMC.

- gradient_tolerance:

  Maximum absolute fitting gradient.

- max_rhat:

  Maximum rank-normalised R-hat.

- min_ess:

  Minimum bulk and tail effective sample sizes.

- stop_on_failure:

  Stop rather than warn when a check fails.

## Value

The object with an attached validation record.

## Details

Fit checks cover optimiser convergence, maximum absolute gradient, a
positive-definite Hessian, parameter bounds, and finite non-negative
numbers at age. MCMC checks require at least two chains, finite R-hat
and effective sample sizes within the thresholds, known sampler
diagnostics, no divergences, and no maximum-tree-depth hits. Missing
sampler diagnostics prevent a passing MCMC check, even if imported
parameter draws are usable. The record is stored in
`x$validation[[scope]]`, including metrics, settings, and the identity
of the checked result. These checks do not establish scientific adequacy
of a model or projection scenario.

## See also

Other assessment workflow:
[`opal_attach_fit()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_attach_fit.md),
[`opal_attach_mcmc()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_attach_mcmc.md),
[`opal_build()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_build.md),
[`opal_fit()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_fit.md),
[`opal_from_fit()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_from_fit.md),
[`opal_io`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_io.md),
[`opal_mcmc()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_mcmc.md),
[`opal_obj()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_obj.md),
[`opal_project()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_project.md),
[`opal_report()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_report.md),
[`opal_update()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_update.md)

## Examples

``` r
assessment <- opal_read(system.file("extdata", "opaka_quickstart_fit.rds",
                                   package = "opal"))
#> Warning: The saved fit runtime identity differs because it was serialized under a different R version; verifying by rebuild.
#> Warning: The saved fit runtime identity differs because it was serialized under a different R version; verifying by rebuild.
assessment <- opal_check(assessment, scope = "fit")
#> All parameters are estimable
#> Warning: The fit check failed; results are retained.
assessment$validation$fit$metrics
#> $convergence
#> [1] 0
#> 
#> $max_gradient
#> [1] 0.003773483
#> 
#> $positive_hessian
#> [1] TRUE
#> 
#> $inside_bounds
#> [1] TRUE
#> 
#> $valid_population
#> [1] TRUE
#> 
```
