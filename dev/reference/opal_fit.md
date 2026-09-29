# Fit an Opal assessment

Pass an
[`opal_obj()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_obj.md)
to optimise its model with
[`nlminb()`](https://rdrr.io/r/stats/nlminb.html), retain the optimiser
passes, and run
[`opal_check()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_check.md).
The default is two sequential passes. A completed optimisation does not
guarantee convergence: inspect `summary(assessment)` and
`assessment$validation$fit` before using results.

## Usage

``` r
opal_fit(
  data,
  obj,
  opt,
  bounds = NULL,
  control = NULL,
  estimability = NULL,
  diagnostics = list(),
  metadata = list(),
  mcmc = NULL,
  mcmc_settings = list(),
  derived = list(),
  makeadfun_args = list(),
  optimizer = "nlminb",
  n_passes = 2L,
  check = TRUE,
  check_args = list()
)
```

## Arguments

- data:

  An `opal_obj` to optimise, or the legacy named model data list used to
  construct `obj`. With an `opal_obj`, returns an updated object.

- obj:

  Legacy only: fitted RTMB objective. Omit for an `opal_obj`.

- opt:

  Legacy only: optimiser result. Omit for an `opal_obj`.

- bounds:

  Legacy only: bounds data frame from
  [`get_bounds()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/get_bounds.md)
  or a list with `lower` and `upper`.

- control:

  Optional control list passed to
  [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html).

- estimability:

  Legacy only: output from
  [`check_estimability()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/check_estimability.md).
  A compact summary is retained.

- diagnostics:

  Legacy only: named list of fit diagnostics.

- metadata:

  Legacy only: named list of user metadata.

- mcmc:

  Legacy only: SparseNUTS-style fit, posterior matrix, or
  iteration-chain-variable array.

- mcmc_settings:

  Legacy only: named list of sampler settings not already present in
  `mcmc`.

- derived:

  Legacy only: named list of portable derived results, such as
  projections or retrospective summaries.

- makeadfun_args:

  Legacy only: named list of additional arguments needed to rebuild
  `RTMB::MakeADFun()`. Core arguments are reserved.

- optimizer:

  Legacy only: non-empty character name of the optimiser used.

- n_passes:

  Number of sequential optimisation passes for an `opal_obj`.

- check:

  Run fitting diagnostics for an `opal_obj`.

- check_args:

  Named arguments passed to
  [`opal_check()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_check.md).

## Value

With an `opal_obj`, an updated `opal_obj`. Legacy constructor calls
return `opal_fit`; use
[`opal_from_fit()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_from_fit.md)
to migrate them.

## Details

For the object workflow, use
`opal_fit(assessment, n_passes = 2L, control = NULL, check = TRUE, check_args = list())`.
Configuration comes from the object; change bounds, maps, or metadata
with
[`opal_update()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_update.md).
Omitting `control` uses the stored controls, or the Opal defaults.
Refitting an unchanged target preserves its posterior with the
provenance of the fit that originally supplied it. Failed checks warn
and retain results.

The legacy call `opal_fit(data, obj, opt, ...)` remains available to
capture an already fitted RTMB model. Arguments marked legacy below
apply only to that interface. New custom-optimiser workflows should use
[`opal_attach_fit()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_attach_fit.md).
Convert an existing legacy fit with
[`opal_from_fit()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_from_fit.md).

## See also

Other assessment workflow:
[`opal_attach_fit()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_attach_fit.md),
[`opal_attach_mcmc()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_attach_mcmc.md),
[`opal_build()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_build.md),
[`opal_check()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_check.md),
[`opal_derived()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_derived.md),
[`opal_example_inputs()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_example_inputs.md),
[`opal_from_fit()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_from_fit.md),
[`opal_io`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_io.md),
[`opal_mcmc()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_mcmc.md),
[`opal_obj()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_obj.md),
[`opal_posterior()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_posterior.md),
[`opal_project()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_project.md),
[`opal_report()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_report.md),
[`opal_update()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_update.md)

## Examples

``` r
inputs <- opaka_quickstart_inputs()
assessment <- opal_obj(inputs$data, inputs$parameters, inputs$map)
# \donttest{
assessment <- opal_fit(assessment)
#> Warning: NA/NaN function evaluation
#> All parameters are estimable
#> Warning: The fit check failed; results are retained.
summary(assessment)
#> <opal_obj> fitted
#> Active parameters: 127 
#> Objective: 1832.412 
#> Fit check: failed  | MCMC check: not run 
assessment$validation$fit
#> $passes
#> [1] FALSE
#> 
#> $metrics
#> $metrics$convergence
#> [1] 0
#> 
#> $metrics$max_gradient
#> [1] 0.001041864
#> 
#> $metrics$positive_hessian
#> [1] TRUE
#> 
#> $metrics$inside_bounds
#> [1] TRUE
#> 
#> $metrics$valid_population
#> [1] TRUE
#> 
#> $metrics$biology
#> $metrics$biology$passes
#> [1] TRUE
#> 
#> $metrics$biology$checks
#>  initial_equilibrium      harvest_penalty           population 
#>                 TRUE                 TRUE                 TRUE 
#>            mortality             maturity   spawning_potential 
#>                 TRUE                 TRUE                 TRUE 
#>               weight          selectivity            steepness 
#>                 TRUE                 TRUE                 TRUE 
#>          recruitment      spawning_output              harvest 
#>                 TRUE                 TRUE                 TRUE 
#>            depletion catch_reconstruction 
#>                 TRUE                 TRUE 
#> 
#> $metrics$biology$metrics
#> $metrics$biology$metrics$initial_penalty
#> [1] 0
#> 
#> $metrics$biology$metrics$total_penalty
#> [1] 0
#> 
#> $metrics$biology$metrics$max_harvest
#> [1] 0.09746108
#> 
#> $metrics$biology$metrics$max_relative_catch_error
#> [1] 6.704801e-10
#> 
#> 
#> 
#> 
#> $error
#> NULL
#> 
#> $identity
#> $identity$target
#> [1] "ced004c827f173c7b89c6399529db21c"
#> 
#> $identity$fit
#> [1] "7459c22d1102869eeae932f6f18cec4b"
#> 
#> 
#> $version
#> [1] "opal_validation_v2"
#> 
#> $settings
#> $settings$gradient_tolerance
#> [1] 0.001
#> 
#> $settings$max_rhat
#> [1] 1.01
#> 
#> $settings$min_ess
#> [1] 100
#> 
#> $settings$penalty_tolerance
#> [1] 1e-10
#> 
#> $settings$catch_tolerance
#> [1] 1e-06
#> 
#> 
#> $checked_at
#> [1] "2026-09-29 03:43:21 UTC"
#> 
#> $payload_id
#> [1] "6db32fb2dc63c56ac9bb06f4c952c3ac"
#> 
# }
```
