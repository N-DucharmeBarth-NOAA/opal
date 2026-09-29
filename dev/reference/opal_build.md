# Build or access an Opal runtime objective

`opal_build()` resolves omitted parameters, maps, and bounds, then
checks that the initial objective and gradient are finite. It does not
optimise.
[`opal_fit()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_fit.md)
and
[`opal_mcmc()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_mcmc.md)
call it automatically.

## Usage

``` r
opal_build(x, silent = TRUE)

opal_rtmb(x, fresh = FALSE)
```

## Arguments

- x:

  An `opal_obj`.

- silent:

  Silence RTMB construction messages.

- fresh:

  Build an isolated objective instead of accessing the cache.

## Value

`opal_build()` returns the updated object; `opal_rtmb()` returns a
transient RTMB objective. Use `fresh = TRUE` for isolated mutable work.

## Details

Call `opal_build()` before requesting the runtime of an unresolved
object. Use `opal_rtmb(assessment, fresh = TRUE)` for external
optimisation, simulation, profiling, or other work that mutates RTMB
state. The runtime is transient; save the assessment with
[`opal_save()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_io.md),
not the RTMB object.

## See also

Other assessment workflow:
[`opal_attach_fit()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_attach_fit.md),
[`opal_attach_mcmc()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_attach_mcmc.md),
[`opal_check()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_check.md),
[`opal_derived()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_derived.md),
[`opal_example_inputs()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_example_inputs.md),
[`opal_fit()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_fit.md),
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
assessment <- opal_build(opal_obj(inputs$data, inputs$parameters, inputs$map))
runtime <- opal_rtmb(assessment, fresh = TRUE)
runtime$fn(runtime$par)
#> [1] 67665.47
```
