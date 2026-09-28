# Project an Opal object and retain the result

Runs
[`project_dynamics()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/project_dynamics.md)
and stores the result, settings, and the identity of its source fit or
posterior. Inputs for future recruitment, selectivity, and catch remain
explicit scientific choices. Updating the model or replacing the source
fit/posterior clears stored projections.

## Usage

``` r
opal_project(x, uncertainty = NULL, name = "projection", seed = NULL, ...)
```

## Arguments

- x:

  A fitted or sampled `opal_obj`.

- uncertainty:

  Either `"mvn"` or `"mcmc"`; required when a posterior exists.

- name:

  Name of the result in `x$derived`.

- seed:

  Optional random seed, restored after the projection.

- ...:

  Arguments passed to
  [`project_dynamics()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/project_dynamics.md),
  including future inputs.

## Value

An updated `opal_obj` with a stored projection.

## Details

Read the stored output from `x$derived[[name]]$result`; `settings`
records the future inputs and seed, and `identity` records the source
result. Set a seed separately before stochastic recruitment or
selectivity helpers; the `seed` argument here controls the dynamics
projection only. For a complete worked example, see
[`vignette("projections")`](https://n-ducharmebarth-noaa.github.io/opal/dev/articles/projections.md).

## See also

[`project_rec_devs()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/project_rec_devs.md),
[`project_selectivity()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/project_selectivity.md),
[`project_dynamics()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/project_dynamics.md)

Other assessment workflow:
[`opal_attach_fit()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_attach_fit.md),
[`opal_attach_mcmc()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_attach_mcmc.md),
[`opal_build()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_build.md),
[`opal_check()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_check.md),
[`opal_fit()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_fit.md),
[`opal_from_fit()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_from_fit.md),
[`opal_io`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_io.md),
[`opal_mcmc()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_mcmc.md),
[`opal_obj()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_obj.md),
[`opal_report()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_report.md),
[`opal_update()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_update.md)
