# Create a portable fitted opal model object

Captures the plain R state needed to reproduce a fitted opal model
without serializing the transient RTMB objective. The objective is
retained only in a session cache and can be rebuilt with
[`rebuild_opal_object()`](https://n-ducharmebarth-noaa.github.io/opal/reference/rebuild_opal_object.md).
Optional SparseNUTS output is normalized to a package-owned `opal_mcmc`
payload.

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
  optimizer = "nlminb"
)
```

## Arguments

- data:

  Named model data list used to construct `obj`.

- obj:

  Fitted RTMB objective created with `opal_model`.

- opt:

  Optimizer result, normally returned by
  [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html).

- bounds:

  Optional bounds data frame from
  [`get_bounds()`](https://n-ducharmebarth-noaa.github.io/opal/reference/get_bounds.md)
  or a list with `lower` and `upper`.

- control:

  Optional optimizer control list.

- estimability:

  Optional output from
  [`check_estimability()`](https://n-ducharmebarth-noaa.github.io/opal/reference/check_estimability.md).
  A compact summary is retained.

- diagnostics:

  Optional named list of fit diagnostics.

- metadata:

  Optional named list of user metadata.

- mcmc:

  Optional SparseNUTS-style fit, posterior matrix, or
  iteration-chain-variable array.

- mcmc_settings:

  Optional named list of sampler settings not already present in `mcmc`.

- derived:

  Optional named list of portable derived results, such as projections
  or retrospective summaries.

- makeadfun_args:

  Optional named list of additional arguments needed to rebuild
  `RTMB::MakeADFun()`. Core arguments are reserved.

- optimizer:

  Non-empty character name of the optimizer used.

## Value

An object inheriting from `opal_fit`.
