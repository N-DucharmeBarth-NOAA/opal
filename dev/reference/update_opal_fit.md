# Update portable results attached to an opal fit

Adds or replaces normalized MCMC output and merges portable diagnostics,
derived results, or metadata. The fitted model state is unchanged.

## Usage

``` r
update_opal_fit(
  x,
  mcmc,
  mcmc_settings = list(),
  diagnostics = list(),
  derived = list(),
  metadata = list()
)
```

## Arguments

- x:

  An `opal_fit` object.

- mcmc:

  Optional replacement MCMC output. If omitted, existing output is
  retained; explicitly supply `NULL` to remove it.

- mcmc_settings:

  Optional settings merged into a replacement MCMC payload.

- diagnostics:

  Optional named list merged into fit diagnostics.

- derived:

  Optional named list merged into derived results.

- metadata:

  Optional named list merged into user metadata.

## Value

An updated `opal_fit`.
