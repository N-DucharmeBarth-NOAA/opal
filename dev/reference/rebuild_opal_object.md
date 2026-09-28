# Rebuild the RTMB objective for a portable opal fit

Recreates the objective from stored data, fitted parameters, map, and
random effects, then verifies the active parameter layout and saved
objective.

## Usage

``` r
rebuild_opal_object(
  x,
  strict = TRUE,
  check_objective = TRUE,
  tolerance = 1e-06,
  silent = FALSE,
  cache = TRUE,
  integrity = c("exact", "portable")
)
```

## Arguments

- x:

  An `opal_fit` object.

- strict:

  Stop on compatibility or parameter-order differences.

- check_objective:

  Compare the rebuilt and saved objective values.

- tolerance:

  Relative tolerance passed to
  [`all.equal()`](https://rdrr.io/r/base/all.equal.html).

- silent:

  Passed to `RTMB::MakeADFun()`.

- cache:

  Cache the rebuilt objective for this R session.

- integrity:

  Runtime-payload verification mode. `"exact"` requires the raw
  serialized payload identity to match. `"portable"` permits a mismatch
  across R versions, provided the rebuilt objective is verified.

## Value

A newly constructed RTMB objective.
