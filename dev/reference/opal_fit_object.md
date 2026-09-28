# Access the runtime objective for an opal fit

Access the runtime objective for an opal fit

## Usage

``` r
opal_fit_object(x, fresh = FALSE, integrity = c("exact", "portable"))
```

## Arguments

- x:

  An `opal_fit` object.

- fresh:

  Construct an isolated objective instead of using the session cache.

- integrity:

  Runtime-payload verification mode. `"exact"` requires the raw
  serialized payload identity to match. `"portable"` permits a mismatch
  across R versions, provided the rebuilt objective is verified.

## Value

A fitted RTMB objective.
