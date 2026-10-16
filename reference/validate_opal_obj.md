# Validate a staged Opal object

Checks structural and model-state integrity without building an
objective. Full result verification additionally checks posterior
payload identities.

## Usage

``` r
validate_opal_obj(x, results = FALSE)
```

## Arguments

- x:

  An `opal_obj`.

- results:

  Verify stored posterior draws and attempt history.

## Value

`x`, invisibly.
