# Save and read portable opal fits

`save_opal_fit()` writes atomically and refuses to replace an existing
file unless requested. `read_opal_fit()` verifies the portable payload
and can rebuild its transient RTMB objective.

## Usage

``` r
save_opal_fit(x, file, compress = "gzip", overwrite = FALSE)

read_opal_fit(
  file,
  strict = FALSE,
  rebuild = strict,
  integrity = c("exact", "portable")
)
```

## Arguments

- x:

  An `opal_fit` object.

- file:

  Path to an RDS file.

- compress:

  Compression passed to
  [`saveRDS()`](https://rdrr.io/r/base/readRDS.html).

- overwrite:

  Replace an existing file.

- strict:

  Treat scientific-contract incompatibility as an error.

- rebuild:

  Rebuild and cache the RTMB objective after reading. Defaults to
  `strict`.

- integrity:

  Runtime-payload verification mode. `"exact"` requires the raw
  serialized payload identity to match. `"portable"` permits a mismatch
  across R versions, provided the rebuilt objective is verified.

## Value

`save_opal_fit()` invisibly returns the normalized path;
`read_opal_fit()` returns an `opal_fit`.
