# Save and read staged Opal objects

Every stage is portable. Saved objects contain a checksum of the
complete payload; their transient RTMB objective is never serialised.
Reads verify the checksum and scientific contract. Legacy fits are
converted with objective verification. No read operation optimises or
samples.

## Usage

``` r
opal_save(x, file, compress = "gzip", overwrite = FALSE)

opal_read(
  file,
  rebuild = FALSE,
  strict = TRUE,
  integrity = c("portable", "exact")
)
```

## Arguments

- x:

  An `opal_obj` or legacy `opal_fit` to save.

- file:

  RDS path.

- compress:

  RDS compression.

- overwrite:

  Allow replacement of an existing file.

- rebuild:

  Rebuild a configured object's runtime after reading.

- strict:

  Reject incompatible scientific contracts. FALSE permits inspection
  only; runtime access still rejects an incompatible model.

- integrity:

  Verification mode for legacy files: `"portable"` checks the rebuilt
  objective across R versions; `"exact"` also requires the original
  runtime checksum. New-format files always verify their payload
  checksum, independently of this legacy option.

## Value

Save returns the path invisibly; read returns an `opal_obj`.

## Details

Assign the return value of `opal_read()` to resume work in another
session. Use `rebuild = TRUE` to verify runtime reconstruction
immediately. Existing files are protected by default; set
`overwrite = TRUE` deliberately when saving an updated assessment to the
same path.

## See also

Other assessment workflow:
[`opal_attach_fit()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_attach_fit.md),
[`opal_attach_mcmc()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_attach_mcmc.md),
[`opal_build()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_build.md),
[`opal_check()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_check.md),
[`opal_fit()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_fit.md),
[`opal_from_fit()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_from_fit.md),
[`opal_mcmc()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_mcmc.md),
[`opal_obj()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_obj.md),
[`opal_project()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_project.md),
[`opal_report()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_report.md),
[`opal_update()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_update.md)

## Examples

``` r
inputs <- opaka_quickstart_inputs()
assessment <- opal_obj(inputs$data, inputs$parameters, inputs$map)
path <- tempfile(fileext = ".rds")
opal_save(assessment, path)
restored <- opal_read(path, rebuild = TRUE)
summary(restored)
#> <opal_obj> built
#> Active parameters: 127 
#> Fit check: not run  | MCMC check: not run 
unlink(path)
```
