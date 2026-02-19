# Calculate recruitment

Computes recruitment based on Beverton-Holt with depensation and
log-normal deviations.

## Usage

``` r
get_recruitment(sbio, rdev, B0, alpha, beta, sigma_r = 0.6)
```

## Arguments

- sbio:

  Spawning biomass.

- rdev:

  Recruitment deviations.

- B0:

  Unfished biomass.

- alpha, beta:

  Beverton-Holt stock recruitment parameters.

- sigma_r:

  Lognormal SD of recruitment deviations.

## Value

Recruitment value (numeric).
