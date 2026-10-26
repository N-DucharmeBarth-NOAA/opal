# Indian Ocean yellowfin tuna (simulated) starting parameters

Starting values for
[io_yft_data](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/io_yft_data.md),
at the SS3 solution: SS3's unfished spawning biomass (`log_B0`),
recruitment deviations (`rdev_y`), steepness, sigma_R and growth
(including the age-specific K multipliers `log_k_mult`). Selectivity
(`par_sel`) starts from opal curves fitted by least squares to SS3's
selectivity-at-age projected onto length.

## Usage

``` r
io_yft_parameters
```

## Format

A named list of 18 parameter vectors and matrices, including `log_B0`,
`log_cpue_q` (4), `rdev_y` (256), `par_sel` (20 x 6) and `log_k_mult`
(9).

## Source

As for
[io_yft_data](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/io_yft_data.md).

## See also

[io_yft_map](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/io_yft_map.md)
