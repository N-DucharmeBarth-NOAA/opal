# Indian Ocean yellowfin tuna (simulated) parameter map

The map for fit "C2" of
[io_yft_data](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/io_yft_data.md):
`log_B0`, the 256 recruitment deviations, the four quarterly
catchabilities and 49 selectivity parameters are estimated (310 in
total). Growth, natural mortality, steepness and sigma_R are fixed at
SS3's values.

## Usage

``` r
io_yft_map
```

## Format

A named list of factors, as for the `map` argument of
[`opal_obj()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_obj.md).

## Source

As for
[io_yft_data](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/io_yft_data.md).

## Details

Selectivity estimates the opal equivalents of the parameters SS3
estimated: columns 1–4 of the double normal for fisheries 1–14 and both
logistic parameters for the index; fisheries 15–16 mirror 14, and the
index series 18–20 mirror 17. Nine parameters the length data do not
inform (the plateau and descending width of fisheries 4, 5, 7 and 13,
and the descending width of fishery 9) are fixed at their starting
values.
