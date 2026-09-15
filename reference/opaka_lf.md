# 'Opakapaka Length Frequency Data

Length frequency composition data extracted from a Stock Synthesis 3
(SS3) model and prepared for RTMB (R Template Model Builder) stock
assessment modeling. Observations are organized by fishery, year, month,
and 5 cm length bins.

## Usage

``` r
opaka_lf
```

## Format

A data.table with 7 columns: year (integer), month (integer), ts
(integer), fishery (numeric), bin (numeric), value (numeric), week
(numeric). Multiple rows per fishery-year (one per 5 cm length bin).
Total rows = 1394 (82 observation-years x 17 bins).

## Source

Extracted from Stock Synthesis 3 data file (`data_echo.ss_new`).
Historical period only (1949-2023).

## Details

### Data workflow

**Preparation steps:**

1.  Parse length composition observations from SS3 data file
    (`data_echo.ss_new`)

2.  Extract proportions and sample sizes by fleet-year

3.  Convert proportions to pseudo-counts (proportion x Nsamp)

4.  Reshape to long format (one row per fishery-year-length bin)

5.  Truncate to historical period (1949-2023)

### Output columns

- year:

  Calendar year (1949-2023). Integer.

- month:

  Observation month (1 for all observations). Integer.

- ts:

  Model timestep number (1-75). Integer.

- fishery:

  Fleet index: 1 = Comm (75 obs), 3 = ResFish (7 obs). Numeric.

- bin:

  Length bin lower edge (cm). 5 cm bins: 5, 10, 15, ..., 85. Numeric.

- value:

  Pseudo-count (proportion x Nsamp) in this length bin. Numeric.

- week:

  Temporal indicator (always 1). Numeric.

### Fleet coverage

- Fleet 1 (Comm):

  75 observation-years (1949-2023). Sample sizes range from ~23 to ~755.

- Fleet 2 (Non_comm):

  No length composition data available.

- Fleet 3 (ResFish):

  7 observation-years (2017-2023). Fixed sample size of 60 per year.

### Notes

- Data length bins are 5 cm wide (5, 10, 15, ..., 85 cm), giving 17 bins

- Values are pseudo-counts derived from SS3 proportions x Nsamp

- Sum of values across bins per fishery-year gives the effective sample
  size

- Very small pseudo-counts (\< 1e-10) represent SS3's numerical zeros
  and should be treated as absent observations

## See also

[opaka_data](https://n-ducharmebarth-noaa.github.io/opal/reference/opaka_data.md)
for the main data object,
[opaka_parameters](https://n-ducharmebarth-noaa.github.io/opal/reference/opaka_parameters.md)
for starting parameter values.

## Examples

``` r
if (FALSE) { # \dontrun{
  data(opaka_lf)
  str(opaka_lf)
  # View sample data
  head(opaka_lf)
  # Get effective sample size per fishery-year
  opaka_lf[, .(total_n = sum(value)), by = .(year, fishery)]
  # Count observations by fleet
  opaka_lf[, uniqueN(paste(year, month)), by = fishery]
} # }
```
