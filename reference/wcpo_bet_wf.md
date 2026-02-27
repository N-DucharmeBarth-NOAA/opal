# West Central Pacific Ocean Bigeye Tuna (WCPO BET) Weight Frequency Data

Weight frequency composition data extracted from MFCL (Multifan-CL)
frequency files and prepared for RTMB (R Template Model Builder) stock
assessment modeling. Observations are organized by fishery, year, month,
and 1 kg weight bins.

## Usage

``` r
wcpo_bet_wf
```

## Format

A data.table with 7 columns: year (integer), month (numeric), ts
(integer), fishery (numeric), bin (integer), value (numeric), week
(numeric). Multiple rows per fishery-year-month (one per 1 kg weight
bin).

## Source

Extracted from MFCL (Multifan-CL) model frequency file outputs and
prepared for RTMB model input.

## Details

### Data workflow

**Preparation steps:**

1.  Extract weight composition data from MFCL frequency file using
    `wtfrq()` - Returns wide-format data with weight bins as columns (1,
    2, 3, ... 200 kg) - Observations grouped by fishery-year-month

2.  Reshape to long format (one row per fishery-year-month-weight bin
    combination)

3.  Merge with timestep table to add model timestep indexing consistent
    with catch/CPUE

4.  Save as long-format data for RTMB model input

### Output columns

- year:

  Calendar year (1952-2018). Integer.

- month:

  Fishing quarter month (2, 5, 8, 11 representing Q1-Q4 approx.).
  Numeric.

- ts:

  Model timestep number (1-268 for 67 years × 4 quarters). Integer.

- fishery:

  MFCL fleet index (1-15). Numeric.

- bin:

  Weight bin (kg). 1 kg bins: 1, 2, 3, ... 200. Integer.

- value:

  Observations in this weight bin. Numeric.

- week:

  Temporal indicator within quarter (typically 1 for aggregated data).
  Numeric.

### Notes

- MFCL data uses 1 kg weight bins (1, 2, 3, ... 200 kg)

- Data from `wtfrq()` represents observed weight frequency

- Sum of values across bins per fishery-year-month instance gives total
  input sample size

- Bin values are left as counts for RTMB to handle (divide by sample
  size if proportions needed)

## Examples

``` r
if (FALSE) { # \dontrun{
  data(wcpo_bet_wf)
  str(wcpo_bet_wf)
  # View sample data
  head(wcpo_bet_wf)
  # Get total sample size per fishery-year-month
  wcpo_bet_wf[, .(total_n = sum(value)), by = .(year, month, fishery)]
} # }
```
