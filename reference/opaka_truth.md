# Opakapaka SS3 OM/EM truth and EM output (extracted)

A small list containing Stock Synthesis (SS3) output extracted from an
operating model (OM) and an estimation model (EM) pair using
`r4ss::SS_output()`. Intended as a truth / comparison object for
simulation-experiment workflows and diagnostics.

## Usage

``` r
opaka_truth
```

## Format

A list with three named elements, `om`, `em`, and `em_no_platoon`, each
containing the SS3 `sprseries` and `cpue` outputs as data.frames
suitable for plotting and comparison of true vs estimated indices and
spawning biomass trajectories.

## Source

Produced by calling `r4ss::SS_output()` on OM and EM model directories
and collecting the `sprseries` and `cpue` components.

## Details

The object is a list with three elements: `om`, `em`, and
`em_no_platoon`. Each element is a sub-list that contains `sprseries`
(spawning biomass time series) and `cpue` (standardized index output).
The contained data frames have the following structures:

- om\$sprseries:

  A data.frame with 100 rows and 2 columns:

  Yr

  :   Integer year (e.g., 1949...)

  SSB

  :   Numeric spawning stock biomass

- om\$cpue:

  A data.frame with 107 rows and 7 columns:

  Yr

  :   Integer year

  Fleet

  :   Integer fleet id

  Fleet_name

  :   Character fleet label

  Vuln_bio

  :   Numeric vulnerable biomass predicted by SS3

  Obs

  :   Observed CPUE (integer in OM)

  Exp

  :   Expected CPUE (numeric)

  SE

  :   Standard error (numeric)

- em\$sprseries:

  A data.frame with 100 rows and 2 columns: same format as
  `om$sprseries`.

- em\$cpue:

  A data.frame with 107 rows and 7 columns: same column names as
  `om$cpue`, with `Obs` shown as numeric in the EM output shown here.

- em_no_platoon\$sprseries:

  A data.frame with 100 rows and 2 columns: same format as
  `om$sprseries`.

- em_no_platoon\$cpue:

  A data.frame with 107 rows and 7 columns: same column names as
  `om$cpue`, with `Obs` shown as numeric in the EM output shown here.

## Examples

``` r
if (FALSE) { # \dontrun{
  data(opaka_truth)
  str(opaka_truth)
  plot(opaka_truth$om$sprseries$Yr, opaka_truth$om$sprseries$SSB, type = "l")
  lines(opaka_truth$em$sprseries$Yr, opaka_truth$em$sprseries$SSB, col = "red")
  lines(opaka_truth$em_no_platoon$sprseries$Yr,
        opaka_truth$em_no_platoon$sprseries$SSB, col = "blue")
  # Compare CPUE
  head(opaka_truth$om$cpue)
  head(opaka_truth$em$cpue)
} # }
```
