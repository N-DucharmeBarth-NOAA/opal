# 'Opakapaka Stock Assessment Data

A comprehensive list containing biological parameters, catch and CPUE
data, and Bayesian prior specifications prepared for RTMB (R Template
Model Builder) stock assessment modeling. Extracted from a Stock
Synthesis 3 (SS3) model and converted to opal data conventions.

## Usage

``` r
opaka_data
```

## Format

A list with 22 named elements combining dimensions, biological
parameters, observed data (catch and CPUE), and Bayesian priors for RTMB
model initialization.

## Source

Extracted from Stock Synthesis 3 (SS3 v3.30.19) model input files
(`data_echo.ss_new`, `control.ss_new`). Historical data only
(1949-2023), forecast period excluded.

## Details

### Data structure (22 list elements)

**Age and time dimensions:**

- age_a:

  Vector of integer ages (1:44). Length = 44. Internal opal ages mapping
  to SS3 ages 0-43 (opal age `i` = SS3 age `i - 1`).

- n_age:

  Number of age classes (44).

- n_season:

  Number of seasons per year (1).

- n_fishery:

  Number of fisheries/fleets (3).

- n_year:

  Number of years in model (75, representing 1949-2023).

- first_yr:

  First year designation in model (1, representing 1949).

- last_yr:

  Last year designation in model (75, representing 2023).

- years:

  Vector of model years (1:75).

**Length structure:**

- len_bin_start:

  First length bin lower edge (5 cm).

- len_bin_width:

  Width of each length bin (5 cm).

- n_len:

  Total number of length bins (17).

**Catch data:**

- first_yr_catch:

  First year with catch data (1).

- catch_units_f:

  Numeric vector of unit codes by fishery (3 elements). All fisheries
  use units = 1 (biomass).

- catch_obs_ysf:

  Three-dimensional array of observed catch (75 years x 1 season x 3
  fisheries). Dimensions follow (year, season, fishery). Fleet 3
  (ResFish survey) has zero catch.

**CPUE data:**

- cpue_switch:

  Indicator for CPUE inclusion (1 = included).

- cpue_data:

  Data.table with 82 rows and 8 columns containing CPUE observations
  from two fleets:

  - year: Calendar year (1949-2023)

  - month: Observation month (7 for Comm fleet, 1 for ResFish)

  - ts: Sequential timestep (1-75)

  - fishery: Fleet identifier (1 = Comm, 3 = ResFish)

  - index: Integer index identifier (1 for Comm, 2 for ResFish)

  - metric: Data type ("cpue")

  - units: Unit code (1 for biomass)

  - value: CPUE observations

  - se: Log-space standard error (0.2 for Comm, 0.15 for ResFish)

  Comm fleet (fleet 1) covers 1949-2023 (75 obs). ResFish survey
  (fleet 3) covers 2017-2023 (7 obs).

**Biological parameters:**

- lw_a:

  Length-weight allometric coefficient (1.75e-05).

- lw_b:

  Length-weight allometric exponent (2.99).

- maturity:

  Maturity at length vector (17 elements). Logistic maturity
  probabilities at each data length bin midpoint (Mat50 = 40.7 cm, slope
  = -2.26).

- fecundity:

  Fecundity at length vector (17 elements). Fecundity equals weight at
  each length bin midpoint (fecundity option: eggs = Wt).

- M:

  Natural mortality at age vector (44 elements). Constant instantaneous
  rate (0.135) across all ages.

**Growth reference ages (Schnute parameterization):**

- A1:

  Reference age for L1 (1). Corresponds to SS3 age 0.

- A2:

  Reference age for L2 (44 = n_age). Corresponds to SS3 age 43.

**Selectivity:**

- sel_type_f:

  Selectivity function type by fishery (3 elements). Integer codes: 1 =
  logistic, 2 = double-normal, 3 = double Richards.

**Priors:**

- priors:

  List of 8 Bayesian prior specifications, each containing:

  - type: Prior distribution type ("none" = uninformative)

  - par1: First parameter (NA when type = "none")

  - par2: Second parameter (NA when type = "none")

  - index: Parameter index in the model

  Priors cover: log_B0, log_cpue_q, par_sel, log_L1, log_L2, log_k,
  log_CV1, log_CV2.

### Fleet definitions

- Fleet 1 (Comm):

  Commercial fishery. Catch fleet with logistic selectivity. CPUE and
  length composition data available 1949-2023.

- Fleet 2 (Non_comm):

  Non-commercial fishery. Catch fleet with logistic selectivity (fixed
  parameters). Catch data only, no CPUE or compositions.

- Fleet 3 (ResFish):

  Research fishery survey. Survey fleet with double-normal selectivity.
  CPUE and length composition data available 2017-2023.

### Growth conversion

Growth parameters were converted from the SS3 standard von Bertalanffy
parameterization (L_at_Amin = 6.0 at age 0, Linf = 67.5 at age 999, k =
0.242) to the opal Schnute parameterization with A1 = 1 and A2 = n_age.
The conversion is exact to machine precision (~1e-14). The VB growth
coefficient k is unchanged between parameterizations.

## See also

[opaka_parameters](https://n-ducharmebarth-noaa.github.io/opal/reference/opaka_parameters.md)
for starting parameter values,
[opaka_lf](https://n-ducharmebarth-noaa.github.io/opal/reference/opaka_lf.md)
for length frequency composition data.

## Examples

``` r
if (FALSE) { # \dontrun{
  data(opaka_data)
  str(opaka_data)
  # Access catch data
  head(opaka_data$cpue_data)
  # Get biological parameters
  opaka_data$M
  # Check dimensions
  cat("Ages:", opaka_data$n_age, "\n")
  cat("Years:", opaka_data$n_year, "\n")
  cat("Fleets:", opaka_data$n_fishery, "\n")
  cat("Length bins:", opaka_data$n_len, "\n")
} # }
```
