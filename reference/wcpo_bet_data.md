# West Central Pacific Ocean Bigeye Tuna (WCPO BET) Assessment Data

A comprehensive list containing biological parameters, catch and CPUE
data, and Bayesian prior specifications prepared for RTMB (R Template
Model Builder) stock assessment modeling of western central Pacific
bigeye tuna.

## Usage

``` r
wcpo_bet_data
```

## Format

A list with 25 named elements combining dimensions, biological
parameters, observed data (catch and CPUE), and Bayesian priors for RTMB
model initialization.

## Source

Prepared from MFCL (Multifan-CL) model outputs and frequency files, with
biological parameters and priors derived from assessment modeling.

## Details

### Data structure (25 list elements)

**Age and time dimensions:**

- age_a:

  Vector of integer ages (1:40). Length = 40.

- n_age:

  Number of age classes (40).

- n_season:

  Number of seasons per year (1).

- n_fishery:

  Number of fisheries/fleets (15).

- n_year:

  Number of years in model (268, equivalent to 67 years × 4 quarters).

- first_yr:

  First year designation in model (1, representing 1952).

- last_yr:

  Last year designation in model (268, representing 2018).

- years:

  Vector of model years (1:268).

**Length structure:**

- len_bin_start:

  First length bin lower edge (10 cm).

- len_bin_width:

  Width of each length bin (2 cm).

- n_len:

  Total number of length bins (95).

**Catch data:**

- first_yr_catch:

  First year with catch data (1).

- catch_units_f:

  Numeric vector of unit codes by fishery (15 elements). 1 = metric tons
  (fisheries 8-14), 2 = thousands of fish (fisheries 1-7, 15).

- catch_obs_ysf:

  Three-dimensional array of observed catch (268 years × 1 season × 15
  fisheries). Dimensions follow (year, season, fishery). Units follow
  catch_units_f.

**CPUE data:**

- cpue_switch:

  Indicator for CPUE inclusion (1 = included, typically for fishery 15).

- cpue_data:

  Tibble with 268 rows and 8 columns containing standardized CPUE
  observations:

  - year: Calendar year (1952-2018)

  - month: Quarter month (2, 5, 8, 11 for Q1-Q4)

  - ts: Sequential timestep (1-268)

  - fishery: Fleet identifier (usually 15 for survey fleet)

  - metric: Data type ("cpue")

  - units: Unit code (typically 2 for fishery 15)

  - value: Normalized CPUE observations

  - se: Log-space standard error of CPUE observations

**Biological parameters:**

- lw_a:

  Length-weight allometric coefficient (6.48e-05).

- lw_b:

  Length-weight allometric exponent (2.78).

- maturity:

  Maturity at length vector (95 elements). Ogive-based maturity
  probabilities by length bin.

- fecundity:

  Fecundity at length vector (95 elements). Relative fecundity by length
  bin.

- M:

  Natural mortality at age vector (40 elements). Age-specific
  instantaneous mortality rates.

**Age bounds:**

- A1:

  Minimum age for credibility interval (1).

- A2:

  Maximum age for credibility interval (40).

**Selectivity:**

- sel_type_f:

  Selectivity function type by fishery (15 elements). Integer codes
  representing selectivity curve shapes (e.g., 2 = double-normal).

**Priors:**

- priors:

  List of 8 Bayesian prior specifications, each containing:

  - type: Prior distribution type ("normal", etc.)

  - par1: First parameter (mean or location)

  - par2: Second parameter (std dev or scale)

  - index: Parameter index in the model

  Priors cover: log_B0, log_cpue_q, par_sel, log_L1, log_L2, log_k,
  log_CV1, log_CV2.

## Examples

``` r
if (FALSE) { # \dontrun{
  data(wcpo_bet_data)
  str(wcpo_bet_data)
  # Access catch data
  head(wcpo_bet_data$cpue_data)
  # Get biological parameters
  wcpo_bet_data$M
} # }
```
