# West Central Pacific Ocean Bigeye Tuna (WCPO BET) Assessment Parameters

A list of estimated parameter values for the WCPO bigeye tuna stock
assessment model. Contains growth parameters, recruitment deviations,
selectivity curves, and other model parameters fitted using RTMB.

## Usage

``` r
wcpo_bet_parameters
```

## Format

A list with 14 named elements containing estimated parameters from model
optimization.

## Details

### Parameter structure (14 list elements)

**Stock parameters:**

- log_B0:

  Initial spawning biomass (log-scale): 12

- log_h:

  Steepness parameter (log-scale): -0.0513

- log_sigma_r:

  Recruitment standard deviation (log-scale): -0.511

**Observation model parameters:**

- log_cpue_q:

  CPUE catchability coefficient (log-scale): 0

- cpue_creep:

  CPUE creep adjustment: 0

- log_cpue_tau:

  CPUE log-standard deviation (log-scale): -Inf

- log_cpue_omega:

  CPUE process error (log-scale): 0

**Recruitment deviations:**

- rdev_y:

  Recruitment deviations by year (268 elements, quarterly timesteps)

**Selectivity parameters:**

- par_sel:

  Selectivity parameters matrix (15 fisheries × 6 columns). Rows
  represent different fisheries/fleets, columns represent selectivity
  curve parameters (e.g., inflection point, slope for logistic curves or
  parameters for double-normal curves).

**Growth parameters:**

- log_L1:

  Log of length at age 1 (von Bertalanffy): 3.431379

- log_L2:

  Log of asymptotic length (von Bertalanffy): 5.03333

- log_k:

  Log of growth rate parameter (von Bertalanffy): -2.32024

- log_CV1:

  Log of CV at young ages: -1.826289

- log_CV2:

  Log of CV at old ages: -2.230151

## See also

[wcpo_bet_data](https://n-ducharmebarth-noaa.github.io/opal/reference/wcpo_bet_data.md)
for the corresponding data object used in model fitting.

## Examples

``` r
if (FALSE) { # \dontrun{
  data(wcpo_bet_parameters)
  str(wcpo_bet_parameters)
  # Access growth parameters
  exp(wcpo_bet_parameters$log_L1)  # L1 on real scale
  exp(wcpo_bet_parameters$log_L2)  # Linf on real scale
  exp(wcpo_bet_parameters$log_k)   # k on real scale
} # }
```
