# 'Opakapaka Stock Assessment Parameters

A list of parameter values for the opakapaka stock assessment model,
extracted from a fitted Stock Synthesis 3 (SS3) model and converted to
opal conventions. Contains growth parameters, recruitment deviations,
selectivity curves, and other model parameters suitable as starting
values for RTMB optimization.

## Usage

``` r
opaka_parameters
```

## Format

A list with 14 named elements containing parameters extracted from SS3
model optimization and converted to opal conventions.

## Source

Extracted from Stock Synthesis 3 parameter file (`ss.par`) and control
file (`control.ss_new`). Growth parameters converted from standard VB to
Schnute parameterization with A1 = 1, A2 = n_age.

## Details

### Parameter structure (14 list elements)

**Stock parameters:**

- log_B0:

  Log initial recruitment/biomass (SS3's LN(R0)): 5.656

- log_h:

  Steepness parameter (log-scale): log(0.76) = -0.274

- log_sigma_r:

  Recruitment standard deviation (log-scale): log(0.52) = -0.654

**Observation model parameters:**

- log_cpue_q:

  CPUE catchability coefficients (log-scale). Vector of length 2: -3.773
  (Comm fleet), -6.253 (ResFish survey).

- cpue_creep:

  CPUE creep adjustment: 0

- log_cpue_tau:

  CPUE observation error (log-scale): -Inf (disabled)

- log_cpue_omega:

  CPUE process error (log-scale): 0

**Recruitment deviations:**

- rdev_y:

  Recruitment deviations by year (75 elements, annual timesteps
  1949-2023). Extracted from SS3 estimated main recruitment deviations.

**Selectivity parameters:**

- par_sel:

  Selectivity parameters matrix (3 fisheries x 6 columns).

  - Fleet 1 (Comm): Logistic. Cols 1-2 = inflection (36.1 cm), 95\\

  - Fleet 2 (Non_comm): Logistic (fixed). Cols 1-2 = inflection (40 cm),
    95\\

  - Fleet 3 (ResFish): Double-normal (SS3 pattern 24). Cols 1-6 = peak
    (20.5), top_logit (-1.06), ascend_se (-0.34), descend_se (4.11),
    start_logit (-999), end_logit (-1.40).

**Growth parameters (Schnute parameterization, log-scale):**

- log_L1:

  Log of length at reference age A1 = 1 (opal internal), corresponding
  to SS3 age 0: log(6.0) = 1.792

- log_L2:

  Log of length at reference age A2 = 44 (opal internal), corresponding
  to SS3 age 43: log(67.498) = 4.212

- log_k:

  Log of von Bertalanffy growth coefficient: log(0.242) = -1.419

- log_CV1:

  Log of CV at young ages: log(0.085) = -2.465

- log_CV2:

  Log of CV at old ages: log(0.085) = -2.465. Equal to log_CV1, giving
  constant CV across all lengths.

## See also

[opaka_data](https://n-ducharmebarth-noaa.github.io/opal/reference/opaka_data.md)
for the corresponding data object,
[opaka_lf](https://n-ducharmebarth-noaa.github.io/opal/reference/opaka_lf.md)
for length frequency composition data.

## Examples

``` r
if (FALSE) { # \dontrun{
  data(opaka_parameters)
  str(opaka_parameters)
  # Access growth parameters on real scale
  exp(opaka_parameters$log_L1)   # L1 = 6.0 cm
  exp(opaka_parameters$log_L2)   # L2 = 67.5 cm
  exp(opaka_parameters$log_k)    # k = 0.242
  # Steepness
  exp(opaka_parameters$log_h)    # h = 0.76
} # }
```
