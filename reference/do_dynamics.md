# Population dynamics

Runs the core age- and season-structured population dynamics loop for
bigeye tuna. Starts from initial equilibrium numbers (derived from B0
and h), applies seasonal harvest, natural mortality, spawning, and
recruitment (Beverton-Holt with log-normal deviates), and computes
predicted catches and harvest rates.

## Usage

``` r
do_dynamics(
  data,
  parameters,
  B0,
  R0,
  alpha,
  beta,
  h = 0.95,
  sigma_r = 0.6,
  M_a,
  spawning_potential_a,
  weight_fya,
  init_number_a,
  sel_fya
)
```

## Arguments

- data:

  A `list` of model data. Must contain at minimum: `first_yr`,
  `first_yr_catch`, `n_year`, `n_season`, `n_fishery`, `n_age`,
  `catch_obs_ysf`, `catch_units_f`.

- parameters:

  A `list` of model parameters. Must contain at minimum: `rdev_y`.

- B0:

  Numeric. Unfished equilibrium spawning biomass.

- R0:

  Numeric. Unfished equilibrium recruitment.

- alpha:

  Numeric. Beverton-Holt stock-recruitment alpha parameter.

- beta:

  Numeric. Beverton-Holt stock-recruitment beta parameter.

- h:

  Numeric (0.2–1). Steepness of the Beverton-Holt stock-recruitment
  relationship.

- sigma_r:

  Numeric \> 0. Standard deviation of log recruitment deviations.

- M_a:

  Numeric vector of length `n_age`. Natural mortality at age. Passed
  explicitly so AD gradients propagate if M is ever estimated.

- spawning_potential_a:

  Numeric vector of length `n_age`. Spawning potential at age (maturity
  × fecundity). Passed explicitly so AD gradients propagate if growth is
  ever estimated.

- weight_fya:

  Numeric array `[n_fishery, n_year, n_age]`. Mean weight at age by
  fishery and year. Passed explicitly so AD gradients propagate if
  growth is ever estimated.

- init_number_a:

  Numeric vector of length `n_age`. Initial equilibrium numbers-at-age
  (from
  [`get_initial_numbers`](https://n-ducharmebarth-noaa.github.io/opal/reference/get_initial_numbers.md)).

- sel_fya:

  Numeric array `[n_fishery, n_year, n_age]`. Fishery-specific
  selectivity at age by year (from
  [`get_selectivity`](https://n-ducharmebarth-noaa.github.io/opal/reference/get_selectivity.md)).

## Value

A named list with:

- number_ysa:

  Numbers-at-age array `[n_year+1, n_season, n_age]`.

- lp_penalty:

  Total penalty from
  [`posfun`](https://n-ducharmebarth-noaa.github.io/opal/reference/posfun.md)
  (harvest rate constraints).

## Details

All derived biology arrays (`M_a`, `spawning_potential_a`, `weight_fya`)
are passed as explicit arguments rather than read from `data`. This
ensures that AD gradients propagate correctly if any of these quantities
carry estimated parameters in the future (e.g., growth parameters
estimated via the PLA, or natural mortality via the Lorenzen equation).
