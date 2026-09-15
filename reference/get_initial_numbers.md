# Initial numbers and Beverton-Holt parameters

Computes the initial equilibrium numbers-at-age, unfished recruitment
(R0), and Beverton-Holt stock-recruitment parameters.

## Usage

``` r
get_initial_numbers(
  B0,
  h,
  M_a,
  spawning_potential_a,
  init_F_f = NULL,
  sel_fa = NULL,
  init_rdev_a = NULL,
  sigma_r = 0.6,
  init_bias_adj_a = NULL
)
```

## Arguments

- B0:

  Unfished spawning biomass.

- h:

  Beverton-Holt steepness parameter.

- M_a:

  a `vector` of natural mortality at age.

- spawning_potential_a:

  a `vector` of spawning potential at age (maturity x fecundity).

- init_F_f:

  an optional `vector` of initial fishing mortality by fishery.

- sel_fa:

  an optional matrix of selectivity-at-age with dimensions
  `[n_fishery, n_age]`.

- init_rdev_a:

  an optional `vector` of initial age deviations.

- sigma_r:

  recruitment standard deviation used in lognormal correction.

- init_bias_adj_a:

  an optional `vector` of bias adjustment scalars for initial age
  deviations. Defaults to zero so fixed zero initial deviations do not
  alter the equilibrium initial age structure.

## Value

A list containing:

- Ninit:

  Initial numbers-at-age (vector).

- Ninit0:

  Initial unfished numbers-at-age (vector).

- R0:

  Unfished recruitment (scalar).

- alpha:

  BH alpha parameter.

- beta:

  BH beta parameter.
