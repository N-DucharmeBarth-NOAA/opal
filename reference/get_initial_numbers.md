# Initial numbers and Beverton-Holt parameters

Computes the initial equilibrium numbers-at-age, unfished recruitment
(R0), and Beverton-Holt stock-recruitment parameters.

## Usage

``` r
get_initial_numbers(B0, h, M_a, spawning_potential_a)
```

## Arguments

- B0:

  Unfished spawning biomass.

- h:

  Beverton-Holt steepness parameter.

- M_a:

  a `vector` of natural mortality at age.

- spawning_potential_a:

  a `vector` of spawning potential at age (maturity × fecundity).

## Value

A list containing:

- Ninit:

  Initial numbers-at-age (vector).

- R0:

  Unfished recruitment (scalar).

- alpha:

  BH alpha parameter.

- beta:

  BH beta parameter.
