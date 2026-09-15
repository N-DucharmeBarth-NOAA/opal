# Compute unfished equilibrium quantities from natural mortality only

Internal helper used by
[`get_initial_numbers`](https://n-ducharmebarth-noaa.github.io/opal/reference/get_initial_numbers.md)
to calculate unfished survivorship-per-recruit, \\R_0\\, and
Beverton-Holt parameters.

## Usage

``` r
get_unfished_init(B0, h, M_a, spawning_potential_a)
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

## Value

A list with `rel_N`, `R0`, `alpha`, and `beta`.
