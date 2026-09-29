# Calculate deterministic equilibrium reference points

Calculates MSY under a fixed fleet allocation, selectivity, and biology.
This is an equilibrium calculation, independent of the projection
engine.

## Usage

``` r
opal_msy(
  x,
  fleet_weights,
  year = x$data$n_year,
  uncertainty = c("fit", "mcmc"),
  draws = NULL,
  u_max = 1,
  grid_size = 201L,
  name = "msy"
)
```

## Arguments

- x:

  A fitted or sampled `opal_obj`.

- fleet_weights:

  Non-negative fleet allocation weights, normalised to sum to one.
  Required because the fleet mix is a scientific choice.

- year:

  Model-year index supplying selectivity and weight at age.

- uncertainty:

  Use the fitted point or complete joint MCMC draws.

- draws:

  Optional retained draw indices when `uncertainty = "mcmc"`.

- u_max:

  Maximum scalar seasonal harvest intensity, in `(0, 1]`.

- grid_size:

  Number of harvest intensities used to bracket maxima.

- name:

  Name of the stored result.

## Value

An updated object.
[`opal_derived()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_derived.md)
returns per-draw MSY, spawning output at MSY, recruitment at MSY,
seasonal intensity at MSY, depletion at MSY, terminal spawning output
relative to MSY, and status probabilities.

## Details

Fleet harvest fractions are `u * fleet_weights[f] * selectivity[f,a]`.
Natural mortality follows seasonal fishing, matching historical
dynamics. Yield is total catch in the units of `weight_fya_mod`,
including fleets whose input catches use numbers. Beverton-Holt
equilibrium recruitment is solved analytically without recruitment
deviations or bias corrections. No positive continuation is accepted as
a biological equilibrium.

These are deterministic, fixed-biology reference points, not stochastic
MSY or management advice. Spawning potential at recruitment age must be
zero to match the dynamics' spawning-before-recruitment convention. An
upper-bound optimum is flagged and should not be treated as a resolved
MSY. Posterior calculations use each draw's own biology and do not
discard invalid draws. Equal draw weights are used in status
probabilities.
