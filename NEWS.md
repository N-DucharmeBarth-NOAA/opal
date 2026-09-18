# opal 0.0.4

- 2026/09/16 ([c5dcdf1](https://github.com/N-DucharmeBarth-NOAA/opal/commit/c5dcdf1025e4be1647b4a0a9c200e4221df5691c)): Dirichlet-multinomial composition preparation now preserves half-up rounded
  effective sample sizes with largest-remainder integer counts. Existing saved
  fits use the previous scientific model contract.
- 2026/09/16 ([bac7a99](https://github.com/N-DucharmeBarth-NOAA/opal/commit/bac7a99dd2f1a416cfca6d1d14991c3167b4e5c5)): Fished equilibrium initialization now uses the same seasonal harvest-fraction
  survival convention as the population dynamics, eliminating drift for partial
  selectivity.