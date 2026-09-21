# Changelog

## opal 0.0.4

- 2026/09/22
  ([2ae7722](https://github.com/N-DucharmeBarth-NOAA/opal/commit/2ae772291b8d3deae5b3f25f172a99124d149091)):
  Use portable integrity in the bundled-fit examples so they run across
  R versions while retaining strict model compatibility and objective
  checks. Build the package website on pull requests to `dev` as well as
  `main`. Model numerics are unchanged.
- 2026/09/22
  ([06f0d7c](https://github.com/N-DucharmeBarth-NOAA/opal/commit/06f0d7ce2cc1c868fe4873b289f244a78a69a451)):
  Constrain and penalise initial seasonal survival before seasonal
  compounding when fisheries overlap, and protect equilibrium
  recruitment from becoming negative under excessive initial fishing.
  The initialisation penalty is included in the model objective and
  reported as `lp_init_penalty`. Safely positive equilibria are
  unchanged. Previously saved fits use the earlier scientific model
  contract; the current contract is v4.
- 2026/09/16
  ([c5dcdf1](https://github.com/N-DucharmeBarth-NOAA/opal/commit/c5dcdf1025e4be1647b4a0a9c200e4221df5691c)):
  Dirichlet-multinomial composition preparation now preserves half-up
  rounded effective sample sizes with largest-remainder integer counts.
  Existing saved fits use the previous scientific model contract.
- 2026/09/16
  ([bac7a99](https://github.com/N-DucharmeBarth-NOAA/opal/commit/bac7a99dd2f1a416cfca6d1d14991c3167b4e5c5)):
  Fished equilibrium initialization now uses the same seasonal
  harvest-fraction survival convention as the population dynamics,
  eliminating drift for partial selectivity.
