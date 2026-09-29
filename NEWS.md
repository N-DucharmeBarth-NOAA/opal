# opal 0.0.4

- 2026/09/29: Add biological fit and posterior checks with versioned validation
  records, posterior report summaries, OSA residuals and dataset-level SDNR
  plots, composition and prior–posterior plots, objective profiles, resumable
  model grids, balanced posterior selection, and deterministic equilibrium
  MSY summaries. Include a simulated posterior and an assessment-tools tutorial.
  OSA covers lognormal indices and multinomial, Dirichlet, and
  Dirichlet–multinomial length and weight compositions, including latent-state
  integration. Dirichlet–multinomial OSA uses sequential beta-binomial densities;
  ordinary fitting and simulation retain the existing density implementation.
  Fix seasonal CPUE indexing, boundary diagnostics, invalid prior handling,
  and figure units. Scientific contract v5 changes CPUE predictions when
  `n_season > 1`; the single-season Opakapaka reference is numerically unchanged.
  Earlier saved model contracts require rebuilding under the new contract.
  Projection behaviour, selectivity time blocks, and close-kin components are
  outside this change.

- 2026/09/29: Make assessment objects the primary documentation workflow, add
  a runnable projection vignette, organise the function reference, and publish
  a separately labelled development pkgdown site with internal-link checks.
  Model numerics and function interfaces are unchanged.

- 2026/09/29: Expand tests for projection uncertainty, recruitment forecasts,
  selectivity, parameter diagnostics, and bundled data access. Tighten optimiser
  stopping tolerances in the Quickstart and projection examples while retaining
  the existing fit-check thresholds. Package model calculations are unchanged.

- 2026/09/29: Introduce the staged `opal_obj` workflow: configuration, two-pass
  fitting, MCMC, diagnostics, direct plots, projections with source identities,
  and portable `opal_save()`/`opal_read()` persistence. Preserve legacy fit
  construction and migrate saved fits with objective verification. Configuration
  changes invalidate dependent results; metadata updates and same-target refits
  preserve compatible samples. Model numerics and scientific contract v4 are unchanged.

- 2026/09/22 ([2ae7722](https://github.com/N-DucharmeBarth-NOAA/opal/commit/2ae772291b8d3deae5b3f25f172a99124d149091)): Use portable integrity in the bundled-fit examples so they run across R
  versions while retaining strict model compatibility and objective checks.
  Build the package website on pull requests to `dev` as well as `main`.
  Model numerics are unchanged.
- 2026/09/22 ([06f0d7c](https://github.com/N-DucharmeBarth-NOAA/opal/commit/06f0d7ce2cc1c868fe4873b289f244a78a69a451)): Constrain and penalise initial seasonal survival before seasonal
  compounding when fisheries overlap, and protect equilibrium recruitment from
  becoming negative under excessive initial fishing. The initialisation penalty
  is included in the model objective and reported as `lp_init_penalty`.
  Safely positive equilibria are unchanged. Previously saved fits use the
  earlier scientific model contract; the current contract is v4.
- 2026/09/16 ([c5dcdf1](https://github.com/N-DucharmeBarth-NOAA/opal/commit/c5dcdf1025e4be1647b4a0a9c200e4221df5691c)): Dirichlet-multinomial composition preparation now preserves half-up rounded
  effective sample sizes with largest-remainder integer counts. Existing saved
  fits use the previous scientific model contract.
- 2026/09/16 ([bac7a99](https://github.com/N-DucharmeBarth-NOAA/opal/commit/bac7a99dd2f1a416cfca6d1d14991c3167b4e5c5)): Fished equilibrium initialization now uses the same seasonal harvest-fraction
  survival convention as the population dynamics, eliminating drift for partial
  selectivity.
