# Changelog

## opal 0.1.1

- 2026/10/15 ([PR
  \#143](https://github.com/N-DucharmeBarth-NOAA/opal/pull/143)): Fix
  fit validation for random-effect models by calculating the marginal
  Hessian with [`optimHess()`](https://rdrr.io/r/stats/optim.html)
  instead of calling RTMB’s unsupported exact Hessian. Converged,
  biologically valid random-effect fits can now pass validation and
  proceed to grid sampling. Diagnostic version v3 invalidates old checks
  and grid checkpoints; rerun
  [`opal_check()`](https://n-ducharmebarth-noaa.github.io/opal/reference/opal_check.md)
  on existing v5 objects. This diagnostic fix does not change model
  objectives or fitted estimates.

- Quantify the existing seasonal CPUE correction in reproducible paired
  simulations: maximum prediction changes at fixed parameters are 8.24%
  and 12.28%, and refitted B0 changes are +4.99% and +11.71%, for the
  specified two- and four-season examples. These are illustrative, not
  real-stock impact estimates or general correction factors.
  Full-coverage examples fail under the legacy indexing but evaluate
  under v5. The single-season control is unchanged, and all bundled
  stock inputs use one season. See the seasonal CPUE impact article for
  observation coverage, diagnostics, and instructions for rebuilding and
  refitting saved contract-v4-or-earlier assessments.

- 2026/10/13
  ([e1fbdbd](https://github.com/N-DucharmeBarth-NOAA/opal/commit/e1fbdbd78eeebff20b0a13e55cb3c8e8148408cc)):
  Add optional SS3 GrowthModel 3 cumulative age-specific K multipliers
  and plus-group length adjustment, credited to Stock Synthesis (Methot
  and Wetzel 2013). Existing model configurations retain their exact
  growth calculation and numerical contract.

## opal 0.1.0

This release makes the staged `opal_obj` workflow the primary way to run
an assessment, from configuration and fitting through MCMC, diagnostics,
projections and saved results. It also adds posterior and fit-checking
tools. Scientific contract v5 changes CPUE predictions when
`n_season > 1`, so saved models from earlier contracts must be rebuilt.

- 2026/09/29
  ([0ea5022](https://github.com/N-DucharmeBarth-NOAA/opal/commit/0ea50229de2a576f7a139c39f9c46f73022c3a9f)):
  Add biological fit and posterior checks with versioned validation
  records, posterior report summaries, OSA residuals and dataset-level
  SDNR plots, composition and prior–posterior plots, objective profiles,
  resumable model grids, balanced posterior selection, and deterministic
  equilibrium MSY summaries. Include a simulated posterior and an
  assessment-tools tutorial. OSA covers lognormal indices and
  multinomial, Dirichlet, and Dirichlet–multinomial length and weight
  compositions, including latent-state integration.
  Dirichlet–multinomial OSA uses sequential beta-binomial densities;
  ordinary fitting and simulation retain the existing density
  implementation. Fix seasonal CPUE indexing, boundary diagnostics,
  invalid prior handling, and figure units. Scientific contract v5
  changes CPUE predictions when `n_season > 1`; the single-season
  Opakapaka reference is numerically unchanged. Earlier saved model
  contracts require rebuilding under the new contract. Projection
  behaviour, selectivity time blocks, and close-kin components are
  outside this change.

- 2026/09/29
  ([fba96ea](https://github.com/N-DucharmeBarth-NOAA/opal/commit/fba96ea2968f6d3016d78bb44b55c862e3a2ee57)):
  Make assessment objects the primary documentation workflow, add a
  runnable projection vignette, organise the function reference, and
  publish a separately labelled development pkgdown site with
  internal-link checks. Model numerics and function interfaces are
  unchanged.

- 2026/09/29
  ([3475a2e](https://github.com/N-DucharmeBarth-NOAA/opal/commit/3475a2ea09640ee10d8a9477b44d6227bf61a2db)):
  Expand tests for projection uncertainty, recruitment forecasts,
  selectivity, parameter diagnostics, and bundled data access. Tighten
  optimiser stopping tolerances in the Quickstart and projection
  examples while retaining the existing fit-check thresholds. Package
  model calculations are unchanged.

- 2026/09/29
  ([17ae7a6](https://github.com/N-DucharmeBarth-NOAA/opal/commit/17ae7a69e5c8450c1dd5691560dd01048708e9f4)):
  Introduce the staged `opal_obj` workflow: configuration, two-pass
  fitting, MCMC, diagnostics, direct plots, projections with source
  identities, and portable
  [`opal_save()`](https://n-ducharmebarth-noaa.github.io/opal/reference/opal_io.md)/[`opal_read()`](https://n-ducharmebarth-noaa.github.io/opal/reference/opal_io.md)
  persistence. Preserve legacy fit construction and migrate saved fits
  with objective verification. Configuration changes invalidate
  dependent results; metadata updates and same-target refits preserve
  compatible samples. Model numerics and scientific contract v4 are
  unchanged.

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
