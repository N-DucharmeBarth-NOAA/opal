# Contributing to opal

Thank you for your interest in contributing to opal. This guide
describes how to propose changes, set up a development environment, and
the checks a pull request must pass. Contributions of all sizes are
welcome, including bug reports, documentation fixes, tests, and new
model features.

opal is under active development and interfaces may change. If you are
unsure whether a change fits the project, open an issue first.

## Where to contribute

| Contribution | Location |
|----|----|
| Package code, tests, vignettes, and function documentation | [N-DucharmeBarth-NOAA/opal](https://github.com/N-DucharmeBarth-NOAA/opal) |
| Project website, WCPFC documents, and case-study write-ups | [N-DucharmeBarth-NOAA/opal-documentation](https://github.com/N-DucharmeBarth-NOAA/opal-documentation) (published at <https://connect.fisheries.noaa.gov/opal/>) |
| Bug reports, feature requests, and design discussion | [opal issues](https://github.com/N-DucharmeBarth-NOAA/opal/issues) |

## Issues

Before opening an issue, search existing issues to avoid duplicates.

A **bug report** should include a minimal reproducible example, the
output of [`sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html),
and the opal version or commit. For numerical problems, state what you
expected, what you observed, and how you checked it (for example,
against an SS3 or MULTIFAN-CL reference model).

A **feature request or substantial change** should be discussed in an
issue before work starts. Issues intended for implementation, including
those assigned to coding agents, should be self-contained, with:

- the files and functions affected, including exact function signatures;
- acceptance criteria;
- whether the change is expected to alter model numerics (see [Numerical
  changes](#numerical-changes));
- backward-compatibility implications and expected effects on existing
  tests;
- notes on AD safety for any code that operates on estimated quantities.

## Branches and pull requests

- **`main`** is the stable release branch.
- **`dev`** is the active development branch.

Create a feature branch from `dev` and open pull requests against `dev`.
Keep each pull request focused on one change. A pull request should:

1.  describe what changed and why, and link the relevant issue;
2.  state whether model numerics change, and if so follow the
    regeneration protocol below;
3.  add or update tests for new or changed behaviour;
4.  update roxygen documentation and vignettes when exported function
    signatures or behaviour change;
5.  add a `NEWS.md` entry for user-visible or numerical changes.

## Development setup

The package environment is pinned with `renv`. From the repository root:

``` r

renv::restore()
```

`compResidual`, used for one-step-ahead residuals on composition data,
requires the multivariate OSA headers to be installed into `TMB`. If you
need it, install `TMB` first, then:

``` r

TMB:::install.contrib(
  "https://github.com/vtrijoulet/OSA_multivariate_dists/archive/main.zip"
)
```

`SparseNUTS` is not on CRAN. If it fails to restore, see the
installation notes in the
[README](https://n-ducharmebarth-noaa.github.io/README.html#installation).

## Development workflow

- Load the package with
  [`devtools::load_all()`](https://devtools.r-lib.org/reference/load_all.html)
  rather than [`source()`](https://rdrr.io/r/base/source.html).
- After editing roxygen comments, run
  [`devtools::document()`](https://devtools.r-lib.org/reference/document.html).
  The `document.yaml` workflow also regenerates `man/` and `NAMESPACE`
  on pushes that touch `R/`, and commits the result, so pull before
  continuing work on a branch.
- Run
  [`devtools::test()`](https://devtools.r-lib.org/reference/test.html)
  before opening a pull request. Run
  [`devtools::check()`](https://devtools.r-lib.org/reference/check.html)
  for changes that affect exports, dependencies, or documentation.
- Prefer small, targeted changes. Refactors that touch several files or
  array dimensions should be proposed in an issue first, with their
  effect on test fixtures described.

## Code guidelines

### AD safety

Model code runs on RTMB’s automatic-differentiation tape. Any function
that may receive AD values must:

- declare the relevant overloads at the top of the function, for example
  `"[<-" <- ADoverload("[<-")` (and `"c" <- ADoverload("c")` where
  needed);
- avoid base R coercions that strip the `advector` class, such as
  [`matrix()`](https://rdrr.io/r/base/matrix.html) or
  [`as.vector()`](https://rdrr.io/r/base/vector.html) on AD values. Fix
  type loss where it occurs, not in the receiving function;
- take derived biological quantities (for example `M_a`,
  `spawning_potential_a`, `weight_fya`) as explicit arguments rather
  than reading them from `data`, so gradients propagate when the
  underlying parameters are estimated.

After a change, confirm that `obj$fn(obj$par)` and `obj$gr(obj$par)`
return finite values.

### Exported functions

Exported functions are documented in `man/` and used in vignettes and
tests. Changing a signature requires updating the roxygen documentation,
`NAMESPACE`, affected vignettes, and tests in the same pull request.

### Model components

- **Selectivity**: changes to selectivity forms must update the
  `par_sel` dimensions,
  [`get_map()`](https://n-ducharmebarth-noaa.github.io/opal/reference/get_map.md),
  and
  [`get_bounds()`](https://n-ducharmebarth-noaa.github.io/opal/reference/get_bounds.md)
  together.
- **Composition likelihoods**: data flow from raw data through
  [`prep_lf_data()`](https://n-ducharmebarth-noaa.github.io/opal/reference/prep_lf_data.md)
  or
  [`prep_wf_data()`](https://n-ducharmebarth-noaa.github.io/opal/reference/prep_wf_data.md)
  to flattened vectors and then to
  [`get_length_like()`](https://n-ducharmebarth-noaa.github.io/opal/reference/get_length_like.md)
  or
  [`get_weight_like()`](https://n-ducharmebarth-noaa.github.io/opal/reference/get_weight_like.md).
  `lf_switch` and `wf_switch` select the likelihood. Changes to one
  stage must be checked against the others.
- **Observations**: mark observed data with `OBS()` so that
  `obj$simulate()` and one-step-ahead residuals remain available.

### Documentation style

Documentation and vignettes should be factual and concise. Avoid
promotional language. State limitations and known differences from other
platforms directly.

## Tests

Tests use `testthat` edition 3.

- Shared fixtures are in `tests/testthat/helper-*.R` (for example
  `opaka_inputs()`, `opaka_obj()`, and `synth_*()` helpers). Do not
  build models at the top level of test files; use the cached helpers.

- Numerical reference outputs are in `tests/testthat/_reference/`.
  **Never create or overwrite reference files from within a test.**

- Simulation self-tests are skipped by default. Run them locally with:

  ``` r

  Sys.setenv(OPAL_SELFTEST = "true", OPAL_SELFTEST_NSIM = "30")
  devtools::test(filter = "selftest")
  ```

  `OPAL_SELFTEST_START` may be `"truth"` (default) or `"default"`, and
  `OPAL_SELFTEST_OUT` saves results to a directory. Self-tests also run
  in `.github/workflows/selftest.yaml`. Calibration limits are
  documented at the top of each self-test file; do not change them
  without rerunning the calibration and recording it there.

## Numerical changes

The golden test (`tests/testthat/test-golden-opaka.R`) compares the
objective, gradient, and every `REPORT()` element of the ’opakapaka
quickstart model against a stored reference, at both the starting values
and the MLE. It also checks that the reference file was produced under
the current `.opal_model_scientific_version`, and that the bundled
example fit is compatible with the current model.
[`read_opal_fit()`](https://n-ducharmebarth-noaa.github.io/opal/reference/opal_fit_io.md)
uses the scientific version to identify fits created under an earlier
model contract.

A numerical change is any change that alters the objective, gradient, or
reported quantities for an existing model configuration. Examples
include changes to dynamics, likelihoods, priors, data preparation,
parameter transformations, or default values.

Any change that alters model numerics follows the regeneration protocol:

1.  Show the golden test failing.
2.  Document the change in `NEWS.md`.
3.  Bump `.opal_model_scientific_version`.
4.  Regenerate `tests/testthat/_reference/opaka-quickstart.rds` with
    `data-raw/make-test-references.R`.
5.  Regenerate `inst/extdata/opaka_quickstart_fit.rds` with
    `data-raw/generate-opaka-quickstart-fit.R`.

All five steps belong in the same pull request. The pull request
description should include the golden test failure from step 1 and
explain why the new values are correct, for example through comparison
with a reference model, an analytical result, or a targeted test of the
changed component.

Additional rules:

- Do not resolve a golden test failure by loosening tolerances.
- If a change is not intended to alter numerics and the golden test
  fails, treat the failure as a bug.
- Adding a new `REPORT()` element does not require regeneration.
  Removing or renaming one does, because the golden test checks that all
  reference elements are still reported.
- For changes that could affect estimation performance, run the
  simulation self-tests and report the results in the pull request.

## NEWS.md

Add entries under the development version heading at the top of
`NEWS.md`, following the existing format: the date, a link to the
commit, and a short description of the change and its consequences. For
numerical changes, state that previously saved fits use the earlier
scientific model contract.

## Vignettes

- Package vignettes in `vignettes/` are built by pkgdown and should run
  quickly. Expensive computations, such as MCMC, should use pre-computed
  results cached in `inst/extdata/`, or be guarded with
  `!identical(Sys.getenv("CI"), "true")`.
- Longer development workflows, such as the WCPO bigeye tuna case study
  and the regression baseline, are in `dev/vignettes/` and are not built
  with the package.

## Continuous integration

GitHub Actions run `R CMD check` on Ubuntu (vignettes are not built
during the check), regenerate documentation, run simulation self-tests,
and build the pkgdown site. The pkgdown site is deployed from `main`.

## Coding agents

Repository conventions for coding agents are in
`.github/copilot-instructions.md`, and apply to any agent used on the
repository. When assigning an issue to GitHub Copilot, select `dev` as
the base branch (the default is `main`). Agent-authored pull requests
are reviewed to the same standard as any other, including the
numerical-change protocol.

## License

opal is released under the GNU General Public License, version 3 or
later. By contributing, you agree that your contributions will be
licensed under the same terms.
