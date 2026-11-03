# Contributing to opal

Thank you for your interest in contributing to opal. This guide describes how to propose changes, set up a development environment, and the checks a pull request must pass. Contributions of all sizes are welcome, including bug reports, documentation fixes, tests, and new model features.

opal is under active development and interfaces may change. If you are unsure whether a change fits the project, open an issue first.

## Where to contribute

| Contribution | Location |
|---|---|
| Package code, tests, vignettes, and function documentation | [N-DucharmeBarth-NOAA/opal](https://github.com/N-DucharmeBarth-NOAA/opal) |
| Project website, WCPFC documents, and case-study write-ups | [N-DucharmeBarth-NOAA/opal-documentation](https://github.com/N-DucharmeBarth-NOAA/opal-documentation) (published at <https://connect.fisheries.noaa.gov/opal/>) |
| Bug reports, feature requests, and design discussion | [opal issues](https://github.com/N-DucharmeBarth-NOAA/opal/issues) |

## Issues

Before opening an issue, search existing issues to avoid duplicates.

A **bug report** should include a minimal reproducible example, the output of `sessionInfo()`, and the opal version or commit. For numerical problems, state what you expected, what you observed, and how you checked it (for example, against an SS3 or MULTIFAN-CL reference model).

A **feature request or substantial change** should be discussed in an issue before work starts. Issues intended for implementation, including those assigned to coding agents, should be self-contained, with:

- the files and functions affected, including exact function signatures;
- acceptance criteria;
- whether the change is expected to alter model numerics (see [Numerical changes](#numerical-changes));
- backward-compatibility implications and expected effects on existing tests;
- notes on AD safety for any code that operates on estimated quantities.

## Branches and pull requests

- **`main`** is the stable release branch.
- **`dev`** is the active development branch.

Create a feature branch from `dev` and open pull requests against `dev`. Keep each pull request focused on one change. A pull request should:

1. describe what changed and why, and link the relevant issue;
2. state whether model numerics change, and if so follow the regeneration protocol below;
3. add or update tests for new or changed behaviour;
4. update roxygen documentation and vignettes when exported function signatures or behaviour change;
5. add a `NEWS.md` entry for user-visible or numerical changes.

## Development environments

Open a GitHub Codespace from **Code → Codespaces** to build the development environment automatically. CI, Codespaces, and the Copilot coding agent use `renv.lock` as the single source of truth for the R version and package records.

All environments run `.github/scripts/setup-r-env.R`, which restores the packages and applies the TMB OSA-header and `compResidual` installation workarounds. To restore locally, first install the R version recorded in `renv.lock` and the required system libraries, then run from the repository root:

``` sh
Rscript .github/scripts/setup-r-env.R
```

In Codespaces, locked GitHub packages are downloaded from public commit archives first, avoiding API rate limits and repository-scoped token restrictions. If that fails, setup tries the original API URL. Package versions and commit SHAs still come from `renv.lock`; do not install newer GitHub versions to work around a download failure.

The Ubuntu policy and Quarto CLI version default are documented in `.github/actions/setup-opal-env/action.yml`; Codespaces must match them. Workflow-specific tools not recorded in the lockfile remain explicit installs after the shared setup.

After intentional dependency updates, run `renv::snapshot()` and commit `renv.lock`. Commit R-version updates in that lockfile too; CI, Codespaces, and Copilot follow it without separate version edits.

Copilot reads its thin setup-workflow wrapper from the default branch, `main`, but uses the working checkout's environment definition. Scheduled runs and manual runs on `main` with `ref: dev` warm caches for dev-based tasks. After merging environment setup changes to `main`, maintainers should run that manual setup, merge `main` into `dev`, and check that a dev-based Copilot task restores the cache.

## Development workflow

- Load the package with `devtools::load_all()` rather than `source()`.
- After editing roxygen comments, run `devtools::document()`. The `document.yaml` workflow also regenerates `man/` and `NAMESPACE` on pushes that touch `R/`, and commits the result, so pull before continuing work on a branch.
- Run `devtools::test()` before opening a pull request. Run `devtools::check()` for changes that affect exports, dependencies, or documentation.
- Prefer small, targeted changes. Refactors that touch several files or array dimensions should be proposed in an issue first, with their effect on test fixtures described.

## Code guidelines

### AD safety

Model code runs on RTMB's automatic-differentiation tape. Any function that may receive AD values must:

- declare the relevant overloads at the top of the function, for example `"[<-" <- ADoverload("[<-")` (and `"c" <- ADoverload("c")` where needed);
- avoid base R coercions that strip the `advector` class, such as `matrix()` or `as.vector()` on AD values. Fix type loss where it occurs, not in the receiving function;
- take derived biological quantities (for example `M_a`, `spawning_potential_a`, `weight_fya`) as explicit arguments rather than reading them from `data`, so gradients propagate when the underlying parameters are estimated.

After a change, confirm that `obj$fn(obj$par)` and `obj$gr(obj$par)` return finite values.

### Exported functions

Exported functions are documented in `man/` and used in vignettes and tests. Changing a signature requires updating the roxygen documentation, `NAMESPACE`, affected vignettes, and tests in the same pull request.

### Model components

- **Selectivity**: changes to selectivity forms must update the `par_sel` dimensions, `get_map()`, and `get_bounds()` together.
- **Composition likelihoods**: data flow from raw data through `prep_lf_data()` or `prep_wf_data()` to flattened vectors and then to `get_length_like()` or `get_weight_like()`. `lf_switch` and `wf_switch` select the likelihood. Changes to one stage must be checked against the others.
- **Observations**: mark observed data with `OBS()` so that `obj$simulate()` and one-step-ahead residuals remain available.

### Documentation style

Documentation and vignettes should be factual and concise. Avoid promotional language. State limitations and known differences from other platforms directly.

## Tests

Tests use `testthat` edition 3.

- Shared fixtures are in `tests/testthat/helper-*.R` (for example `opaka_inputs()`, `opaka_obj()`, and `synth_*()` helpers). Do not build models at the top level of test files; use the cached helpers.
- Numerical reference outputs are in `tests/testthat/_reference/`. **Never create or overwrite reference files from within a test.**
- Simulation self-tests are skipped by default. Run them locally with:

  ``` r
  Sys.setenv(OPAL_SELFTEST = "true", OPAL_SELFTEST_NSIM = "30")
  devtools::test(filter = "selftest")
  ```

  `OPAL_SELFTEST_START` may be `"truth"` (default) or `"default"`, and `OPAL_SELFTEST_OUT` saves results to a directory. Self-tests also run in `.github/workflows/selftest.yaml`. Calibration limits are documented at the top of each self-test file; do not change them without rerunning the calibration and recording it there.

### Performance benchmarks

`.github/workflows/benchmark.yaml` runs independently of the numerical and
simulation tests. Its `memory` job runs on every pull request to `dev` or
`main`, weekly against `dev`, and on manual dispatch. The PR `timing-ab`
job installs the exact PR BASE and HEAD in separate libraries outside the
checkouts and measures them interleaved on the same runner. PR timing uses five
paired rounds of six indicative cases: objectives and gradients for opaka, BET
and IO YFT. Tape construction, reports, repeated fits, trajectories and Hessians
remain in the full weekly suite and the unchanged memory suite. PR timing targets 10–20 minutes including
cached setup. The measurement step has a 20-minute PR timeout; setup time and
uncached dependency restoration can exceed the target, so this is not yet a
hosted runtime guarantee.

Weekly `timing-ab` retains the full 17-case, ten-round suite, comparing two
separately installed copies of the same exact `dev` SHA. Manual dispatch also
compares BASE=HEAD at the dispatched revision; select `timing_profile: smoke`
for the PR profile or `full` (default) for the weekly profile. This checks
same-runner stability, not cumulative timing creep; weekly instruction counts
provide the committed-baseline creep gate. The full memory suite is unchanged,
including opaka/BET fit counts on PRs and IO YFT fit counts weekly/manually.
The `instructions`
job uses Valgrind on weekly and manual runs only. `OPAL_BENCHMARK_FULL_FITS`
controls only the costly `io_yft_fit_evals` case: it is enabled on weekly/manual
runs and disabled on PRs. The opaka and BET fit-count gates still run on every
PR, alongside objective construction and evaluation. All jobs use the shared
Ubuntu 22.04 environment and one OpenMP, OpenBLAS, and MKL thread.

Benchmarks are opt-in locally. From the repository root, install the package
into a separate library, then run memory measurements:

``` sh
mkdir -p benchmark-out/head-lib
R CMD INSTALL --library="$PWD/benchmark-out/head-lib" .
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export OPAL_BENCHMARK=true OPAL_BENCHMARK_ROOT="$PWD"
export OPAL_BENCHMARK_OUT="$PWD/benchmark-out"
export OPAL_BENCHMARK_LIBRARY="$PWD/benchmark-out/head-lib"
export OPAL_BENCHMARK_TOL=0.10
OPAL_BENCHMARK_LAYER=memory Rscript -e \
  'testthat::test_local(filter = "benchmark-models", stop_on_failure = TRUE)'
```

For instruction counts, install Valgrind and run
`OPAL_BENCHMARK_LAYER=instructions Rscript dev/run-benchmark-instructions.R`
with the same environment. To reproduce PR timing, install BASE and HEAD in
separate libraries, set `OPAL_BENCHMARK_BASE_LIBRARY` and
`OPAL_BENCHMARK_HEAD_LIBRARY` to their absolute paths, and run
`OPAL_BENCHMARK_TIMING_PROFILE=smoke OPAL_BENCHMARK_ROUNDS=5 Rscript dev/run-benchmark-ab.R`.
For the full suite use `OPAL_BENCHMARK_TIMING_PROFILE=full` and
`OPAL_BENCHMARK_ROUNDS=10`. Both retain the 10% threshold and fail timing only
when the lower 95% ratio confidence bound exceeds 1.10. Do not reuse a loaded package between
revisions. Workflow artifacts retain measurements and environment metadata;
the job summaries report results even when a regression fails.

The default regression tolerance is `0.10` (10%). Manual dispatch accepts a
`tolerance` fraction and an `update_baseline` boolean. Set `update_baseline` to generate candidate
baselines, or set `OPAL_BENCHMARK_UPDATE=true` locally. Download the artifacts,
review the measured changes and metadata, and explicitly merge the memory
artifact's `baseline.rds` and instruction artifact's `instructions.rds` before
committing. The committed `tests/testthat/_benchmark/baseline.rds` stores all
layers together in its `results` table; `instructions.rds` is an intermediate
artifact, not a second tracked baseline. The workflow never commits or pushes.
Baselines are measurements,
not universal constants: refresh them on the documented Ubuntu environment
after changes to the runner image, R, dependencies, or benchmark definitions.
Do not treat measurements from another machine as calibrated hosted baselines.

With the same environment as above, explicit regeneration commands are:

``` sh
# Regenerate measured baselines (20 timing samples for applicable cases).
OPAL_BENCHMARK_FULL_FITS=true Rscript dev/update-benchmark-baseline.R 20
# Regenerate deterministic replay trajectories after reviewing model changes.
Rscript dev/update-benchmark-trajectories.R opaka bet io_yft
# Produce reviewable instruction measurements without gating against old values.
OPAL_BENCHMARK_UPDATE=true OPAL_BENCHMARK_LAYER=instructions \
  Rscript dev/run-benchmark-instructions.R
# Regenerate memory/timing baselines and merge reviewed instruction measurements.
OPAL_BENCHMARK_FULL_FITS=true \
  Rscript dev/update-benchmark-baseline.R 20 benchmark-out/instructions.rds
```

For a dispatched update, extract the reviewed memory and instruction artifacts
into `benchmark-out/`, then explicitly combine them:

``` sh
Rscript - <<'RS'
baseline <- readRDS("benchmark-out/baseline.rds")
instructions <- readRDS("benchmark-out/instructions.rds")
baseline$results <- rbind(
  subset(baseline$results, metric != "instructions"),
  instructions$results
)
baseline$metadata$instruction_measurement <- instructions$metadata
# Before saving, record reviewed run/attempt/artifact URLs in
# baseline$metadata$hosted_artifacts alongside each layer's provenance.
saveRDS(baseline, "tests/testthat/_benchmark/baseline.rds", version = 3)
RS
```

For an intentional A/B time or memory regression, explain the trade-off in a nonempty
`Performance rationale: ...` line in the PR description. A repository maintainer
with `maintain` or `admin` permission may apply `performance-accepted`.
Label additions/removals and PR description edits retrigger the workflow, so
removing acceptance re-enables the gate. Before checking out PR code, the workflow checks the live PR
label, rationale, latest label actor, and that actor's current repository
permission through the read-only GitHub API. Missing evidence or API failure
keeps the gate enabled. This exception covers the same-runner A/B comparison
only; it does not suppress measurements or waive the committed-baseline memory
gate, fit-count regressions, or instruction regressions. Fork PRs use `pull_request`, not
`pull_request_target`, and receive no elevated token.

Before making benchmark checks required, maintainers must record at least five
hosted BASE=HEAD timing runs and verify the gate detects an injected 15%
slowdown. Calibrate memory and instruction checks on that environment as well.
Branch protection is configured separately, not by this workflow. The
`instructions` job is deliberately skipped on PRs: its skipped status is not
evidence of an instruction measurement. Require the applicable PR checks
(`memory`, `timing-ab`) only after calibration; evaluate `instructions` through
weekly/manual runs rather than claiming it measures every PR.

#### Calibration evidence and pending rollout (2026-11-03)

At SHA `f0037b83831d1428f5158974f5e1674d151d9967`, hosted
[run 41724716265, attempt 2](https://github.com/N-DucharmeBarth-NOAA/opal/actions/runs/41724716265)
completed the [memory job](https://github.com/N-DucharmeBarth-NOAA/opal/actions/runs/41724716265/job/126228220003)
successfully: measurement took 8 minutes, total job approximately 12 minutes.
[Artifact 13275492072](https://github.com/N-DucharmeBarth-NOAA/opal/actions/runs/41724716265/artifacts/13275492072)
contains measurements and comparisons for 17 cases, all passing. Its provenance
records Ubuntu 22.04, AMD EPYC 9V45, R 4.5.2, RTMB 1.9, TMB 1.9.21 and one
OMP/OpenBLAS/MKL thread. These reviewed rows replace their initial sandbox
candidates; IO YFT full-fit and instruction rows still retain explicitly marked
sandbox provenance pending a full hosted calibration.

The original full A/B job was still running after 90 minutes when inspected;
this is neither a pass nor BASE=HEAD evidence. Run 41721585609 was cancelled;
do not count it as a successful noise run. The shorter profile's hosted runtime
and false-failure rate are not yet validated. Existing unit tests demonstrate
memory/count regression rejection and the 5 MiB RSS floor, but are not hosted
injected-regression evidence. No slowdown injection is present in the code.

The updated [run 41739587874](https://github.com/N-DucharmeBarth-NOAA/opal/actions/runs/41739587874)
returned `action_required`, not a hosted smoke-test result. The agent's
available GitHub tools can read runs/logs/artifacts but cannot dispatch or
approve workflows. A maintainer must approve the current PR workflow and use **Actions → benchmarks →
Run workflow**, select this PR branch, and:

1. Run with `update_baseline: true`, `timing_profile: full` to obtain all 18
   memory/timing/count workloads and the four instruction measurements. Review
   and combine artifacts using the commands above; keep exact SHA, run attempt
   and artifact links. Do not present retained sandbox rows as hosted results.
2. At one unchanged SHA, dispatch at least five independent runs with
   `timing_profile: smoke`. Preserve every run, including failures; each has
   five paired rounds but is only **one** independent hosted run. Record
   case-level ratios and 95% intervals from `benchmark-timing-ab/results.csv`.
3. On this same branch, add an explicitly opt-in, temporary computational
   injection and demonstrate approximately 15% slowdown rejection by timing
   and instructions. Preserve failure summaries/artifacts and injection
   settings, then remove the mechanism and revalidate the clean branch.
   Do not use sleep-only injection or `performance-accepted`.
4. Verify clean measurements against the reviewed hosted baseline before
   configuring required PR checks `memory` and `timing-ab` on `dev` and `main`.
   `instructions` remains weekly/manual, never a required skipped PR check.

## Numerical changes

The golden test (`tests/testthat/test-golden-opaka.R`) compares the objective, gradient, and every `REPORT()` element of the 'opakapaka quickstart model against a stored reference, at both the starting values and the MLE. It also checks that the reference file was produced under the current `.opal_model_scientific_version`, and that the bundled example fit is compatible with the current model. `read_opal_fit()` uses the scientific version to identify fits created under an earlier model contract.

A numerical change is any change that alters the objective, gradient, or reported quantities for an existing model configuration. Examples include changes to dynamics, likelihoods, priors, data preparation, parameter transformations, or default values.

Any change that alters model numerics follows the regeneration protocol:

1. Show the golden test failing.
2. Document the change in `NEWS.md`.
3. Bump `.opal_model_scientific_version`.
4. Regenerate `tests/testthat/_reference/opaka-quickstart.rds` with `data-raw/make-test-references.R`.
5. Regenerate `inst/extdata/opaka_quickstart_fit.rds` with `data-raw/generate-opaka-quickstart-fit.R`.

All five steps belong in the same pull request. The pull request description should include the golden test failure from step 1 and explain why the new values are correct, for example through comparison with a reference model, an analytical result, or a targeted test of the changed component.

Additional rules:

- Do not resolve a golden test failure by loosening tolerances.
- If a change is not intended to alter numerics and the golden test fails, treat the failure as a bug.
- Adding a new `REPORT()` element does not require regeneration. Removing or renaming one does, because the golden test checks that all reference elements are still reported.
- For changes that could affect estimation performance, run the simulation self-tests and report the results in the pull request.

## NEWS.md

Add entries under the development version heading at the top of `NEWS.md`, following the existing format: the date, a link to the commit, and a short description of the change and its consequences. For numerical changes, state that previously saved fits use the earlier scientific model contract.

## Vignettes

- Package vignettes in `vignettes/` are built by pkgdown and should run quickly. Expensive computations, such as MCMC, should use pre-computed results cached in `inst/extdata/`, or be guarded with `!identical(Sys.getenv("CI"), "true")`.
- Longer development workflows, such as the WCPO bigeye tuna case study and the regression baseline, are in `dev/vignettes/` and are not built with the package.

## Continuous integration

GitHub Actions run `R CMD check` on Ubuntu (vignettes are not built during the check), regenerate documentation, run simulation self-tests, and build the pkgdown site. The release pkgdown site is deployed from `main`; pushes to `dev` build and
deploy development documentation under `/opal/dev/`. Both builds check local
links and retain the rendered site as a workflow artifact. Development
vignettes under `dev/vignettes/` are run deliberately, separately from the
three package vignettes.

## Coding agents

Repository conventions for coding agents are in `.github/copilot-instructions.md`, and apply to any agent used on the repository. When assigning an issue to GitHub Copilot, select `dev` as the base branch (the default is `main`). Agent-authored pull requests are reviewed to the same standard as any other, including the numerical-change protocol.

## License

opal is released under the GNU General Public License, version 3 or later. By contributing, you agree that your contributions will be licensed under the same terms.
