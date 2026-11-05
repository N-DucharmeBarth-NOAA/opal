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

## Development environments

Open a GitHub Codespace from **Code → Codespaces** to build the
development environment automatically. CI, Codespaces, and the Copilot
coding agent use `renv.lock` as the single source of truth for the R
version and package records.

All environments run `.github/scripts/setup-r-env.R`, which restores the
packages and applies the TMB OSA-header and `compResidual` installation
workarounds. To restore locally, first install the R version recorded in
`renv.lock` and the required system libraries, then run from the
repository root:

``` sh
Rscript .github/scripts/setup-r-env.R
```

In Codespaces, locked GitHub packages are downloaded from public commit
archives first, avoiding API rate limits and repository-scoped token
restrictions. If that fails, setup tries the original API URL. Package
versions and commit SHAs still come from `renv.lock`; do not install
newer GitHub versions to work around a download failure.

The Ubuntu policy and Quarto CLI version default are documented in
`.github/actions/setup-opal-env/action.yml`; Codespaces must match them.
Workflow-specific tools not recorded in the lockfile remain explicit
installs after the shared setup.

After intentional dependency updates, run
[`renv::snapshot()`](https://rstudio.github.io/renv/reference/snapshot.html)
and commit `renv.lock`. Commit R-version updates in that lockfile too;
CI, Codespaces, and Copilot follow it without separate version edits.

Copilot reads its thin setup-workflow wrapper from the default branch,
`main`, but uses the working checkout’s environment definition.
Scheduled runs and manual runs on `main` with `ref: dev` warm caches for
dev-based tasks. After merging environment setup changes to `main`,
maintainers should run that manual setup, merge `main` into `dev`, and
check that a dev-based Copilot task restores the cache.

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
  [`get_map()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/get_map.md),
  and
  [`get_bounds()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/get_bounds.md)
  together.
- **Composition likelihoods**: data flow from raw data through
  [`prep_lf_data()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/prep_lf_data.md)
  or
  [`prep_wf_data()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/prep_wf_data.md)
  to flattened vectors and then to
  [`get_length_like()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/get_length_like.md)
  or
  [`get_weight_like()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/get_weight_like.md).
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

### Performance benchmarks

`.github/workflows/benchmark.yaml` runs independently of the numerical
and simulation tests. Its `memory` job runs on every pull request to
`dev` or `main`, weekly against `dev`, and on manual dispatch. The PR
`timing-ab` job installs the exact PR BASE and HEAD in separate
libraries outside the checkouts and measures them interleaved on the
same runner. PR timing uses five paired rounds of four indicative cases:
objectives and gradients for Opakapaka and IO YFT. BET timing, tape
construction, reports, repeated fits, trajectories and Hessians remain
in the full weekly suite and the unchanged memory suite. In hosted
validation the PR timing measurement step took 7–11 minutes, within its
20-minute PR timeout; setup and uncached dependency restoration add to
the job time. Each measured case/round still runs in a fresh subprocess.
A/B reuses one separate, fresh idle reference per installed library
across the job, rather than reloading the same idle package before every
round; this matches how the multi-case memory runner subtracts one idle
reference.

PR timing measures GitHub’s test merge commit against its first parent,
the target-branch tip it merged, so target-branch changes are not
attributed to the PR. Rounds alternate whether BASE or HEAD runs first.

Optimizer-path metrics, meaning fit counts (`fn_evals`, `gr_evals`,
`iterations`) and the fit cases’ `mem_alloc_bytes`, are compared BASE
versus HEAD within `timing-ab`, one fit per revision on the same runner,
not against the committed baseline. They are deterministic on one
machine but differ between hosted CPU models: identical code gave 1225
BET fn evaluations on an AMD EPYC 9V45 and 1651 on an Intel Xeon 8370C,
and allocated 30% more on the latter. PRs gate the well-conditioned
`opaka_fit` (+10% fails; `performance-accepted` does not waive count
failures). The full profile also compares `bet_fit` and
`io_yft_fit_evals`, but only warns, because their paths change with
harmless floating-point rounding. The committed baseline records these
metrics as trends; it still gates the fit cases’ peak RSS and structural
size, which do not depend on the path.

Weekly `timing-ab` retains the full 17-case, ten-round suite, comparing
two separately installed copies of the same exact `dev` SHA. This checks
same-runner stability, not cumulative timing creep; daily instruction
counts provide the committed-baseline creep gate. Manual dispatch also
compares BASE=HEAD at the dispatched revision; select
`timing_profile: smoke` for the PR profile or `full` (default) for the
weekly profile. For BASE-versus-HEAD coverage of BET, tape construction,
trajectories, reports, Hessians and the BET/IO YFT fit counts, dispatch
on the PR branch with `base_ref` set to the target branch or SHA and
`timing_profile: full`; this on-demand comparison is not a required PR
check. Each scheduled or manual run has its own concurrency group, so
queued calibration runs are never cancelled; superseded PR runs are.

The full memory suite runs on PRs and records every workload against the
committed baseline. The `instructions` job uses Valgrind on daily and
manual runs only. `OPAL_BENCHMARK_FULL_FITS` controls only the costly
`io_yft_fit_evals` case: it is enabled on scheduled and manual runs and
disabled on PRs. All jobs use the shared Ubuntu 22.04 environment and
one OpenMP, OpenBLAS, and MKL thread.

Benchmarks are opt-in locally. From the repository root, install the
package into a separate library, then run memory measurements:

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
with the same environment. To reproduce PR timing, install BASE and HEAD
in separate libraries, set `OPAL_BENCHMARK_BASE_LIBRARY` and
`OPAL_BENCHMARK_HEAD_LIBRARY` to their absolute paths, and run
`OPAL_BENCHMARK_TIMING_PROFILE=smoke OPAL_BENCHMARK_ROUNDS=5 Rscript dev/run-benchmark-ab.R`.
For the full suite use `OPAL_BENCHMARK_TIMING_PROFILE=full` and
`OPAL_BENCHMARK_ROUNDS=10`. Both retain the 10% threshold and fail
timing only when the lower 95% ratio confidence bound exceeds 1.10. Do
not reuse a loaded package between revisions. Workflow artifacts retain
measurements and environment metadata; the job summaries report results
even when a regression fails. Re-running a workflow deletes the previous
attempt’s artifacts (its job summary remains), so download them before
re-running.

The default regression tolerance is `0.10` (10%). Manual dispatch
accepts a `tolerance` fraction and an `update_baseline` boolean. Set
`update_baseline` to generate candidate baselines, or set
`OPAL_BENCHMARK_UPDATE=true` locally. Download the artifacts, review the
measured changes and metadata, and explicitly merge the memory
artifact’s `baseline.rds` and instruction artifact’s `instructions.rds`
before committing. The committed
`tests/testthat/_benchmark/baseline.rds` stores all layers together in
its `results` table; `instructions.rds` is an intermediate artifact, not
a second tracked baseline. The workflow never commits or pushes.
Baselines are measurements, not universal constants: refresh them on the
documented Ubuntu environment after changes to the runner image, R,
dependencies, or benchmark definitions. Do not treat measurements from
another machine as calibrated hosted baselines.

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

For a dispatched update, extract the reviewed `benchmark-memory` and
`benchmark-instructions` artifacts, then explicitly combine them,
recording the reviewed run/attempt/artifact URLs:

``` sh
Rscript dev/combine-benchmark-baseline.R \
  benchmark-out/memory/baseline.rds benchmark-out/instructions/instructions.rds \
  --hosted https://github.com/N-DucharmeBarth-NOAA/opal/actions/runs/<run>/attempts/<n> \
  https://github.com/N-DucharmeBarth-NOAA/opal/actions/runs/<run>/artifacts/<id>
```

The script requires all four instruction cases and every memory, timing
and fit-count metric. A partial artifact cannot replace layer-wide
provenance.

For an intentional A/B time or memory regression, explain the trade-off
in a nonempty `Performance rationale: ...` line in the PR description. A
repository maintainer with `maintain` or `admin` permission may apply
`performance-accepted`. Label additions/removals and PR description
edits retrigger the workflow, so removing acceptance re-enables the
gate. Before checking out PR code, the workflow checks the live PR
label, rationale, latest label actor, and that actor’s current
repository permission through the read-only GitHub API. Missing evidence
or API failure keeps the gate enabled. This exception covers the
same-runner A/B comparison only; it does not suppress measurements or
waive the committed-baseline memory gate, fit-count regressions, or
instruction regressions. Fork PRs use `pull_request`, not
`pull_request_target`, and receive no elevated token.

#### Hosted calibration and validation (2026-11-04)

Validation ran on the PR branch through a temporary push trigger, since
`workflow_dispatch` needs the workflow on the default branch. Each
pushed commit selected jobs and settings in its message. A temporary,
opt-in injection mechanism provided the regressions. The trigger and the
injection code were removed afterwards (`498ef62` and the following
commit). All runs used Ubuntu 22.04, R 4.5.2 and one OMP/OpenBLAS/MKL
thread. Hosted runners varied between Intel Xeon 8370C/8573C and AMD
EPYC 7763/9V45/9V74.

| Scenario | Commit | Run | Result |
|----|----|----|----|
| A1 calibration (`update` mode) | `7ea4c1f` | [41895537099](https://github.com/N-DucharmeBarth-NOAA/opal/actions/runs/41895537099) | `memory` and `instructions` succeeded; [memory artifact 13356907664](https://github.com/N-DucharmeBarth-NOAA/opal/actions/runs/41895537099/artifacts/13356907664) (Intel Xeon 8370C), [instruction artifact 13356677878](https://github.com/N-DucharmeBarth-NOAA/opal/actions/runs/41895537099/artifacts/13356677878) (AMD EPYC 7763). Combined into the committed 136-row baseline (`2ec851e`). |
| A2 noise, BASE=HEAD smoke, 5 rounds | `9db1de8` | [41904602801](https://github.com/N-DucharmeBarth-NOAA/opal/actions/runs/41904602801), attempts 1–5 | 5/5 passed; see the table below. |
| B1 +15% computational slowdown on HEAD | `0c22771` | [41922666176](https://github.com/N-DucharmeBarth-NOAA/opal/actions/runs/41922666176) | `instructions` failed on all four cases (ratios 1.150–1.157). `timing-ab` failed on `opaka_fn` and `opaka_gr`; IO YFT was slower but not detected (below). |
| B2 memory: hold max(15% of baseline, 6 MiB) | `4dd74b2` | [41922799259](https://github.com/N-DucharmeBarth-NOAA/opal/actions/runs/41922799259) | `memory` failed only on the four injected cases. |
| B2b perturbed `opaka_fit` start on HEAD | `02b76eb` | [41922920727](https://github.com/N-DucharmeBarth-NOAA/opal/actions/runs/41922920727) | `timing-ab` failed on `opaka_fit`: `fn_evals` 315→917, `gr_evals` 258→814, `iterations` 255→812, `mem_alloc_bytes` 1.84×. Timing passed. |
| B3 clean, injection removed | `ff3e082` | [42069709377](https://github.com/N-DucharmeBarth-NOAA/opal/actions/runs/42069709377) | All three jobs passed: timing ratios 0.98–1.00 (highest upper bound 1.026), `opaka_fit` counts identical, instructions within 0.15%. |

PR-event runs also passed at `40b3db6`
([41891944073](https://github.com/N-DucharmeBarth-NOAA/opal/actions/runs/41891944073),
first verified `pr-merge` comparison), `9db1de8`
([41904613006](https://github.com/N-DucharmeBarth-NOAA/opal/actions/runs/41904613006),
first with the `opaka_fit` count gate) and `ff3e082`
([42069724617](https://github.com/N-DucharmeBarth-NOAA/opal/actions/runs/42069724617)).
PR runs superseded by a later push were cancelled, as designed, and are
not evidence.

A2 timing ratios (HEAD/BASE median, with 95% paired-bootstrap interval):

| Attempt | Runner | `opaka_fn` | `opaka_gr` | `io_yft_fn` | `io_yft_gr` |
|----|----|----|----|----|----|
| 1 | Intel Xeon 8573C | 1.000 (0.997–1.000) | 1.002 (0.993–1.002) | 0.946 (0.906–0.985) | 0.986 (0.976–1.054) |
| 2 | AMD EPYC 9V74 | 0.999 (0.994–1.000) | 0.999 (0.997–1.020) | 0.999 (0.999–1.000) | 1.000 (0.995–1.005) |
| 3 | not recovered | 0.999 (0.995–1.002) | 1.000 (0.992–1.005) | 1.004 (0.998–1.006) | 0.995 (0.990–0.998) |
| 4 | not recovered | 0.916 (0.901–1.088) | 0.997 (0.966–1.103) | 0.999 (0.929–1.094) | 1.005 (0.939–1.092) |
| 5 | AMD EPYC 9V74 | 0.999 (0.999–1.009) | 1.005 (0.996–1.008) | 0.995 (0.873–1.007) | 0.996 (0.993–1.002) |

No attempt came near a false failure (highest lower bound 1.004).
Attempt 4 was noisy, with upper bounds up to 1.103. Re-running deleted
the artifacts of attempts 3 and 4, so their rows come from the retained
job summaries.

B1 timing results on an AMD EPYC 9V45 (20 injected measurements, HEAD
only):

| Case        | Ratio | 95% interval | Gate            |
|-------------|-------|--------------|-----------------|
| `opaka_fn`  | 1.152 | 1.133–1.238  | failed          |
| `opaka_gr`  | 1.143 | 1.116–1.203  | failed          |
| `io_yft_fn` | 1.134 | 1.066–1.251  | passed (warned) |
| `io_yft_gr` | 1.140 | 1.092–1.180  | passed (warned) |

In B2, `mem_alloc_bytes` and `structural_bytes` failed for all four
injected cases. Peak RSS failed for `opaka_fn` (2.12×), `opaka_gr`
(2.13×) and `io_yft_gr` (+42 MB, 1.115×). For `io_yft_fn` it rose 32 MB
(1.086×) and passed. The rise exceeded the 5 MiB floor but not 10% of
the 368 MB baseline.

Unit tests (`test-performance-gates.R`, `test-benchmark-models.R`) cover
the comparison logic: the 10% threshold, the 5 MiB RSS floor,
same-runner fit counts, paired-bootstrap intervals, baseline merging and
A/B provenance. The hosted runs above show the same gates working end to
end on GitHub runners.

#### Known limitations

- **PR timing sensitivity.** Timing fails only when the lower interval
  bound exceeds 1.10. Five rounds reliably catch about 15% on the
  Opakapaka cases, but on the noisier IO YFT cases only slowdowns of
  roughly 20–25% or more. Doubling the rounds would narrow intervals by
  about 30%: still short of reliable 15% detection for IO YFT, and over
  the 20-minute PR step limit. Smaller slowdowns are caught by the
  instruction counts, which run daily (below).
- **Peak RSS on large models.** Peak RSS must exceed both 5 MiB and 10%,
  so on the 368 MB IO YFT cases it misses increases below about 37 MB.
  Allocation and structural size still gate those cases.
- **Runner hardware varies.** Fit counts and fit allocations depend on
  the CPU, so they are compared BASE versus HEAD on one runner, never
  against the committed baseline. Committed timing rows are reference
  values, not gates.
- **The harness is part of the PR.** `timing-ab` runs the benchmark
  harness from the HEAD checkout, so a PR can change its own gate.
  Review changes to `tests/testthat/helper-benchmark.R`,
  `helper-performance.R` and `dev/run-benchmark-*.R` as changes to the
  gate.
- **Untested dispatch path.** `workflow_dispatch`, including `base_ref`,
  cannot run until the workflow reaches the default branch. The same
  jobs ran through the temporary push trigger, but the dispatch inputs
  themselves are not hosted-tested.

#### Required checks (maintainer handoff)

Branch protection is configured in repository settings, not by this
workflow.

- Require `memory` and `timing-ab` on `dev` and `main`.
- Do not require `instructions`. It is skipped on PRs, and a skipped
  status is not evidence of a measurement. It runs daily against `dev`
  (all jobs run on Mondays; the other days run only `instructions`). For
  a PR expected to affect speed, a maintainer can dispatch the workflow
  on the PR branch to get instruction counts before merging.
- `performance-accepted` waives only A/B timing and memory failures. It
  never waives fit-count or instruction regressions; those need a
  reviewed baseline refresh or a fix.
- When a daily instruction run fails, the regression is already on
  `dev`. Bisect the merges since the last passing run.

## Numerical changes

The golden test (`tests/testthat/test-golden-opaka.R`) compares the
objective, gradient, and every `REPORT()` element of the ’opakapaka
quickstart model against a stored reference, at both the starting values
and the MLE. It also checks that the reference file was produced under
the current `.opal_model_scientific_version`, and that the bundled
example fit is compatible with the current model.
[`read_opal_fit()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_fit_io.md)
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
and build the pkgdown site. The release pkgdown site is deployed from
`main`; pushes to `dev` build and deploy development documentation under
`/opal/dev/`. Both builds check local links and retain the rendered site
as a workflow artifact. Development vignettes under `dev/vignettes/` are
run deliberately, separately from the three package vignettes.

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
