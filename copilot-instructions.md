# Copilot instructions for opal

Purpose: help AI coding agents be immediately productive working on the
**opal** R package — the **o**pen **p**opulation **a**ssessment
**l**ibrary for fisheries stock assessment.

## Big picture

`opal` is an R package that implements age- and season-structured
population dynamics models for fisheries stock assessment. It is built
on [RTMB](https://github.com/kaskr/RTMB) for automatic differentiation,
enabling gradient-based optimization (`nlminb`), the Laplace
approximation for random effects, and MCMC sampling (`SparseNUTS`).

Two case studies drive development:

- **’Opakapaka** — the bundled quickstart model. It is the basis for
  `vignettes/`, the bundled example fit, the golden regression test, and
  the simulation self-tests. Start here.
- **WCPO bigeye tuna (BET)** — a high-dimensional quarterly model with
  length and weight compositions. Its workflows live in `dev/vignettes/`
  and are not built with the package.

The design is general and is not specific to either case study.

## Repo layout

| Path | Contents |
|----|----|
| `R/` | All package source code. Core logic lives here. |
| `man/` | roxygen2-generated `.Rd` help files. **Do not edit by hand.** |
| `data/` | Bundled `.rda` datasets (`wcpo_bet_*`, `opaka_*`, and `opal_baseline*`) |
| `inst/extdata/` | Bundled portable fitted-model objects and other package files. |
| `data-raw/` | Scripts used to generate package data and test references. |
| `dev/vignettes/` | Development-only analysis and baseline workflows. |
| `tests/testthat/` | `testthat` edition 3 test suite |
| `tests/testthat/_reference/` | Numerical reference outputs. **Never write these from a test.** |
| `vignettes/` | Quarto vignettes including `quickstart.qmd` and `opakapaka-fit-review.qmd`. |
| `renv/` | `renv` lockfile and library for reproducible dependencies |
| `.github/workflows/` | CI: R CMD check, pkgdown deployment, documentation, and self-tests. |
| `.github/CONTRIBUTING.md` | Contribution guide: branches, PR requirements, numerical-change protocol. |
| `DESCRIPTION`, `NAMESPACE` | Standard R package metadata (roxygen2-managed) |

## Where to start

1.  **README.md** — package overview and installation
2.  **`.github/CONTRIBUTING.md`** — what a pull request must satisfy,
    including the numerical-change protocol
3.  **`vignettes/quickstart.qmd`** — worked Opakapaka tutorial: inputs →
    model fitting → diagnostics → plots
4.  **`R/example-opaka.R`**
    ([`opaka_quickstart_inputs()`](https://n-ducharmebarth-noaa.github.io/opal/reference/opaka_quickstart_inputs.md))
    — the shared quickstart configuration
5.  **`R/model.R`**
    ([`opal_model()`](https://n-ducharmebarth-noaa.github.io/opal/reference/opal_model.md))
    — the central model function that orchestrates all components
6.  **`R/dynamics.R`**
    ([`do_dynamics()`](https://n-ducharmebarth-noaa.github.io/opal/reference/do_dynamics.md))
    — age-season forward population simulation
7.  **`R/likelihoods.R`** — CPUE, length-composition, and
    weight-composition likelihood components

## Architecture

### Model function: `opal_model(parameters, data)`

The core function in `R/model.R`. It is an RTMB-compatible closure
that: 1. Unpacks data + parameters via `getAll(data, parameters)` 2.
Runs modular steps in sequence: growth → PLA → weight-at-age → biology →
selectivity → initial numbers → dynamics → priors → likelihoods 3.
Aggregates NLL:
`nll = lp_prior + lp_penalty + lp_rec + lp_init_rec + sum(lp_cpue) + sum(lp_lf) + sum(lp_wf)`
4. Reports derived quantities via `REPORT()` (spawning biomass,
numbers-at-age, etc.)

Usage:
`RTMB::MakeADFun(func = cmb(opal_model, data), parameters = params, map = map)`
where `cmb(f, d)` creates `function(p) f(p, d)`.

### Key components (each is an exported function in `R/`)

| Function | File | Role |
|----|----|----|
| [`get_growth()`](https://n-ducharmebarth-noaa.github.io/opal/reference/get_growth.md), [`get_sd_at_age()`](https://n-ducharmebarth-noaa.github.io/opal/reference/get_sd_at_age.md) | `growth.R` | Schnute VB growth + SD-at-age |
| [`get_pla()`](https://n-ducharmebarth-noaa.github.io/opal/reference/get_pla.md) | `growth.R` | Probability-of-length-at-age matrix (age-length key) |
| [`resolve_bio_vector()`](https://n-ducharmebarth-noaa.github.io/opal/reference/resolve_bio_vector.md), [`get_maturity_at_age()`](https://n-ducharmebarth-noaa.github.io/opal/reference/get_maturity_at_age.md), [`get_weight_at_length()`](https://n-ducharmebarth-noaa.github.io/opal/reference/get_weight_at_length.md) | `growth.R` | Age- or length-based biological inputs, converted on the AD tape |
| [`get_selectivity()`](https://n-ducharmebarth-noaa.github.io/opal/reference/get_selectivity.md) | `selectivity.R` | Dispatches to [`sel_logistic()`](https://n-ducharmebarth-noaa.github.io/opal/reference/sel_logistic.md), [`sel_double_normal()`](https://n-ducharmebarth-noaa.github.io/opal/reference/sel_double_normal.md), [`sel_double_richards()`](https://n-ducharmebarth-noaa.github.io/opal/reference/sel_double_richards.md), or [`sel_length()`](https://n-ducharmebarth-noaa.github.io/opal/reference/sel_length.md) (length-based → age via PLA) |
| [`convert_ss3_selex_to_rtmb()`](https://n-ducharmebarth-noaa.github.io/opal/reference/convert_ss3_selex_to_rtmb.md), [`convert_rtmb_selex_to_ss3()`](https://n-ducharmebarth-noaa.github.io/opal/reference/convert_rtmb_selex_to_ss3.md) | `selectivity.R` | Translate selectivity parameterizations to and from SS3 |
| [`get_initial_numbers()`](https://n-ducharmebarth-noaa.github.io/opal/reference/get_initial_numbers.md) | `dynamics.R` | Equilibrium N-at-age, R0, BH alpha/beta |
| [`do_dynamics()`](https://n-ducharmebarth-noaa.github.io/opal/reference/do_dynamics.md) | `dynamics.R` | Forward age-season simulation, conditioned on observed catch |
| [`get_harvest_rate()`](https://n-ducharmebarth-noaa.github.io/opal/reference/get_harvest_rate.md) | `dynamics.R` | Per-fishery harvest rates (catch / vulnerable abundance) with `posfun` penalty |
| [`get_recruitment()`](https://n-ducharmebarth-noaa.github.io/opal/reference/get_recruitment.md) | `recruitment.R` | Beverton-Holt SRR with log-normal deviations |
| [`get_cpue_like()`](https://n-ducharmebarth-noaa.github.io/opal/reference/get_cpue_like.md) | `likelihoods.R` | Log-normal CPUE likelihood, `n_index` indices |
| [`get_length_like()`](https://n-ducharmebarth-noaa.github.io/opal/reference/get_length_like.md) | `likelihoods.R` | Multinomial / Dirichlet / Dirichlet-multinomial LF likelihood |
| [`get_weight_like()`](https://n-ducharmebarth-noaa.github.io/opal/reference/get_weight_like.md) | `likelihoods.R` | Same three options for weight compositions (predicted by rebinning from length) |
| [`get_priors()`](https://n-ducharmebarth-noaa.github.io/opal/reference/get_priors.md), [`evaluate_priors()`](https://n-ducharmebarth-noaa.github.io/opal/reference/evaluate_priors.md) | `priors.R` | Prior specification and evaluation (normal, lognormal, beta, t) |
| [`get_parameters()`](https://n-ducharmebarth-noaa.github.io/opal/reference/get_parameters.md) | `parameters.R` | Default parameter list |
| [`get_map()`](https://n-ducharmebarth-noaa.github.io/opal/reference/get_map.md) | `parameters.R` | Default map (which params to fix) |
| [`get_bounds()`](https://n-ducharmebarth-noaa.github.io/opal/reference/get_bounds.md), [`check_bounds()`](https://n-ducharmebarth-noaa.github.io/opal/reference/check_bounds.md) | `parameters.R` | Optimization bounds for `nlminb` |
| [`opal_fit()`](https://n-ducharmebarth-noaa.github.io/opal/reference/opal_fit.md), [`save_opal_fit()`](https://n-ducharmebarth-noaa.github.io/opal/reference/opal_fit_io.md), [`read_opal_fit()`](https://n-ducharmebarth-noaa.github.io/opal/reference/opal_fit_io.md) | `opal-fit.R` | Portable fitted-model, MCMC, and derived-result storage |
| [`get_data()`](https://n-ducharmebarth-noaa.github.io/opal/reference/get_data.md) | `get-data.R` | Legacy data builder (BET-specific) |
| [`prep_lf_data()`](https://n-ducharmebarth-noaa.github.io/opal/reference/prep_lf_data.md) | `prep-lf.R` | Prepare length-frequency data for model |
| [`prep_wf_data()`](https://n-ducharmebarth-noaa.github.io/opal/reference/prep_wf_data.md) | `prep-wf.R` | Prepare weight-frequency data for model |
| [`project_dynamics()`](https://n-ducharmebarth-noaa.github.io/opal/reference/project_dynamics.md), [`project_rec_devs()`](https://n-ducharmebarth-noaa.github.io/opal/reference/project_rec_devs.md), [`project_selectivity()`](https://n-ducharmebarth-noaa.github.io/opal/reference/project_selectivity.md) | `projections.R` | Forward projections from estimated or sampled parameters |
| `plot_*()` | `plots.R` | ggplot2 visualization functions |

`NAMESPACE` is the authoritative list of exports; check it before
assuming a function exists.

### Typical workflow (from `vignettes/quickstart.qmd`)

``` r

library(opal)
inputs <- opaka_quickstart_inputs()
data <- inputs$data
params <- inputs$parameters
map <- inputs$map

# Build AD object
obj <- MakeADFun(func = cmb(opal_model, data), parameters = params, map = map)
bounds <- get_bounds(obj, params)

# Optimize (double run for convergence)
opt <- nlminb(obj$par, obj$fn, obj$gr, lower = bounds$lower, upper = bounds$upper)
opt <- nlminb(opt$par, obj$fn, obj$gr, lower = bounds$lower, upper = bounds$upper)

# Diagnostics
check_estimability(obj)
get_cor_pairs(obj)
sdreport(obj)

# Visualize
plot_cpue(data, obj)
plot_biomass_spawning(list(data), list(obj))

# Store fitted state and posterior output without serializing the RTMB object
fit <- opal_fit(data, obj, opt, bounds = bounds, mcmc = mcmc_fit)
save_opal_fit(fit, "fit.rds")
fit <- read_opal_fit("fit.rds", strict = TRUE)
```

### Portable fitted-model objects

`opal_fit` is the durable boundary for model results. It stores plain-R
data, fitted parameters, the parameter map, optimizer output, normalized
MCMC draws, diagnostics, arbitrary derived results (for example
projections), and provenance. RTMB objectives contain session-specific
environments and external pointers, so they are cached in memory but
never serialized. Use
[`opal_fit_object()`](https://n-ducharmebarth-noaa.github.io/opal/reference/opal_fit_object.md)
or
[`opal_fit_report()`](https://n-ducharmebarth-noaa.github.io/opal/reference/opal_fit_report.md)
to rebuild/access runtime state and
[`update_opal_fit()`](https://n-ducharmebarth-noaa.github.io/opal/reference/update_opal_fit.md)
to attach later MCMC or derived results.

Saved fits record `.opal_model_scientific_version`.
[`read_opal_fit()`](https://n-ducharmebarth-noaa.github.io/opal/reference/opal_fit_io.md)
uses it to identify fits created under an earlier model contract, which
is why a numerical change requires a version bump (see below).

## AD-safe coding patterns

These patterns are **critical** — breaking them will cause silent errors
or tape corruption:

- **`ADoverload` at function top**: Every function touching AD values
  must call `"[<-" <- ADoverload("[<-")` (and sometimes
  `"c" <- ADoverload("c")`) before any array subset-assignment.
- **`getAll()` unpacking**: Both
  [`opal_model()`](https://n-ducharmebarth-noaa.github.io/opal/reference/opal_model.md)
  and
  [`do_dynamics()`](https://n-ducharmebarth-noaa.github.io/opal/reference/do_dynamics.md)
  use `getAll(data, parameters, warn = FALSE)` to unpack list elements
  into local scope.
- **No class-stripping coercions**: Base R functions such as
  [`matrix()`](https://rdrr.io/r/base/matrix.html) and
  [`as.vector()`](https://rdrr.io/r/base/vector.html) strip the
  `advector` class. Fix the coercion where the type is lost, not in the
  receiving function.
- **Pass derived biology explicitly**:
  [`do_dynamics()`](https://n-ducharmebarth-noaa.github.io/opal/reference/do_dynamics.md)
  and the likelihoods take `M_a`, `spawning_potential_a`, and
  `weight_fya` as arguments rather than reading them from `data`, so
  gradients propagate when the underlying parameters are estimated.
- **`REPORT()` / `OBS()`**: RTMB’s mechanisms for reporting derived
  quantities and marking simulation-capable observations.
- **PLA as the bridge**: The probability-of-length-at-age matrix
  connects age-based dynamics to length-based observations and
  selectivity. It sits on the AD tape when growth parameters are
  estimated.
- **Flat vector storage**: Length/weight composition observations are
  stored as flat vectors (`lf_obs_flat`, `lf_obs_ints`, `lf_obs_prop`)
  rather than ragged lists, for RTMB compatibility.

## Parameter structure

Key parameters returned by `get_parameters(data)`:

| Parameter | Description |
|----|----|
| `log_B0` | Log unfished spawning biomass |
| `log_h` | Log steepness (BH SRR) |
| `log_sigma_r` | Log recruitment SD |
| `log_cpue_q`, `cpue_creep`, `log_cpue_tau`, `log_cpue_omega` | CPUE observation model |
| `log_L1`, `log_L2`, `log_k` | Growth (Schnute VB) |
| `log_CV1`, `log_CV2` | Growth variability |
| `par_sel` | Selectivity matrix `[n_fishery, 6]` |
| `log_lf_tau` | LF variance adjustment per fishery |
| `log_wf_tau` | WF variance adjustment per fishery |
| `rdev_y` | Annual recruitment deviations |

All log-transformed for unconstrained optimization.
[`get_map()`](https://n-ducharmebarth-noaa.github.io/opal/reference/get_map.md)
fixes many by default (steepness, sigma_r, growth, CPUE
creep/sigma/omega).

## Data structure

The `data` list contains 25+ elements defining the model dimensions and
observations:

- Dimensions: `n_year`, `n_season`, `n_age`, `n_fishery`, `n_len`,
  `n_wt`, `n_index`
- Observations: `catch_obs_ysf`, `cpue_*`, `lf_*`, `wf_*`
- Biology: `M_a` or M parameters, `maturity_at_length`, `lw_a`, `lw_b`
- Switches: `lf_switch` (1=multinomial, 2=Dirichlet, 3=DM), `wf_switch`
  (0=off, or 1/2/3), `catch_units_f` (1=weight, 2=numbers),
  `removal_switch_f`
- Array naming convention uses dimension suffixes: `_ysf`
  (year-season-fishery), `_fya` (fishery-year-age), `_ysa`
  (year-season-age)

Composition data is prepared via
[`prep_lf_data()`](https://n-ducharmebarth-noaa.github.io/opal/reference/prep_lf_data.md)
and
[`prep_wf_data()`](https://n-ducharmebarth-noaa.github.io/opal/reference/prep_wf_data.md),
which convert wide-format data frames to the flat vector format the
model expects.

## Environment & development

- **R packages** are managed with `renv/`. Run
  [`renv::activate()`](https://rstudio.github.io/renv/reference/activate.html)
  then
  [`renv::restore()`](https://rstudio.github.io/renv/reference/restore.html)
  to set up.
- **Key dependencies**: `RTMB`, `RTMBdist`, `SparseNUTS`, `ggplot2`,
  `dplyr`, `forecast`
- **Documentation**: roxygen2-based. After editing `R/*.R` files,
  regenerate with
  [`devtools::document()`](https://devtools.r-lib.org/reference/document.html).
  The `document.yaml` workflow also regenerates and commits `man/` and
  `NAMESPACE` on pushes touching `R/`, so pull before continuing on a
  branch.
- **Tests**: `testthat` edition 3, run with
  [`devtools::test()`](https://devtools.r-lib.org/reference/test.html).
  Shared fixtures live in `helper-*.R` (`opaka_inputs()`, `opaka_obj()`,
  `synth_*()`); do not build models at the top level of test files.
  `tests/testthat/_reference/` holds numerical reference outputs —
  **never create or overwrite them from a test**; regenerate them only
  through the protocol below. Simulation self-tests are gated by
  `OPAL_SELFTEST=true` (see also `OPAL_SELFTEST_NSIM`,
  `OPAL_SELFTEST_START`, `OPAL_SELFTEST_OUT`) and run in
  `.github/workflows/selftest.yaml`. Calibration limits documented at
  the top of a self-test file must not be changed without rerunning the
  calibration and recording it there.
- **CI**: GitHub Actions run `R CMD check` on Ubuntu only and ignore
  vignettes; they also build pkgdown and regenerate documentation.
- **Branch workflow**: PRs go against `dev`; `main` is the stable
  release branch.

## Numerical changes

A change is numerical if it alters the objective, gradient, or any
`REPORT()` element for an existing model configuration — for example
changes to dynamics, likelihoods, priors, data preparation, parameter
transformations, or defaults. `tests/testthat/test-golden-opaka.R`
compares all three against a stored reference at the starting values and
the MLE, and will fail when this happens.

Do **not** loosen tolerances, edit a reference file directly, or
regenerate references from inside a test. If the change was not intended
to alter numerics, treat the failure as a bug and find the cause.

If the change is intended, follow the regeneration protocol, all within
one PR:

1.  Show the golden test failing.
2.  Document the change in `NEWS.md`.
3.  Bump `.opal_model_scientific_version`.
4.  Regenerate `tests/testthat/_reference/opaka-quickstart.rds` with
    `data-raw/make-test-references.R`.
5.  Regenerate `inst/extdata/opaka_quickstart_fit.rds` with
    `data-raw/generate-opaka-quickstart-fit.R`.

The PR description must include the step 1 failure and explain why the
new values are correct (comparison with a reference model, an analytical
result, or a targeted test of the changed component). Adding a new
`REPORT()` element does not require regeneration; removing or renaming
one does. See `.github/CONTRIBUTING.md` for the full contribution guide.

## Debugging tips

- Run `check_estimability(obj)` after fitting to detect non-identifiable
  parameters via Hessian eigenvalue analysis.
- Run `get_cor_pairs(obj, threshold = 0.95)` to find highly correlated
  parameter pairs.
- Use
  [`get_par_table()`](https://n-ducharmebarth-noaa.github.io/opal/reference/get_par_table.md)
  to review initial vs estimated values, gradients, and bounds
  proximity.
- Use `obj$simulate()` with RTMB’s `OBS()` mechanism for
  simulation-based diagnostics.

## AI edit guidance

- **AD safety first**: When adding or modifying any function that
  receives AD values, include `ADoverload("[<-")` at the top. Test that
  `obj$fn(obj$par)` returns a finite value after changes.
- **Preserve exported API signatures**: All exported functions are
  documented in `man/` and used in vignettes. Changing signatures
  requires updating roxygen docs, NAMESPACE, and vignettes.
- **Run tests after changes**:
  [`devtools::test()`](https://devtools.r-lib.org/reference/test.html)
  covers the core components. If adding new functionality, add
  corresponding tests in `tests/testthat/`.
- **Regenerate docs**: After editing roxygen comments in `R/` files, run
  [`devtools::document()`](https://devtools.r-lib.org/reference/document.html)
  to update `man/` and `NAMESPACE`.
- **Keep patches minimal**: Modify one `R/` file at a time, run
  [`devtools::test()`](https://devtools.r-lib.org/reference/test.html)
  and
  [`devtools::check()`](https://devtools.r-lib.org/reference/check.html)
  to validate.
- **Composition data flow**: When modifying likelihood functions,
  understand the full pipeline: raw data →
  [`prep_lf_data()`](https://n-ducharmebarth-noaa.github.io/opal/reference/prep_lf_data.md)/[`prep_wf_data()`](https://n-ducharmebarth-noaa.github.io/opal/reference/prep_wf_data.md)
  → flat vectors → likelihood function. The `lf_switch`/`wf_switch`
  controls which distribution is used.
- **Selectivity changes**: When modifying selectivity, update `par_sel`
  dimensions,
  [`get_map()`](https://n-ducharmebarth-noaa.github.io/opal/reference/get_map.md)
  (which elements are fixed), and
  [`get_bounds()`](https://n-ducharmebarth-noaa.github.io/opal/reference/get_bounds.md)
  simultaneously.
- **Flag numerical impact**: State in the PR description whether a
  change alters model numerics, and follow the protocol above if it
  does.
- **Documentation tone**: Keep documentation and vignettes factual and
  concise. No promotional language; state limitations and known
  differences from other platforms directly.
- **Version control**: Commit changes with clear messages, referencing
  related issues or PRs. Ensure that `.opal_model_scientific_version` is
  updated appropriately to reflect the changes made.
