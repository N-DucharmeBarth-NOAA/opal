# Copilot instructions for opal

Purpose: help AI coding agents be immediately productive working on the
**opal** R package — the **o**pen **p**opulation **a**ssessment
**l**ibrary for fisheries stock assessment.

## Big picture

`opal` is an R package that implements age- and season-structured
population dynamics models for fisheries stock assessment. It is built
on [RTMB](https://github.com/kaskr/RTMB) for automatic differentiation,
enabling gradient-based optimization (`nlminb`) and MCMC sampling
(`SparseNUTS`). The primary application is WCPO bigeye tuna (BET), but
the design is general.

## Repo layout

| Path                       | Contents                                                                                       |
|----------------------------|------------------------------------------------------------------------------------------------|
| `R/`                       | All package source code (23 files). Core logic lives here.                                     |
| `man/`                     | roxygen2-generated `.Rd` help files. **Do not edit by hand.**                                  |
| `data/`                    | Bundled `.rda` datasets (`wcpo_bet_data`, `wcpo_bet_lf`, `wcpo_bet_wf`, `wcpo_bet_parameters`) |
| `tests/testthat/`          | `testthat` edition 3 test suite                                                                |
| `vignettes/`               | `quickstart.Rmd` (overview) and `bet.Rmd` (full worked example)                                |
| `renv/`                    | `renv` lockfile and library for reproducible dependencies                                      |
| `.github/workflows/`       | CI: R-CMD-check, pkgdown deploy, roxygen2, Rmd rendering                                       |
| `DESCRIPTION`, `NAMESPACE` | Standard R package metadata (roxygen2-managed)                                                 |

## Where to start

1.  **README.md** — package overview and installation
2.  **`vignettes/bet.Rmd`** — full worked tutorial: data prep → model
    fitting → diagnostics → plots
3.  **`R/model.R`**
    ([`opal_model()`](https://n-ducharmebarth-noaa.github.io/opal/reference/opal_model.md))
    — the central model function that orchestrates all components
4.  **`R/dynamics.R`**
    ([`do_dynamics()`](https://n-ducharmebarth-noaa.github.io/opal/reference/do_dynamics.md))
    — age-season forward population simulation
5.  **`R/likelihoods.R`** — CPUE, length-composition, and
    weight-composition likelihood components

## Architecture

### Model function: `opal_model(parameters, data)`

The core function in `R/model.R`. It is an RTMB-compatible closure
that: 1. Unpacks data + parameters via `getAll(data, parameters)` 2.
Runs modular steps in sequence: growth → PLA → weight-at-age → biology →
selectivity → initial numbers → dynamics → priors → likelihoods 3.
Aggregates NLL:
`nll = lp_prior + lp_rec + sum(lp_cpue) + sum(lp_lf) + sum(lp_wf)` 4.
Reports derived quantities via `REPORT()` (spawning biomass,
numbers-at-age, etc.)

Usage:
`RTMB::MakeADFun(func = cmb(opal_model, data), parameters = params, map = map)`
where `cmb(f, d)` creates `function(p) f(p, d)`.

### Key components (each is an exported function in `R/`)

| Function                                                                                                                                                                           | File                  | Role                                                                    |
|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-----------------------|-------------------------------------------------------------------------|
| [`get_growth()`](https://n-ducharmebarth-noaa.github.io/opal/reference/get_growth.md), [`get_sd_at_age()`](https://n-ducharmebarth-noaa.github.io/opal/reference/get_sd_at_age.md) | `growth.R`            | Schnute VB growth + SD-at-age                                           |
| [`get_pla()`](https://n-ducharmebarth-noaa.github.io/opal/reference/get_pla.md)                                                                                                    | `growth.R`            | Probability-of-length-at-age matrix (age-length key)                    |
| [`get_selectivity()`](https://n-ducharmebarth-noaa.github.io/opal/reference/get_selectivity.md)                                                                                    | `selectivity.R`       | Logistic or double-normal selectivity (length-based → age via PLA)      |
| [`get_initial_numbers()`](https://n-ducharmebarth-noaa.github.io/opal/reference/get_initial_numbers.md)                                                                            | `dynamics.R`          | Equilibrium N-at-age, R0, BH alpha/beta                                 |
| [`do_dynamics()`](https://n-ducharmebarth-noaa.github.io/opal/reference/do_dynamics.md)                                                                                            | `dynamics.R`          | Forward age-season simulation with Baranov catch equation               |
| [`get_harvest_rate()`](https://n-ducharmebarth-noaa.github.io/opal/reference/get_harvest_rate.md)                                                                                  | `dynamics.R`          | Per-fishery harvest rates with `posfun` penalty                         |
| [`get_recruitment()`](https://n-ducharmebarth-noaa.github.io/opal/reference/get_recruitment.md)                                                                                    | `recruitment.R`       | Beverton-Holt SRR with log-normal deviations                            |
| [`get_cpue_like()`](https://n-ducharmebarth-noaa.github.io/opal/reference/get_cpue_like.md)                                                                                        | `likelihoods.R`       | Log-normal CPUE likelihood                                              |
| [`get_length_like()`](https://n-ducharmebarth-noaa.github.io/opal/reference/get_length_like.md)                                                                                    | `likelihoods.R`       | Multinomial / Dirichlet / Dirichlet-multinomial LF likelihood           |
| [`get_weight_like()`](https://n-ducharmebarth-noaa.github.io/opal/reference/get_weight_like.md)                                                                                    | `likelihoods.R`       | Same three options for weight compositions                              |
| [`get_M()`](https://n-ducharmebarth-noaa.github.io/opal/reference/get_M.md), [`get_M_length()`](https://n-ducharmebarth-noaa.github.io/opal/reference/get_M_length.md)             | `natural-mortality.R` | Piecewise or Lorenzen natural mortality                                 |
| [`evaluate_priors()`](https://n-ducharmebarth-noaa.github.io/opal/reference/evaluate_priors.md)                                                                                    | `priors.R`            | Prior evaluation (normal, lognormal, beta, t)                           |
| [`get_parameters()`](https://n-ducharmebarth-noaa.github.io/opal/reference/get_parameters.md)                                                                                      | `parameters.R`        | Default parameter list                                                  |
| [`get_map()`](https://n-ducharmebarth-noaa.github.io/opal/reference/get_map.md)                                                                                                    | `parameters.R`        | Default map (which params to fix)                                       |
| [`get_bounds()`](https://n-ducharmebarth-noaa.github.io/opal/reference/get_bounds.md)                                                                                              | `parameters.R`        | Optimization bounds for `nlminb`                                        |
| [`get_data()`](https://n-ducharmebarth-noaa.github.io/opal/reference/get_data.md)                                                                                                  | `get-data.R`          | Legacy data builder (BET-specific)                                      |
| [`prep_lf_data()`](https://n-ducharmebarth-noaa.github.io/opal/reference/prep_lf_data.md)                                                                                          | `prep-lf.R`           | Prepare length-frequency data for model                                 |
| [`prep_wf_data()`](https://n-ducharmebarth-noaa.github.io/opal/reference/prep_wf_data.md)                                                                                          | `prep-wf.R`           | Prepare weight-frequency data for model                                 |
| [`run_grid()`](https://n-ducharmebarth-noaa.github.io/opal/reference/run_grid.md)                                                                                                  | `grid.R`              | Grid-based sensitivity analysis (parallel optimization)                 |
| [`sample_grid()`](https://n-ducharmebarth-noaa.github.io/opal/reference/sample_grid.md)                                                                                            | `grid.R`              | Sample grid cells proportional to likelihood for integrated uncertainty |
| [`get_posterior()`](https://n-ducharmebarth-noaa.github.io/opal/reference/get_posterior.md)                                                                                        | `get-posterior.R`     | Extract derived quantities from MCMC draws                              |
| [`project_dynamics()`](https://n-ducharmebarth-noaa.github.io/opal/reference/project_dynamics.md)                                                                                  | `projections.R`       | Forward projections from posterior                                      |
| [`opalprofile()`](https://n-ducharmebarth-noaa.github.io/opal/reference/opalprofile.md)                                                                                            | `profile.R`           | 1D likelihood profiling                                                 |
| `plot_*()`                                                                                                                                                                         | `plots.R`             | ggplot2 visualization functions                                         |

### Typical workflow (from `vignettes/bet.Rmd`)

``` r
library(opal)
data    <- wcpo_bet_data
lf      <- wcpo_bet_lf
params  <- wcpo_bet_parameters

# Prepare composition data
data <- prep_lf_data(data, lf_wide = pivot_wider(...))

# Configure
params  <- get_parameters(data)
priors  <- get_priors(params, data)
map     <- get_map(params)
bounds  <- get_bounds(obj, params)

# Build AD object
obj <- MakeADFun(func = cmb(opal_model, data), parameters = params, map = map)

# Optimize (double run for convergence)
opt <- nlminb(obj$par, obj$fn, obj$gr, lower = bounds$lower, upper = bounds$upper)
opt <- nlminb(opt$par, obj$fn, obj$gr, lower = bounds$lower, upper = bounds$upper)

# Diagnostics
check_estimability(obj)
get_cor_pairs(obj)
sdreport(obj)

# Visualize
rep <- obj$report()
plot_lf(data, rep)
plot_cpue(data, rep)
plot_biomass_spawning(data, rep)
```

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

| Parameter                                                    | Description                         |
|--------------------------------------------------------------|-------------------------------------|
| `log_B0`                                                     | Log unfished spawning biomass       |
| `log_h`                                                      | Log steepness (BH SRR)              |
| `log_sigma_r`                                                | Log recruitment SD                  |
| `log_cpue_q`, `cpue_creep`, `log_cpue_tau`, `log_cpue_omega` | CPUE observation model              |
| `log_L1`, `log_L2`, `log_k`                                  | Growth (Schnute VB)                 |
| `log_CV1`, `log_CV2`                                         | Growth variability                  |
| `par_sel`                                                    | Selectivity matrix `[n_fishery, 6]` |
| `log_lf_tau`                                                 | LF variance adjustment per fishery  |
| `log_wf_tau`                                                 | WF variance adjustment per fishery  |
| `rdev_y`                                                     | Annual recruitment deviations       |

All log-transformed for unconstrained optimization.
[`get_map()`](https://n-ducharmebarth-noaa.github.io/opal/reference/get_map.md)
fixes many by default (steepness, sigma_r, growth, CPUE
creep/sigma/omega).

## Data structure

The `data` list contains 25+ elements defining the model dimensions and
observations:

- Dimensions: `n_year`, `n_season`, `n_age`, `n_fishery`, `n_len`
- Observations: `catch_obs_ysf`, `cpue_*`, `lf_*`, `wf_*`
- Biology: `M_a` or M parameters, `maturity_at_length`, `lw_a`, `lw_b`
- Switches: `lf_switch` (1=multinomial, 2=Dirichlet, 3=DM), `wf_switch`
  (0=off, or 1/2/3)
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
  `dplyr`, `tidyr`, `foreach`, `doParallel`, `loo`, `rstan`, `forecast`,
  `mgcv`
- **Documentation**: roxygen2-based. After editing `R/*.R` files,
  regenerate with `devtools::document()`.
- **Tests**: `testthat` edition 3. Run with `devtools::test()`. Tests
  cover dynamics, growth, all three likelihood types, selectivity, data
  prep, rebinning, and utilities.
- **CI**: GitHub Actions run `R CMD check`, pkgdown site builds,
  roxygen2 re-generation, and Rmd rendering.
- **Branch workflow**: PRs go against `dev`; `main` is the stable
  release branch.

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
- If optimization fails, try
  [`run_grid()`](https://n-ducharmebarth-noaa.github.io/opal/reference/run_grid.md)
  which does triple `nlminb` restarts for robustness.

## AI edit guidance

- **AD safety first**: When adding or modifying any function that
  receives AD values, include `ADoverload("[<-")` at the top. Test that
  `obj$fn(obj$par)` returns a finite value after changes.
- **Preserve exported API signatures**: All exported functions are
  documented in `man/` and used in vignettes. Changing signatures
  requires updating roxygen docs, NAMESPACE, and vignettes.
- **Run tests after changes**: `devtools::test()` covers the core
  components. If adding new functionality, add corresponding tests in
  `tests/testthat/`.
- **Regenerate docs**: After editing roxygen comments in `R/` files, run
  `devtools::document()` to update `man/` and `NAMESPACE`.
- **Keep patches minimal**: Modify one `R/` file at a time, run
  `devtools::test()` and `devtools::check()` to validate.
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
- **Grid/sensitivity**:
  [`get_grid()`](https://n-ducharmebarth-noaa.github.io/opal/reference/get_grid.md)
  creates parameter combinations;
  [`run_grid()`](https://n-ducharmebarth-noaa.github.io/opal/reference/run_grid.md)
  optimizes each. Changes to the parameter or data structure must be
  reflected in both.
