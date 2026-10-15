# Seasonal CPUE correction: impact and rebuilding saved fits

Scientific model contract v5 corrects how CPUE observations select
population abundance in models with more than one season. The size of
the effect depends on the data and model. The paired examples below
quantify it for controlled simulations; they are **not estimates of the
effect on a real stock assessment**.

## What changed

CPUE time steps use the flattened convention

``` r

ts <- (year_index - 1L) * n_season + season_index
```

The old calculation used `number_ysa[ts, 1, ]` and selectivity from year
`ts`. The corrected calculation decodes `year_index` and `season_index`,
then uses `number_ysa[year_index, season_index, ]` and selectivity from
that year. For example, with four seasons, `ts = 4` means year 1, season
4; the old code used year 4, season 1. With one season, both
calculations use exactly the same cells.

For a fixed parameter vector, this correction changes CPUE predictions
and the CPUE likelihood. It does not change the underlying population
trajectory or the other likelihood components. Refitting can then change
parameters, population trajectories, and derived quantities. Reusing an
old optimum or posterior does not account for that change.

## Quantified examples

The comparison uses the complete model with 12 years, two fisheries, two
CPUE indices, and multinomial length and weight compositions. It
estimates `log_B0` and the two log catchabilities; other parameters
remain fixed. Annual catches are the same across the one-, two-, and
four-season scenarios, with seasonal shares proportional to
`1:n_season`. A normal prior on `log_B0` is centred on `log(100000)`,
with SD 0.3. CPUE and composition observations are simulated with seed
90210 under the corrected model.

Within each pair, observations, initial parameters, maps, priors,
bounds, and optimiser controls are identical. Only the CPUE indexing
function changes. The legacy function is frozen from `main` at commit
`2dd8e3d` (contract v4), while the rest of the model uses current code.
This isolates the CPUE correction from other changes in the release.
Each fit uses two [`nlminb()`](https://rdrr.io/r/stats/nlminb.html)
passes, with `rel.tol = 1e-10`, `x.tol = 1e-8`, and iteration/evaluation
limits of 10,000.

**Coverage limitation:** to make both versions evaluable, the paired
fits retain `ts = 1:12` for each index: 24 observations in total. These
span all 12 years in the single-season control, the first six years with
two seasons, and the first three years with four seasons. Composition
data cover all 12 years. Therefore, differences between the two- and
four-season examples are not a general estimate of how bias scales with
the number of seasons.

At the same starting parameters, the maximum absolute CPUE prediction
changes were approximately **8.24%** with two seasons and **12.28%**
with four seasons. Percentage changes use the legacy prediction as the
denominator. The objective differences below arise entirely from CPUE.

| Seasons | Median absolute prediction change (%) | Maximum absolute prediction change (%) | NLL change (v5 − v4) |
|---:|---:|---:|---:|
| 1 | 0.000 | 0.000 | 0.000 |
| 2 | 7.389 | 8.243 | -5.656 |
| 4 | 7.222 | 12.280 | -4.481 |

Table 1: Effect of correcting CPUE indexing at identical parameter
values. Negative NLL changes mean a smaller negative log-likelihood
under the corrected calculation.

After refitting, estimated unfished spawning output (`B0`) increased by
approximately **4.99%** and **11.71%**, respectively. Terminal spawning
output increased by **5.32%** and **12.57%**. These are conditional
results from the specified simulated examples, not a universal
correction factor or a bound on the possible effect. Different
observation coverage, recruitment, priors, selectivity, or estimated
parameter sets can give different magnitudes or signs.

| Seasons | Quantity | Legacy indexing | Corrected indexing | Change (%) |
|---:|:---|---:|---:|---:|
| 1 | B0 | 84607.5704 | 84607.5704 | 0.0000 |
| 1 | q1 | 0.7899 | 0.7899 | 0.0000 |
| 1 | q2 | 1.1508 | 1.1508 | 0.0000 |
| 1 | terminal spawning output | 80056.3102 | 80056.3102 | 0.0000 |
| 1 | terminal depletion | 0.9462 | 0.9462 | 0.0000 |
| 2 | B0 | 81518.8183 | 85587.7986 | 4.9915 |
| 2 | q1 | 0.7877 | 0.7899 | 0.2817 |
| 2 | q2 | 1.1475 | 1.1508 | 0.2814 |
| 2 | terminal spawning output | 76470.7126 | 80541.5717 | 5.3234 |
| 2 | terminal depletion | 0.9381 | 0.9410 | 0.3162 |
| 4 | B0 | 78050.5707 | 87187.2527 | 11.7061 |
| 4 | q1 | 0.7872 | 0.7899 | 0.3502 |
| 4 | q2 | 1.1467 | 1.1508 | 0.3496 |
| 4 | terminal spawning output | 72737.2999 | 81878.8540 | 12.5679 |
| 4 | terminal depletion | 0.9319 | 0.9391 | 0.7715 |

Table 2: Paired fitted quantities. Changes are 100 × (v5 / v4 − 1);
depletion changes are relative percentages, not percentage points.
Spawning output uses the simulated model’s units.

All six fits returned optimiser convergence code zero, a maximum
absolute gradient below 0.001, a positive-definite Hessian, and passing
biological checks. In the single-season control, objectives, gradients,
reports, and fitted quantities matched exactly. Existing Opakapaka
golden tests separately check the bundled single-season assessment.

| Seasons | Indexing | Objective | Convergence | Max. gradient | Min. Hessian eigenvalue | Biology passes |
|---:|:---|---:|---:|:---|---:|:---|
| 1 | v4 | 265.3073 | 0 | 3.92e-07 | 10.1096 | TRUE |
| 1 | v5 | 265.3073 | 0 | 3.92e-07 | 10.1096 | TRUE |
| 2 | v4 | 270.5169 | 0 | 3.98e-05 | 10.0509 | TRUE |
| 2 | v5 | 265.0644 | 0 | 2.66e-07 | 10.2878 | TRUE |
| 4 | v4 | 269.6137 | 0 | 6.21e-05 | 9.8717 | TRUE |
| 4 | v5 | 265.4842 | 0 | 1.55e-04 | 10.5249 | TRUE |

Table 3: Numerical and biological checks for each paired fit. Passing
checks do not make the legacy time indexing scientifically correct.

## Full seasonal coverage can fail under the old code

The same experiment also tries observations at every seasonal time step.
When `ts > n_year`, the old calculation attempts to select a nonexistent
year from the selectivity array and raises a subscript error. The
corrected models evaluate successfully. A percentage difference cannot
be assigned to these failed legacy evaluations.

| Seasons | Observations | Legacy error            | Corrected objective finite |
|--------:|-------------:|:------------------------|:---------------------------|
|       2 |           48 | subscript out of bounds | TRUE                       |
|       4 |           96 | subscript out of bounds | TRUE                       |

Table 4: Evaluation with full seasonal CPUE coverage.

## Which bundled assessments are affected?

The three bundled input datasets all use `n_season = 1`. The BET dataset
has quarterly model steps represented along its `n_year` axis, rather
than four seasons within each year. The CPUE indexing correction
therefore makes no numerical change to these datasets under that
convention. Calendar labels such as “quarterly” alone do not identify an
affected model; inspect its dimensions and how `cpue_data$ts` was
constructed.

| Dataset       | Year-axis steps | Seasons | First ts | Last ts |
|:--------------|----------------:|--------:|---------:|--------:|
| opal_baseline |             268 |       1 |        1 |     268 |
| opakapaka     |              75 |       1 |        1 |      75 |
| wcpo_bet      |             268 |       1 |        1 |     268 |

Table 5: Dimensions and CPUE time-step ranges in the bundled input
datasets.

The packaged `opaka_quickstart_fit.rds` and `assessment_example.rds`
already record scientific contract v5. The older `mcmc_fit.rds` and
`opaka_mcmc_fit.rds` are raw `tmbfit` objects without an Opal
scientific-contract record; they cannot be certified as current from
their filenames alone.

## Rebuilding fits from contract v4 or earlier

Opal deliberately rejects old scientific contracts for runtime use,
including single-season fits whose numerical calculations may be
unchanged. Changing a version string or disabling strict checks does not
rebuild or refit a model. Start from its original data, map, parameters,
random-effect specification, bounds, and construction settings, then fit
a new object under v5.

Before rebuilding, verify `cpue_data$ts` against the intended year and
season. If historical inputs used year numbers instead of flattened
seasonal steps, reconstruct `ts` from the original observation dates. Do
not infer missing seasons automatically. This affects the interpretation
of the data, not just the saved-file format.

The following recipe reads a legacy `opal_fit` or a saved `opal_obj`
envelope for its configuration, constructs a new object, refits it, and
writes a new file. It leaves the original file intact and does not copy
old posterior draws or derived results. Inspect the fit diagnostics
before using it.

``` r

saved <- readRDS("assessment-v4.rds")
old <- if (is.list(saved) && identical(saved$format, "opal_obj")) saved$object else saved
stopifnot(inherits(old, "opal_fit") || inherits(old, "opal_obj"))
start <- if (inherits(old, "opal_obj") && !is.null(old$fit$parameters)) {
  old$fit$parameters
} else old$parameters

# Verify old$data$cpue_data$ts before constructing the new object.
rebuilt <- opal_obj(
  data = old$data, parameters = start, map = old$map, random = old$random,
  bounds = old$bounds, control = old$control,
  makeadfun_args = old$makeadfun_args,
  metadata = list(rebuilt_from = "assessment-v4.rds")
)
rebuilt <- opal_fit(rebuilt)
rebuilt <- opal_check(rebuilt, "fit", stop_on_failure = TRUE)
opal_save(rebuilt, "assessment-v5.rds")
```

If posterior inference is needed, rerun
[`opal_mcmc()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_mcmc.md)
against the new target, check the retained chains, and recompute
residuals, summaries, and reference points. This note does not fit or
certify any historical user assessment or its posterior.

The random-effect validation correction in this release is separate: it
uses a numerical Hessian of the marginal objective when RTMB’s exact
Hessian is unavailable. It changes diagnostic handling, not the
objective or estimates. Validation records now use `opal_validation_v3`;
rerun
[`opal_check()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_check.md)
when an existing v5 object’s checks show as stale. Old grid checkpoints
are invalidated so the previous validation failure is not reused.

## Reproduce the comparison

From the repository root, run:

``` sh
Rscript --vanilla dev/seasonal-cpue/compare.R
```

The [comparison
script](https://github.com/N-DucharmeBarth-NOAA/opal/blob/dev/dev/seasonal-cpue/compare.R)
and [frozen legacy
function](https://github.com/N-DucharmeBarth-NOAA/opal/blob/dev/dev/seasonal-cpue/legacy-cpue-v4.R)
write the CSV evidence in `inst/extdata/seasonal-cpue-impact/`.
Predictions at both the starting values and fitted points are retained
in `predictions.csv`. The numerical comparison was run with the
environment recorded below; it was not a fit of a real multi-season
stock assessment.

    Legacy implementation: 2dd8e3de2bbfb4f00a9489c1ae456d8dd2ce5e14
    Current scientific contract: opal_model_contract_v5
    R version 4.6.1 (2026-06-24)
    Platform: aarch64-apple-darwin25.4.0
    Running under: macOS Golden Gate 27.0.1

    Matrix products: default
    BLAS:   /opt/homebrew/Cellar/openblas/0.3.34/lib/libopenblasp-r0.3.34.dylib
    LAPACK: /opt/homebrew/Cellar/r/4.6.1/lib/R/lib/libRlapack.dylib;  LAPACK version 3.12.1

    locale:
    [1] C.UTF-8/C.UTF-8/C.UTF-8/C/C.UTF-8/C.UTF-8

    time zone: Asia/Tokyo
    tzcode source: internal

    attached base packages:
    [1] stats     graphics  grDevices utils     datasets  methods   base

    other attached packages:
    [1] opal_0.1.1     testthat_3.3.2

    loaded via a namespace (and not attached):
     [1] dotCall64_1.2        spam_2.11-4          gtable_0.3.6
     [4] ellipse_0.5.0        tensorA_0.36.2.1     TMB_1.9.25
     [7] ggplot2_4.0.3        lattice_0.23-1       numDeriv_2016.8-1.1
    [10] vctrs_0.7.3          tools_4.6.1          generics_0.1.4
    [13] stats4_4.6.1         parallel_4.6.1       tibble_3.3.1
    [16] pkgconfig_2.0.3      Matrix_1.7-6         checkmate_2.3.4
    [19] RColorBrewer_1.1-3   S7_0.2.2             desc_1.4.3
    [22] ggridges_0.5.7       distributional_0.9.0 lifecycle_1.0.5
    [25] compiler_4.6.1       farver_2.1.2         stringr_1.6.0
    [28] circular_0.5-2       brio_1.1.5           statmod_1.5.2
    [31] mnormt_2.1.2         codetools_0.2-20     movMF_0.2-11
    [34] pillar_1.11.1        RTMB_2.0             MASS_7.3-66
    [37] boot_1.3-32          abind_1.4-8          nlme_3.1-171
    [40] fracdiff_1.5-4       posterior_1.7.0      tidyselect_1.2.1
    [43] mvtnorm_1.4-2        stringi_1.8.9        slam_0.1-56
    [46] dplyr_1.2.1          reshape2_1.4.5       rprojroot_2.1.1
    [49] grid_4.6.1           colorspace_2.1-4     cli_3.6.6
    [52] magrittr_2.0.5       pkgbuild_1.4.8       RTMBdist_1.1.0
    [55] withr_3.0.3          sn_2.1.3             scales_1.4.0
    [58] backports_1.5.1      forecast_9.0.2       sparseinv_0.1.4
    [61] otel_0.2.0           timeDate_4052.112    zoo_1.9-1
    [64] urca_1.3-4           rlang_1.3.0          SparseNUTS_1.0.2
    [67] Rcpp_1.1.2           glue_1.8.1           pkgload_1.5.3
    [70] R6_2.6.1             plyr_1.8.9
