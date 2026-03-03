# Save & Compare Refactoring Baseline

## Introduction

This vignette saves a complete numerical baseline of the `opal` BET
model before a refactoring change, then compares a post-refactoring
build against it. The workflow has two modes controlled by the `compare`
flag:

- **`compare = FALSE`**: Build the model, evaluate `fn()` and `gr()`,
  save all outputs to `data/opal_baseline.rda`, and update the dataset
  that is exposed via `data(opal_baseline)`.
- **`compare = TRUE`** (default): Load the `opal_baseline` dataset via
  `data(opal_baseline)` (or the local `data/opal_baseline.rda` file),
  rebuild the same inputs, and check numerical equivalence and speedup.

Set the flag here before knitting:

``` r
compare       <- TRUE

# Use here::here() to resolve paths from the project root reliably
baseline_file            <- here::here("data", "opal_baseline.rda")
baseline_data_file       <- here::here("data", "opal_baseline_data.rda")
baseline_parameters_file <- here::here("data", "opal_baseline_parameters.rda")
baseline_map_file        <- here::here("data", "opal_baseline_map.rda")
bench_iters              <- 20
```

## Load inputs

``` r
library(dplyr)
library(RTMB)
library(bench)
library(here)
library(opal)
```

``` r
data(wcpo_bet_data)
data <- wcpo_bet_data
```

### Length composition data

``` r
data(wcpo_bet_lf)

lf_wide <- wcpo_bet_lf %>%
  tidyr::pivot_wider(
    id_cols     = c(fishery, year, month, ts),
    names_from  = bin,
    values_from = value,
    values_fill = 0
  ) %>%
  dplyr::arrange(fishery, ts)

var_adjust_scalars <- 1 / rep(1, data$n_fishery)

data <- prep_lf_data(data, lf_wide,
                     lf_keep_fisheries = c(8, 9),
                     lf_var_adjust     = var_adjust_scalars)
data$lf_switch <- 1L
```

### Weight composition data

``` r
data(wcpo_bet_wf)

data$wt_bin_start <- 1
data$wt_bin_width <- 1
data$n_wt         <- 200L

wf_wide <- wcpo_bet_wf |>
  tidyr::pivot_wider(
    id_cols     = c(fishery, year, month, ts),
    names_from  = bin,
    values_from = value,
    values_fill = 0
  ) |>
  dplyr::arrange(fishery, ts)

data <- prep_wf_data(data, wf_wide,
                     wf_keep_fisheries = c(2),
                     wf_switch = 1L)
```

## Model setup

### Parameters

``` r
data(wcpo_bet_parameters)

parameters <- list(
  log_B0         = 20,
  log_h          = as.numeric(wcpo_bet_parameters$log_h),
  log_sigma_r    = as.numeric(wcpo_bet_parameters$log_sigma_r),
  log_cpue_q     = as.numeric(wcpo_bet_parameters$log_cpue_q),
  cpue_creep     = as.numeric(wcpo_bet_parameters$cpue_creep),
  log_cpue_tau   = as.numeric(wcpo_bet_parameters$log_cpue_tau),
  log_cpue_omega = as.numeric(wcpo_bet_parameters$log_cpue_omega),
  log_lf_tau     = log(rep(0.1, data$n_fishery)),
  log_wf_tau     = rep(0, data$n_fishery),
  log_L1         = as.numeric(wcpo_bet_parameters$log_L1),
  log_L2         = as.numeric(wcpo_bet_parameters$log_L2),
  log_k          = as.numeric(wcpo_bet_parameters$log_k),
  log_CV1        = as.numeric(wcpo_bet_parameters$log_CV1),
  log_CV2        = as.numeric(wcpo_bet_parameters$log_CV2),
  par_sel        = as.matrix(wcpo_bet_parameters$par_sel),
  rdev_y         = as.numeric(wcpo_bet_parameters$rdev_y)
)
```

### Priors

``` r
data$priors <- get_priors(parameters = parameters, data = data)
```

### Parameter map

``` r
map_sel <- matrix(NA, nrow(parameters$par_sel), ncol(parameters$par_sel))

map <- list(
  log_h          = factor(NA),
  log_sigma_r    = factor(NA),
  cpue_creep     = factor(NA),
  log_cpue_tau   = factor(NA),
  log_cpue_omega = factor(NA),
  log_lf_tau     = factor(rep(NA, length(parameters$log_lf_tau))),
  log_wf_tau     = factor(rep(NA, length(parameters$log_wf_tau))),
  log_L1         = factor(NA),
  log_L2         = factor(NA),
  log_k          = factor(NA),
  log_CV1        = factor(NA),
  log_CV2        = factor(NA),
  par_sel        = factor(map_sel)
)
```

### Build the AD object

``` r
t_build <- system.time({
  obj <- MakeADFun(func = cmb(opal_model, data),
                   parameters = parameters, map = map)
})
cat("MakeADFun build time:", round(t_build["elapsed"], 2), "sec\n")
```

    ## MakeADFun build time: 16.9 sec

``` r
cat("Estimated parameters:", length(obj$par), "\n")
```

    ## Estimated parameters: 270

## Evaluate model

Warm up the AD tape, then time `fn()` and `gr()`:

``` r
invisible(obj$fn(obj$par))
invisible(obj$gr(obj$par))
```

    ## outer mgc:  13629.29

``` r
t_fn <- system.time(nll <- obj$fn(obj$par))
t_gr <- system.time(gr  <- obj$gr(obj$par))
```

    ## outer mgc:  13629.29

``` r
rep  <- obj$report()

cat("obj$fn():", round(t_fn["elapsed"], 4), "sec  (NLL =", round(nll, 4), ")\n")
```

    ## obj$fn(): 0 sec  (NLL = 852582 )

``` r
cat("obj$gr():", round(t_gr["elapsed"], 4), "sec  (max|gr| =", round(max(abs(gr)), 6), ")\n")
```

    ## obj$gr(): 0.025 sec  (max|gr| = 13629.29 )

``` r
cat("gr/fn ratio:", round(t_gr["elapsed"] / max(t_fn["elapsed"], 1e-6), 1), "x\n")
```

    ## gr/fn ratio: 25000 x

Stable timing via
[`bench::mark`](https://bench.r-lib.org/reference/mark.html):

``` r
bm <- bench::mark(
  fn = obj$fn(obj$par),
  gr = obj$gr(obj$par),
  iterations = bench_iters,
  check = FALSE
)
```

    ## outer mgc:  13629.29 
    ## outer mgc:  13629.29 
    ## outer mgc:  13629.29 
    ## outer mgc:  13629.29 
    ## outer mgc:  13629.29 
    ## outer mgc:  13629.29 
    ## outer mgc:  13629.29 
    ## outer mgc:  13629.29 
    ## outer mgc:  13629.29 
    ## outer mgc:  13629.29 
    ## outer mgc:  13629.29 
    ## outer mgc:  13629.29 
    ## outer mgc:  13629.29 
    ## outer mgc:  13629.29 
    ## outer mgc:  13629.29 
    ## outer mgc:  13629.29 
    ## outer mgc:  13629.29 
    ## outer mgc:  13629.29 
    ## outer mgc:  13629.29 
    ## outer mgc:  13629.29 
    ## outer mgc:  13629.29

``` r
bm[, c("expression", "min", "median", "itr/sec","total_time")]
```

    ## # A tibble: 2 × 4
    ##   expression      min   median `itr/sec`
    ##   <bch:expr> <bch:tm> <bch:tm>     <dbl>
    ## 1 fn           13.5µs     14µs   61592. 
    ## 2 gr           24.3ms   24.5ms      40.8

## Likelihood component breakdown

``` r
cat("lp_prior:  ", round(rep$lp_prior, 4), "\n")
```

    ## lp_prior:   0

``` r
cat("lp_penalty:", round(rep$lp_penalty, 4), "\n")
```

    ## lp_penalty: 0

``` r
cat("lp_rec:    ", round(rep$lp_rec, 4), "\n")
```

    ## lp_rec:     205.8579

``` r
cat("lp_cpue:   ", round(sum(rep$lp_cpue), 4),
    " (sum of", length(rep$lp_cpue), "obs)\n")
```

    ## lp_cpue:    490.8378  (sum of 268 obs)

``` r
cat("lp_lf:     ", round(sum(rep$lp_lf), 4),
    " (sum of", length(rep$lp_lf), "obs)\n")
```

    ## lp_lf:      133665.5  (sum of 210 obs)

``` r
cat("lp_wf:     ", round(sum(rep$lp_wf), 4),
    " (sum of", length(rep$lp_wf), "obs)\n")
```

    ## lp_wf:      718157.9  (sum of 112 obs)

``` r
cat("Total NLL: ", round(nll, 4), "\n")
```

    ## Total NLL:  852582

## Save or compare

``` r
if (!compare) {
  # ---- Save baseline ----
  baseline <- list(
    timestamp   = Sys.time(),
    description = "Pre-refactoring baseline",
    opal_version = as.character(packageVersion("opal")),
    dimensions  = list(
      n_year    = data$n_year,
      n_age     = data$n_age,
      n_fishery = data$n_fishery,
      n_season  = data$n_season,
      n_len     = data$n_len,
      n_lf      = data$n_lf,
      n_wt      = data$n_wt,
      n_wf      = data$n_wf
    ),
    nll         = nll,
    gradient    = gr,
    max_gr      = max(abs(gr)),
    timing      = list(
      fn_elapsed = t_fn["elapsed"],
      gr_elapsed = t_gr["elapsed"],
      fn_median  = as.numeric(bm$median[1]),
      gr_median  = as.numeric(bm$median[2])
    ),
    report      = list(
      number_ysa         = rep$number_ysa,
      spawning_biomass_y = rep$spawning_biomass_y,
      catch_pred_fya     = rep$catch_pred_fya,
      catch_pred_ysf     = rep$catch_pred_ysf,
      hrate_ysa          = rep$hrate_ysa,
      hrate_ysfa         = rep$hrate_ysfa,
      sel_fya            = rep$sel_fya,
      lp_prior           = rep$lp_prior,
      lp_penalty         = rep$lp_penalty,
      lp_rec             = rep$lp_rec,
      lp_cpue            = rep$lp_cpue,
      lp_lf              = rep$lp_lf,
      lp_wf              = rep$lp_wf
    ),
    par_values  = obj$par
  )

  opal_baseline <- baseline
  save(opal_baseline, file = baseline_file)

  opal_baseline_data <- data
  save(opal_baseline_data, file = baseline_data_file)

  opal_baseline_parameters <- parameters
  save(opal_baseline_parameters, file = baseline_parameters_file)

  opal_baseline_map <- map
  save(opal_baseline_map, file = baseline_map_file)

  cat("Baseline saved to:", baseline_file, "\n")
  cat("  NLL:       ", round(nll, 6), "\n")
  cat("  max|gr|:   ", round(max(abs(gr)), 6), "\n")
  cat("  fn median: ", round(as.numeric(bm$median[1]) * 1000, 1), "ms\n")
  cat("  gr median: ", round(as.numeric(bm$median[2]) * 1000, 1), "ms\n")

} else {
  # ---- Compare against baseline ----
  stopifnot(file.exists(baseline_file))
  data(opal_baseline)
  old <- opal_baseline

  cat("Baseline from:", format(old$timestamp), "\n")
  cat("  Description:", old$description, "\n")
  if (!is.null(old$opal_version)) {
    cat("  opal version:", old$opal_version, "\n")
  }
  cat("\n")

  # --- Numerical equivalence ---
  tol <- 1e-10  # report anything above this

  checks <- list(
    nll                = abs(nll - old$nll),
    number_ysa         = max(abs(rep$number_ysa - old$report$number_ysa)),
    spawning_biomass_y = max(abs(rep$spawning_biomass_y - old$report$spawning_biomass_y)),
    catch_pred_fya     = max(abs(rep$catch_pred_fya - old$report$catch_pred_fya)),
    catch_pred_ysf     = max(abs(rep$catch_pred_ysf - old$report$catch_pred_ysf)),
    hrate_ysa          = max(abs(rep$hrate_ysa - old$report$hrate_ysa)),
    hrate_ysfa         = max(abs(rep$hrate_ysfa - old$report$hrate_ysfa)),
    lp_penalty         = abs(rep$lp_penalty - old$report$lp_penalty),
    lp_rec             = abs(rep$lp_rec - old$report$lp_rec),
    lp_cpue            = max(abs(rep$lp_cpue - old$report$lp_cpue)),
    lp_lf              = max(abs(rep$lp_lf - old$report$lp_lf)),
    lp_wf              = max(abs(rep$lp_wf - old$report$lp_wf)),
    max_gradient       = abs(max(abs(gr)) - old$max_gr)
  )

  cat("--- Numerical equivalence ---\n")
  all_pass <- TRUE
  for (nm in names(checks)) {
    val  <- checks[[nm]]
    pass <- val < tol
    flag <- ifelse(pass, "OK", "FAIL")
    cat(sprintf("  %-22s %12.2e  [%s]\n", nm, val, flag))
    if (!pass) all_pass <- FALSE
  }
  cat("\nOverall:", ifelse(all_pass, "PASS", "FAIL"), "\n\n")

  # --- Timing comparison ---
  cat("--- Timing comparison ---\n")
  cat(sprintf("  fn: %.1f ms -> %.1f ms  (%.1fx)\n",
              old$timing$fn_median * 1000,
              as.numeric(bm$median[1]) * 1000,
              old$timing$fn_median / as.numeric(bm$median[1])))
  cat(sprintf("  gr: %.1f ms -> %.1f ms  (%.1fx)\n",
              old$timing$gr_median * 1000,
              as.numeric(bm$median[2]) * 1000,
              old$timing$gr_median / as.numeric(bm$median[2])))

  # --- Likelihood component comparison ---
  cat("\n--- Likelihood components ---\n")
  cat(sprintf("  %-12s %12s %12s %12s\n", "Component", "Baseline", "Current", "Diff"))
  components <- c("lp_prior", "lp_penalty", "lp_rec")
  for (nm in components) {
    old_val <- old$report[[nm]]
    new_val <- rep[[nm]]
    cat(sprintf("  %-12s %12.4f %12.4f %12.2e\n", nm, old_val, new_val, abs(new_val - old_val)))
  }
  sum_components <- c("lp_cpue", "lp_lf", "lp_wf")
  for (nm in sum_components) {
    old_val <- sum(old$report[[nm]])
    new_val <- sum(rep[[nm]])
    cat(sprintf("  %-12s %12.4f %12.4f %12.2e\n", nm, old_val, new_val, abs(new_val - old_val)))
  }
  cat(sprintf("  %-12s %12.4f %12.4f %12.2e\n", "Total NLL", old$nll, nll, abs(nll - old$nll)))
}
```

    ## Baseline from: 2026-03-03 21:20:32 
    ##   Description: Pre-refactoring baseline 
    ##   opal version: 0.0.3 
    ## 
    ## --- Numerical equivalence ---
    ##   nll                        0.00e+00  [OK]
    ##   number_ysa                 3.73e-09  [FAIL]
    ##   spawning_biomass_y         1.19e-07  [FAIL]
    ##   catch_pred_fya             2.27e-12  [OK]
    ##   catch_pred_ysf             1.46e-11  [OK]
    ##   hrate_ysa                  5.42e-19  [OK]
    ##   hrate_ysfa                 5.42e-19  [OK]
    ##   lp_penalty                 0.00e+00  [OK]
    ##   lp_rec                     0.00e+00  [OK]
    ##   lp_cpue                    6.93e-14  [OK]
    ##   lp_lf                      1.46e-11  [OK]
    ##   lp_wf                      2.91e-11  [OK]
    ##   max_gradient               0.00e+00  [OK]
    ## 
    ## Overall: FAIL 
    ## 
    ## --- Timing comparison ---
    ##   fn: 0.0 ms -> 0.0 ms  (0.7x)
    ##   gr: 16.1 ms -> 24.5 ms  (0.7x)
    ## 
    ## --- Likelihood components ---
    ##   Component        Baseline      Current         Diff
    ##   lp_prior           0.0000       0.0000     0.00e+00
    ##   lp_penalty        -0.0000      -0.0000     0.00e+00
    ##   lp_rec           205.8579     205.8579     0.00e+00
    ##   lp_cpue          490.8378     490.8378     2.27e-13
    ##   lp_lf         133665.5125  133665.5125     2.91e-11
    ##   lp_wf         718157.9275  718157.9275     1.16e-10
    ##   Total NLL     852582.0320  852582.0320     0.00e+00
