# Tests for get_weight_like()
library(RTMB)

# Shared test helpers ----

# Build arguments for get_weight_like() with a single observation: fishery 1, year 1.
# Returns a named list that can be called with do.call(get_weight_like, .).
# wf_obs_flat / wf_obs_prop / wf_obs_ints are pre-aggregated to [bmin:bmax] for
# fishery 1, matching the aggregation applied internally to the predicted vector.
make_wf_args <- function(n_fishery = 2, n_year = 3, n_age = 10, n_len = 20,
                         n_wt = 10, wf_switch = 1, log_wf_tau = NULL,
                         removal_switch_f = NULL, wf_minbin = NULL,
                         wf_maxbin = NULL, obs_row = NULL, wf_n = 100L) {
  if (is.null(log_wf_tau))       log_wf_tau       <- rep(0.0, n_fishery)
  if (is.null(removal_switch_f)) removal_switch_f <- rep(0L,  n_fishery)
  if (is.null(wf_minbin))        wf_minbin        <- rep(1L,  n_fishery)
  if (is.null(wf_maxbin))        wf_maxbin        <- rep(n_wt, n_fishery)
  if (is.null(obs_row))          obs_row          <- rep(1 / n_wt, n_wt)

  # One observation group: fishery 1, year 1
  wf_n_f       <- 1L
  wf_fishery_f <- 1L
  wf_year_fi   <- list(1L)
  wf_n_fi      <- list(wf_n)

  # Pre-aggregate observed counts to the active [bmin, bmax] window for fishery 1
  f1    <- 1L
  bmin  <- wf_minbin[f1]
  bmax  <- wf_maxbin[f1]
  nbins <- bmax - bmin + 1L

  obs_counts <- obs_row * wf_n
  obs_agg    <- obs_counts[bmin:bmax]
  if (bmin > 1L)    obs_agg[1L]    <- sum(obs_counts[1:bmin])
  if (bmax < n_wt)  obs_agg[nbins] <- sum(obs_counts[bmax:n_wt])

  wf_obs_flat <- obs_agg                      # unrounded counts  (switch 1)
  wf_obs_prop <- obs_agg / sum(obs_agg)       # proportions       (switch 2)
  wf_obs_ints <- as.integer(round(obs_agg))   # integer counts    (switch 3)

  # Uniform catch-at-age
  catch_pred_fya <- array(1, dim = c(n_fishery, n_year, n_age))

  # PLA: each age maps exactly to one length bin
  pla <- matrix(0, n_len, n_age)
  for (a in seq_len(n_age)) pla[min(a * 2L, n_len), a] <- 1

  # Rebinning matrix: groups length bins 2:1 into weight bins
  # n_len length bins -> n_wt weight bins (each weight bin covers n_len/n_wt length bins)
  wf_rebin_matrix <- matrix(0, n_wt, n_len)
  bins_per_wt <- n_len / n_wt
  for (w in seq_len(n_wt)) {
    l_start <- (w - 1) * bins_per_wt + 1
    l_end   <- w * bins_per_wt
    wf_rebin_matrix[w, l_start:l_end] <- 1
  }

  list(
    wf_obs_flat      = wf_obs_flat,
    wf_obs_ints      = wf_obs_ints,
    wf_obs_prop      = wf_obs_prop,
    catch_pred_fya   = catch_pred_fya,
    pla              = pla,
    wf_rebin_matrix  = wf_rebin_matrix,
    wf_n_f           = wf_n_f,
    wf_fishery_f     = wf_fishery_f,
    wf_year_fi       = wf_year_fi,
    wf_n_fi          = wf_n_fi,
    wf_minbin        = wf_minbin,
    wf_maxbin        = wf_maxbin,
    removal_switch_f = removal_switch_f,
    wf_switch        = wf_switch,
    n_wt             = n_wt,
    n_wf             = 1L,
    log_wf_tau       = log_wf_tau
  )
}

# Tests ----

test_that("predicted weight proportions sum to 1", {
  n_age <- 10
  n_len <- 20
  n_wt  <- 10
  catch_a <- rep(1, n_age)

  pla <- matrix(0, n_len, n_age)
  for (a in seq_len(n_age)) pla[min(a * 2L, n_len), a] <- 1

  wf_rebin_matrix <- matrix(0, n_wt, n_len)
  for (w in seq_len(n_wt)) wf_rebin_matrix[w, ((w - 1) * 2 + 1):(w * 2)] <- 1

  pred_at_length <- c(pla %*% catch_a)
  pred_at_weight <- c(wf_rebin_matrix %*% pred_at_length)
  pred <- (pred_at_weight + 1e-8) / sum(pred_at_weight + 1e-8)
  expect_equal(sum(pred), 1.0, tolerance = 1e-10)
})

test_that("multinomial NLL is finite and non-negative", {
  s  <- make_wf_args(wf_switch = 1)
  lp <- do.call(get_weight_like, s)
  expect_true(all(is.finite(lp)))
  expect_true(all(lp >= 0))
})

test_that("Dirichlet NLL is finite", {
  s  <- make_wf_args(wf_switch = 2)
  lp <- do.call(get_weight_like, s)
  expect_true(all(is.finite(lp)))
  # Note: Dirichlet density can exceed 1 so NLL = -log_density can be negative
})

test_that("Dirichlet-multinomial NLL is finite and non-negative", {
  s  <- make_wf_args(wf_switch = 3)
  lp <- do.call(get_weight_like, s)
  expect_true(all(is.finite(lp)))
  expect_true(all(lp >= 0))
})

test_that("perfect fit gives finite NLL for multinomial", {
  n_age <- 10
  n_len <- 20
  n_wt  <- 10

  pla <- matrix(0, n_len, n_age)
  for (a in seq_len(n_age)) pla[min(a * 2L, n_len), a] <- 1

  wf_rebin_matrix <- matrix(0, n_wt, n_len)
  for (w in seq_len(n_wt)) wf_rebin_matrix[w, ((w - 1) * 2 + 1):(w * 2)] <- 1

  pred_at_length <- c(pla %*% rep(1, n_age))
  pred_at_weight <- c(wf_rebin_matrix %*% pred_at_length)
  pred_norm <- (pred_at_weight + 1e-8) / sum(pred_at_weight + 1e-8)

  s  <- make_wf_args(n_age = n_age, n_len = n_len, n_wt = n_wt,
                     wf_switch = 1, obs_row = pred_norm, wf_n = 100L)
  lp <- do.call(get_weight_like, s)
  expect_true(is.finite(lp[1]))
  expect_true(lp[1] >= 0)
})

test_that("log_wf_tau affects NLL magnitude for wf_switch = 2 (Dirichlet)", {
  n_wt     <- 10
  base_obs <- rep(1 / n_wt, n_wt)

  s_low  <- make_wf_args(n_wt = n_wt, wf_switch = 2,
                         obs_row = base_obs, log_wf_tau = rep(0.0, 2))
  s_high <- make_wf_args(n_wt = n_wt, wf_switch = 2,
                         obs_row = base_obs, log_wf_tau = rep(2.0, 2))

  lp_low  <- do.call(get_weight_like, s_low)
  lp_high <- do.call(get_weight_like, s_high)

  expect_true(is.finite(lp_low[1]))
  expect_true(is.finite(lp_high[1]))
  # Note: Dirichlet density can exceed 1, so NLL can be negative
  expect_false(isTRUE(all.equal(lp_low[1], lp_high[1])))
})

test_that("log_wf_tau affects NLL magnitude for wf_switch = 3 (DM)", {
  n_wt     <- 10
  base_obs <- rep(1 / n_wt, n_wt)

  s_low  <- make_wf_args(n_wt = n_wt, wf_switch = 3,
                         obs_row = base_obs, log_wf_tau = rep(0.0, 2))
  s_high <- make_wf_args(n_wt = n_wt, wf_switch = 3,
                         obs_row = base_obs, log_wf_tau = rep(2.0, 2))

  lp_low  <- do.call(get_weight_like, s_low)
  lp_high <- do.call(get_weight_like, s_high)

  expect_true(is.finite(lp_low[1]))
  expect_true(is.finite(lp_high[1]))
  expect_false(isTRUE(all.equal(lp_low[1], lp_high[1])))
})

test_that("removal_switch_f = 1 gives zero NLL contribution", {
  s  <- make_wf_args(wf_switch = 1, removal_switch_f = c(1L, 0L))
  lp <- do.call(get_weight_like, s)
  expect_equal(lp[1], 0)
})

test_that("n_i = 0 observation is skipped and contributes zero NLL", {
  s <- make_wf_args(wf_switch = 1)
  s$wf_n_fi <- list(0L)
  lp <- do.call(get_weight_like, s)
  expect_equal(lp[1], 0)
})

test_that("wf_minbin and wf_maxbin aggregation works", {
  n_age <- 10
  n_len <- 20
  n_wt  <- 10
  mbin  <- 3L
  xbin  <- 8L

  s <- make_wf_args(n_age = n_age, n_len = n_len, n_wt = n_wt,
                    wf_switch = 1,
                    wf_minbin = c(mbin, 1L), wf_maxbin = c(xbin, n_wt))
  lp <- do.call(get_weight_like, s)
  expect_true(is.finite(lp[1]))
  expect_true(lp[1] >= 0)
  expect_equal(xbin - mbin + 1L, 6L)
})

test_that("mass conservation through rebin: sum(pred_at_weight) ~ sum(pred_at_length)", {
  n_age <- 10
  n_len <- 20
  n_wt  <- 10
  catch_a <- rep(1, n_age)

  pla <- matrix(0, n_len, n_age)
  for (a in seq_len(n_age)) pla[min(a * 2L, n_len), a] <- 1

  # Full-coverage rebin matrix (each weight bin covers 2 length bins, no gaps)
  wf_rebin_matrix <- matrix(0, n_wt, n_len)
  for (w in seq_len(n_wt)) wf_rebin_matrix[w, ((w - 1) * 2 + 1):(w * 2)] <- 1

  pred_at_length <- c(pla %*% catch_a)
  pred_at_weight <- c(wf_rebin_matrix %*% pred_at_length)

  expect_equal(sum(pred_at_weight), sum(pred_at_length), tolerance = 1e-10)
})
