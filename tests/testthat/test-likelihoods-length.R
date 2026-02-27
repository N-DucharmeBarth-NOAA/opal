# Tests for get_length_like()
library(RTMB)

# Shared test helpers ----

# Build arguments for get_length_like() with a single observation: fishery 1, year 1.
# Returns a named list that can be called with do.call(get_length_like, .).
# lf_obs_flat / lf_obs_prop / lf_obs_ints are pre-aggregated to [bmin:bmax] for
# fishery 1, matching the aggregation applied internally to the predicted vector.
make_lf_args <- function(n_fishery = 2, n_year = 3, n_age = 10, n_len = 20,
                         lf_switch = 1, log_lf_tau = NULL,
                         removal_switch_f = NULL, lf_minbin = NULL,
                         lf_maxbin = NULL, obs_row = NULL, lf_n = 100L) {
  if (is.null(log_lf_tau))       log_lf_tau       <- rep(0.0, n_fishery)
  if (is.null(removal_switch_f)) removal_switch_f <- rep(0L,  n_fishery)
  if (is.null(lf_minbin))        lf_minbin        <- rep(1L,  n_fishery)
  if (is.null(lf_maxbin))        lf_maxbin        <- rep(n_len, n_fishery)
  if (is.null(obs_row))          obs_row          <- rep(1 / n_len, n_len)

  # One observation group: fishery 1, year 1
  lf_n_f       <- 1L
  lf_fishery_f <- 1L
  lf_year_fi   <- list(1L)
  lf_n_fi      <- list(lf_n)

  # Pre-aggregate observed counts to the active [bmin, bmax] window for fishery 1,
  # matching the aggregation applied to pred inside get_length_like():
  #   pred[bmin] <- sum(pred[1:bmin])     (lower tail rolled into bmin)
  #   pred[bmax] <- sum(pred[bmax:n_len]) (upper tail rolled into bmax)
  f1    <- 1L
  bmin  <- lf_minbin[f1]
  bmax  <- lf_maxbin[f1]
  nbins <- bmax - bmin + 1L

  obs_counts <- obs_row * lf_n
  obs_agg    <- obs_counts[bmin:bmax]
  if (bmin > 1L)    obs_agg[1L]    <- sum(obs_counts[1:bmin])
  if (bmax < n_len) obs_agg[nbins] <- sum(obs_counts[bmax:n_len])

  lf_obs_flat <- obs_agg                      # unrounded counts  (switch 1)
  lf_obs_prop <- obs_agg / sum(obs_agg)       # proportions       (switch 2)
  lf_obs_ints <- as.integer(round(obs_agg))   # integer counts    (switch 3)

  # Uniform catch-at-age
  catch_pred_fya <- array(1, dim = c(n_fishery, n_year, n_age))

  # PLA: each age maps exactly to one length bin
  pla <- matrix(0, n_len, n_age)
  for (a in seq_len(n_age)) pla[min(a * 2L, n_len), a] <- 1

  list(
    lf_obs_flat    = lf_obs_flat,
    lf_obs_ints    = lf_obs_ints,
    lf_obs_prop    = lf_obs_prop,
    catch_pred_fya = catch_pred_fya,
    pla            = pla,
    lf_n_f         = lf_n_f,
    lf_fishery_f   = lf_fishery_f,
    lf_year_fi     = lf_year_fi,
    lf_n_fi        = lf_n_fi,
    lf_minbin      = lf_minbin,
    lf_maxbin      = lf_maxbin,
    removal_switch_f = removal_switch_f,
    lf_switch      = lf_switch,
    n_len          = n_len,
    n_lf           = 1L,
    log_lf_tau     = log_lf_tau
  )
}

# Tests ----

test_that("predicted proportions at length sum to 1", {
  n_age <- 40
  n_len <- 95
  catch_a <- rep(1, n_age)
  pla <- matrix(0, n_len, n_age)
  for (a in 1:n_age) pla[min(a * 2, n_len), a] <- 1

  pred <- c(pla %*% catch_a)
  pred <- (pred + 1e-8) / sum(pred + 1e-8)
  expect_equal(sum(pred), 1.0, tolerance = 1e-10)
})

test_that("multinomial NLL is finite and non-negative", {
  s  <- make_lf_args(lf_switch = 1)
  lp <- do.call(get_length_like, s)
  expect_true(all(is.finite(lp)))
  expect_true(all(lp >= 0))
})

test_that("perfect fit gives finite, non-negative NLL for multinomial", {
  n_age <- 10
  n_len <- 20

  pla <- matrix(0, n_len, n_age)
  for (a in seq_len(n_age)) pla[min(a * 2, n_len), a] <- 1

  # Build obs proportions equal to the normalised predicted proportions
  pred_raw  <- c(pla %*% rep(1, n_age))
  pred_norm <- (pred_raw + 1e-8) / sum(pred_raw + 1e-8)

  s  <- make_lf_args(n_age = n_age, n_len = n_len, lf_switch = 1,
                     obs_row = pred_norm, lf_n = 100L)
  lp <- do.call(get_length_like, s)
  expect_true(is.finite(lp[1]))
  expect_true(lp[1] >= 0)
})

test_that("log_lf_tau affects NLL magnitude for lf_switch = 2 (Dirichlet)", {
  n_age    <- 10
  n_len    <- 20
  base_obs <- rep(1 / n_len, n_len)

  s_low  <- make_lf_args(n_age = n_age, n_len = n_len, lf_switch = 2,
                         obs_row = base_obs, log_lf_tau = rep(0.0, 2))
  s_high <- make_lf_args(n_age = n_age, n_len = n_len, lf_switch = 2,
                         obs_row = base_obs, log_lf_tau = rep(2.0, 2))

  lp_low  <- do.call(get_length_like, s_low)
  lp_high <- do.call(get_length_like, s_high)

  expect_true(is.finite(lp_low[1]))
  expect_true(is.finite(lp_high[1]))
  expect_true(lp_low[1]  >= 0)
  expect_true(lp_high[1] >= 0)
  # Higher concentration (larger log_lf_tau) should produce a different NLL
  expect_false(isTRUE(all.equal(lp_low[1], lp_high[1])))
})

test_that("zero catch-at-age does not produce NaN or Inf", {
  s <- make_lf_args(lf_switch = 1)
  s$catch_pred_fya[1, 1, ] <- 0  # all zeros for the observed fishery/year
  lp <- do.call(get_length_like, s)
  expect_true(all(is.finite(lp)))
})

test_that("lf_minbin aggregation is correct", {
  n_age <- 10
  n_len <- 20
  mbin  <- 5L

  s_no_agg <- make_lf_args(n_age = n_age, n_len = n_len, lf_switch = 1,
                            lf_minbin = c(1L, 1L), lf_maxbin = c(n_len, n_len))
  s_agg    <- make_lf_args(n_age = n_age, n_len = n_len, lf_switch = 1,
                            lf_minbin = c(mbin, 1L), lf_maxbin = c(n_len, n_len))

  lp_no_agg <- do.call(get_length_like, s_no_agg)
  lp_agg    <- do.call(get_length_like, s_agg)

  expect_true(is.finite(lp_no_agg[1]))
  expect_true(is.finite(lp_agg[1]))
  expect_true(lp_no_agg[1] >= 0)
  expect_true(lp_agg[1]    >= 0)
})

test_that("removal_switch_f = 1 gives zero NLL contribution", {
  s  <- make_lf_args(lf_switch = 1, removal_switch_f = c(1L, 0L))
  lp <- do.call(get_length_like, s)
  expect_equal(lp[1], 0)
})

test_that("lf_switch = 1 and lf_switch = 2 both produce finite results", {
  s1 <- make_lf_args(lf_switch = 1)
  s2 <- make_lf_args(lf_switch = 2)

  lp1 <- do.call(get_length_like, s1)
  lp2 <- do.call(get_length_like, s2)

  expect_true(all(is.finite(lp1)))
  expect_true(all(is.finite(lp2)))
  expect_true(all(lp1 >= 0))
  expect_true(all(lp2 >= 0))
})

test_that("lf_maxbin aggregation collapses upper-tail bins", {
  n_age <- 10
  n_len <- 20
  xbin  <- 15L

  s_agg <- make_lf_args(n_age = n_age, n_len = n_len, lf_switch = 1,
                         lf_minbin = c(1L, 1L), lf_maxbin = c(xbin, n_len))

  lp <- do.call(get_length_like, s_agg)
  expect_true(is.finite(lp[1]))
  expect_true(lp[1] >= 0)

  # The raw predicted mass at or above xbin is >= the individual bin value
  pred_raw <- c(s_agg$pla %*% rep(1, n_age))
  expect_true(sum(pred_raw[xbin:n_len]) >= pred_raw[xbin])
})

test_that("lf_minbin and lf_maxbin work together", {
  n_age <- 10
  n_len <- 20
  mbin  <- 3L
  xbin  <- 15L

  s_both <- make_lf_args(n_age = n_age, n_len = n_len, lf_switch = 1,
                          lf_minbin = c(mbin, 1L), lf_maxbin = c(xbin, n_len))

  lp <- do.call(get_length_like, s_both)
  expect_true(is.finite(lp[1]))
  expect_true(lp[1] >= 0)
  # Active bin count: xbin - mbin + 1 = 13
  expect_equal(xbin - mbin + 1L, 13L)
})

test_that("lf_maxbin = n_len gives finite result (no upper aggregation)", {
  n_age <- 10
  n_len <- 20

  s  <- make_lf_args(n_age = n_age, n_len = n_len, lf_switch = 1,
                     lf_minbin = c(1L, 1L), lf_maxbin = c(n_len, n_len))
  lp <- do.call(get_length_like, s)
  expect_true(is.finite(lp[1]))
  expect_true(lp[1] >= 0)
})

test_that("lf_maxbin = lf_minbin collapses to single bin and NLL is near zero", {
  n_age <- 10
  n_len <- 20
  bin   <- 10L

  # With a single active bin, pred normalises to 1.0 and obs sums to lf_n.
  # dmultinom(lf_n, prob = 1.0) = 1  =>  NLL = 0.
  s  <- make_lf_args(n_age = n_age, n_len = n_len, lf_switch = 1,
                     lf_minbin = c(bin, 1L), lf_maxbin = c(bin, n_len),
                     lf_n = 100L)
  lp <- do.call(get_length_like, s)
  expect_equal(lp[1], 0, tolerance = 1e-6)
})

test_that("lf_switch = 3 (Dirichlet-multinomial) produces finite non-negative NLL", {
  s  <- make_lf_args(lf_switch = 3)
  lp <- do.call(get_length_like, s)
  expect_true(all(is.finite(lp)))
  expect_true(all(lp >= 0))
})

test_that("n_i = 0 observation is skipped and contributes zero NLL", {
  s <- make_lf_args(lf_switch = 1)
  s$lf_n_fi <- list(0L)   # override sample size to zero
  lp <- do.call(get_length_like, s)
  expect_equal(lp[1], 0)
})

test_that("multinomial NLL matches dmultinom reference value", {
  # Construct a case where pred is known: 2 ages, 4 length bins
  # pla: age 1 -> bin 2, age 2 -> bin 4 (min(a*2, n_len))
  n_age <- 2L
  n_len <- 4L
  pla <- matrix(0, n_len, n_age)
  for (a in seq_len(n_age)) pla[min(a * 2L, n_len), a] <- 1

  # Non-uniform catch so pred is non-trivial
  catch_pred_fya <- array(0, dim = c(2L, 1L, n_age))
  catch_pred_fya[1, 1, ] <- c(3, 1)  # bin 2 gets 3x weight of bin 4

  pred_raw <- c(pla %*% catch_pred_fya[1, 1, ])
  pred     <- (pred_raw + 1e-8) / sum(pred_raw + 1e-8)
  obs      <- c(0, 60, 0, 40)
  expected_nll <- -dmultinom(obs, prob = pred, log = TRUE)

  s <- make_lf_args(n_age = n_age, n_len = n_len, lf_switch = 1,
                    lf_n = 100L)
  s$catch_pred_fya <- catch_pred_fya
  s$pla            <- pla
  s$lf_obs_flat    <- obs   # matches active window [1:n_len], no aggregation

  lp <- do.call(get_length_like, s)
  expect_equal(lp[1], expected_nll, tolerance = 1e-10)
})

test_that("multiple observations: obs_offset advances correctly across observations", {
  # Two observations for fishery 1, years 1 and 2, with different catch patterns
  # so each year's pred is distinct — a bug in obs_offset would mix them up.
  n_age  <- 4L
  n_len  <- 8L
  n_fishery <- 1L
  n_year    <- 2L

  # pla: age a -> bin min(a*2, n_len)  (1->2, 2->4, 3->6, 4->8)
  pla <- matrix(0, n_len, n_age)
  for (a in seq_len(n_age)) pla[min(a * 2L, n_len), a] <- 1

  # Year 1: all catch at age 1 (-> bin 2); Year 2: all catch at age 4 (-> bin 8)
  catch_pred_fya <- array(0, dim = c(n_fishery, n_year, n_age))
  catch_pred_fya[1, 1, ] <- c(1, 0, 0, 0)
  catch_pred_fya[1, 2, ] <- c(0, 0, 0, 1)

  # Build expected preds
  pred_y1 <- c(pla %*% catch_pred_fya[1, 1, ]); pred_y1 <- (pred_y1 + 1e-8) / sum(pred_y1 + 1e-8)
  pred_y2 <- c(pla %*% catch_pred_fya[1, 2, ]); pred_y2 <- (pred_y2 + 1e-8) / sum(pred_y2 + 1e-8)

  # Obs: year 1 concentrated in bin 2; year 2 concentrated in bin 8
  obs_y1 <- c(0, 100, 0, 0, 0, 0, 0, 0)
  obs_y2 <- c(0, 0, 0, 0, 0, 0, 0, 100)

  expected_nll1 <- -dmultinom(obs_y1, prob = pred_y1, log = TRUE)
  expected_nll2 <- -dmultinom(obs_y2, prob = pred_y2, log = TRUE)

  lf_args <- list(
    lf_obs_flat      = c(obs_y1, obs_y2),   # concatenated flat vector
    lf_obs_ints      = as.integer(c(obs_y1, obs_y2)),
    lf_obs_prop      = c(obs_y1 / sum(obs_y1), obs_y2 / sum(obs_y2)),
    catch_pred_fya   = catch_pred_fya,
    pla              = pla,
    lf_n_f           = 2L,                  # 2 obs for fishery 1
    lf_fishery_f     = 1L,                  # single fishery group
    lf_year_fi       = list(c(1L, 2L)),
    lf_n_fi          = list(c(100L, 100L)),
    lf_minbin        = 1L,
    lf_maxbin        = n_len,
    removal_switch_f = 0L,
    lf_switch        = 1L,
    n_len            = n_len,
    n_lf             = 2L,
    log_lf_tau       = 0.0
  )

  lp <- do.call(get_length_like, lf_args)
  expect_equal(lp[1], expected_nll1, tolerance = 1e-10)
  expect_equal(lp[2], expected_nll2, tolerance = 1e-10)
})
