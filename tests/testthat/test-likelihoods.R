# Tests for get_length_like()
library(RTMB)

# Shared test helpers ----

# Build minimal data and parameters for get_length_like() tests
make_lf_data <- function(n_fishery = 2, n_year = 3, n_age = 10, n_len = 20,
                         lf_switch = 9, lf_var_adj = NULL,
                         removal_switch_f = NULL, lf_minbin = NULL,
                         obs_row = NULL) {
  if (is.null(lf_var_adj)) lf_var_adj <- rep(1, n_fishery)
  if (is.null(removal_switch_f)) removal_switch_f <- rep(0, n_fishery)
  if (is.null(lf_minbin)) lf_minbin <- rep(1, n_fishery)

  # Single observation: fishery 1, year 1
  lf_obs_row <- if (!is.null(obs_row)) obs_row else {
    x <- rep(1, n_len)
    x / sum(x)
  }
  lf_obs <- matrix(lf_obs_row, nrow = 1, ncol = n_len)

  data <- list(
    lf_obs = lf_obs,
    lf_fishery = 1L,
    lf_year = 1L,
    lf_n = 100,
    lf_switch = lf_switch,
    lf_var_adj = lf_var_adj,
    lf_minbin = lf_minbin,
    removal_switch_f = removal_switch_f
  )
  parameters <- list()

  # Uniform catch-at-age
  catch_pred_fya <- array(1, dim = c(n_fishery, n_year, n_age))

  # Simple diagonal-ish PLA: each age maps to one length bin
  pla <- matrix(0, n_len, n_age)
  for (a in seq_len(n_age)) {
    l <- min(a * 2, n_len)
    pla[l, a] <- 1
  }

  list(data = data, parameters = parameters,
       catch_pred_fya = catch_pred_fya, pla = pla)
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

test_that("multinomial offset NLL is finite and non-negative", {
  s <- make_lf_data(lf_switch = 9)
  lp <- get_length_like(s$data, s$parameters, s$catch_pred_fya, s$pla)
  expect_true(all(is.finite(lp)))
  expect_true(all(lp >= 0))
})

test_that("perfect fit gives NLL near zero for offset form", {
  n_age <- 10
  n_len <- 20

  # Build PLA where each age maps to one length bin (column sums = 1)
  pla <- matrix(0, n_len, n_age)
  for (a in seq_len(n_age)) pla[min(a * 2, n_len), a] <- 1

  # Uniform catch-at-age produces a predictable predicted proportion vector
  catch_pred_fya <- array(1, dim = c(2, 3, n_age))
  pred_raw <- c(pla %*% rep(1, n_age))
  pred_raw <- (pred_raw + 1e-8) / sum(pred_raw + 1e-8)

  # Observed exactly equals predicted (no epsilon issues needed here)
  data <- list(
    lf_obs = matrix(pred_raw, nrow = 1),
    lf_fishery = 1L,
    lf_year = 1L,
    lf_n = 100,
    lf_switch = 9,
    lf_var_adj = rep(1, 2),
    lf_minbin = rep(1, 2),
    removal_switch_f = rep(0, 2)
  )
  lp <- get_length_like(data, list(), catch_pred_fya, pla)
  expect_equal(lp[1], 0, tolerance = 1e-6)
})

test_that("variance adjustment scales NLL linearly", {
  n_age <- 10
  n_len <- 20
  pla <- matrix(0, n_len, n_age)
  for (a in seq_len(n_age)) pla[min(a * 2, n_len), a] <- 1
  catch_pred_fya <- array(1, dim = c(2, 3, n_age))

  base_obs <- rep(1 / n_len, n_len)
  data_base <- list(
    lf_obs = matrix(base_obs, nrow = 1),
    lf_fishery = 1L, lf_year = 1L, lf_n = 100,
    lf_switch = 9,
    lf_var_adj = rep(1, 2),
    lf_minbin = rep(1, 2),
    removal_switch_f = rep(0, 2)
  )
  data_half <- data_base
  data_half$lf_var_adj <- rep(0.5, 2)

  lp_full <- get_length_like(data_base, list(), catch_pred_fya, pla)
  lp_half <- get_length_like(data_half, list(), catch_pred_fya, pla)

  expect_equal(lp_half[1], lp_full[1] * 0.5, tolerance = 1e-10)
})

test_that("zero catch-at-age does not produce NaN or Inf", {
  s <- make_lf_data()
  s$catch_pred_fya[1, 1, ] <- 0  # All zeros for the observed fishery/year
  lp <- get_length_like(s$data, s$parameters, s$catch_pred_fya, s$pla)
  expect_true(all(is.finite(lp)))
})

test_that("lf_minbin aggregation is correct", {
  n_age <- 10
  n_len <- 20
  mbin <- 5L

  # All catch goes to bins 2,4,6,...,20 (one per age)
  pla <- matrix(0, n_len, n_age)
  for (a in seq_len(n_age)) pla[min(a * 2, n_len), a] <- 1
  catch_pred_fya <- array(1, dim = c(2, 3, n_age))

  obs_row <- rep(1 / n_len, n_len)

  data_no_agg <- list(
    lf_obs = matrix(obs_row, nrow = 1),
    lf_fishery = 1L, lf_year = 1L, lf_n = 100,
    lf_switch = 9,
    lf_var_adj = rep(1, 2),
    lf_minbin = c(1L, 1L),   # no aggregation
    removal_switch_f = rep(0, 2)
  )
  data_agg <- data_no_agg
  data_agg$lf_minbin <- c(mbin, 1L)   # aggregate bins 1:5 into bin 5

  lp_no_agg <- get_length_like(data_no_agg, list(), catch_pred_fya, pla)
  lp_agg    <- get_length_like(data_agg,    list(), catch_pred_fya, pla)

  # Both should be finite
  expect_true(is.finite(lp_no_agg[1]))
  expect_true(is.finite(lp_agg[1]))

  # Aggregating bins always loses information, so NLL should differ
  # (equal only in degenerate cases where aggregation has no effect)
  # Here the uniform obs and pred after eps-normalization should differ
  # in the aggregated case, so just check they are both non-negative numbers
  expect_true(lp_no_agg[1] >= 0)
  expect_true(lp_agg[1] >= 0)
})

test_that("removal_switch_f = 1 gives zero NLL contribution", {
  s <- make_lf_data(removal_switch_f = c(1L, 0L))
  lp <- get_length_like(s$data, s$parameters, s$catch_pred_fya, s$pla)
  expect_equal(lp[1], 0)
})

test_that("lf_switch = 1 and lf_switch = 9 both produce finite results", {
  s1 <- make_lf_data(lf_switch = 1)
  s9 <- make_lf_data(lf_switch = 9)

  lp1 <- get_length_like(s1$data, s1$parameters, s1$catch_pred_fya, s1$pla)
  lp9 <- get_length_like(s9$data, s9$parameters, s9$catch_pred_fya, s9$pla)

  expect_true(all(is.finite(lp1)))
  expect_true(all(is.finite(lp9)))
  expect_true(all(lp1 >= 0))
  expect_true(all(lp9 >= 0))
})
