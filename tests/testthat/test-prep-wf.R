# Tests for prep_wf_data()

# Shared fixture ----------------------------------------------------------------
# Mirrors make_lf_data() in test-prep-lf.R: load wcpo_bet_data + wcpo_bet_wf,
# add weight bin scalars, pivot to wide format, and call prep_wf_data().

make_wf_data <- function(wf_keep_fisheries = c(8, 9), ...) {
  data(wcpo_bet_data, package = "opal", envir = environment())
  data(wcpo_bet_wf,   package = "opal", envir = environment())
  d <- wcpo_bet_data

  # Add weight bin scalars (prerequisite for prep_wf_data)
  d$wt_bin_start <- 1
  d$wt_bin_width <- 1
  d$n_wt         <- 200L

  wf_wide <- wcpo_bet_wf |>
    tidyr::pivot_wider(
      id_cols     = c(fishery, year, month, ts),
      names_from  = bin,
      values_from = value,
      values_fill = 0
    ) |>
    dplyr::arrange(fishery, ts)

  prep_wf_data(d, wf_wide, wf_keep_fisheries = wf_keep_fisheries, ...)
}

# ---- 1. Output fields exist and have correct types / dimensions ---------------

test_that("prep_wf_data attaches expected fields to data", {
  d <- make_wf_data()

  expect_true(is.integer(d$wf_switch))
  expect_true(is.numeric(d$wt_lower))
  expect_true(is.numeric(d$wt_upper))
  expect_true(is.numeric(d$wt_mid))
  expect_true(is.numeric(d$wt_bin_edges))
  expect_true(is.matrix(d$wf_rebin_matrix))
  expect_true(is.integer(d$n_wf))
  expect_true(is.matrix(d$wf_obs_in))
  expect_true(is.numeric(d$wf_obs_flat))
  expect_true(is.integer(d$wf_obs_ints))
  expect_true(is.numeric(d$wf_obs_prop))
  expect_true(is.numeric(d$wf_n))
  expect_true(is.integer(d$wf_fishery))
  expect_true(is.integer(d$wf_fishery_f))
  expect_true(is.integer(d$wf_n_f))
  expect_true(is.integer(d$wf_year))
  expect_true(is.integer(d$wf_minbin))
  expect_true(is.integer(d$wf_maxbin))
  expect_true(is.numeric(d$wf_var_adjust))

  # Row counts must be consistent
  expect_equal(nrow(d$wf_obs_in),     d$n_wf)
  expect_equal(length(d$wf_n),        d$n_wf)
  expect_equal(length(d$wf_fishery),  d$n_wf)
  expect_equal(length(d$wf_year),     d$n_wf)
  expect_equal(sum(d$wf_n_f),         d$n_wf)

  # Only the two requested fisheries
  expect_setequal(d$wf_fishery_f, c(8L, 9L))

  # wf_var_adjust default: rep(1, n_fishery)
  expect_equal(d$wf_var_adjust, rep(1, d$n_fishery))
})

# ---- 2. Fishery filtering --------------------------------------------------------

test_that("wf_keep_fisheries retains only the requested fisheries", {
  d <- make_wf_data(wf_keep_fisheries = c(8, 9))
  expect_setequal(d$wf_fishery_f, c(8L, 9L))
  expect_true(all(d$wf_fishery %in% c(8L, 9L)))
})

test_that("wf_keep_fisheries = single fishery returns only that fishery", {
  d <- make_wf_data(wf_keep_fisheries = 9)
  expect_setequal(d$wf_fishery_f, 9L)
  expect_true(all(d$wf_fishery == 9L))
})

test_that("wf_keep_fisheries reduces n_wf relative to keeping both fisheries", {
  d_both <- make_wf_data(wf_keep_fisheries = c(8, 9))
  d_one  <- make_wf_data(wf_keep_fisheries = 9)
  expect_lt(d_one$n_wf, d_both$n_wf)
  expect_equal(d_one$n_wf, sum(d_both$wf_fishery == 9L))
})

test_that("wf_keep_fisheries = NULL keeps all fisheries present in wf_wide", {
  data(wcpo_bet_data, package = "opal", envir = environment())
  data(wcpo_bet_wf,   package = "opal", envir = environment())
  d <- wcpo_bet_data
  d$wt_bin_start <- 1
  d$wt_bin_width <- 1
  d$n_wt         <- 200L

  wf_wide <- wcpo_bet_wf |>
    tidyr::pivot_wider(
      id_cols     = c(fishery, year, month, ts),
      names_from  = bin,
      values_from = value,
      values_fill = 0
    ) |>
    dplyr::arrange(fishery, ts)

  d_all <- prep_wf_data(d, wf_wide, wf_keep_fisheries = NULL)
  expected_fisheries <- sort(unique(wcpo_bet_wf$fishery))
  expect_setequal(d_all$wf_fishery_f, expected_fisheries)
})

# ---- 3. Zero-row removal --------------------------------------------------------

test_that("rows with zero total sample size are dropped before processing", {
  data(wcpo_bet_data, package = "opal", envir = environment())
  data(wcpo_bet_wf,   package = "opal", envir = environment())
  d <- wcpo_bet_data
  d$wt_bin_start <- 1
  d$wt_bin_width <- 1
  d$n_wt         <- 200L

  wf_wide <- wcpo_bet_wf |>
    tidyr::pivot_wider(
      id_cols     = c(fishery, year, month, ts),
      names_from  = bin,
      values_from = value,
      values_fill = 0
    ) |>
    dplyr::arrange(fishery, ts) |>
    dplyr::filter(fishery %in% c(8, 9))

  # Inject a zero row for fishery 8 at an unused ts
  zero_row            <- wf_wide[1, ]
  zero_row$ts         <- 9999L
  bin_cols            <- setdiff(names(zero_row), c("fishery", "year", "month", "ts"))
  zero_row[, bin_cols] <- 0
  wf_wide_with_zero   <- rbind(wf_wide, zero_row)

  d_normal <- prep_wf_data(d, wf_wide,           wf_keep_fisheries = c(8, 9))
  d_extras <- prep_wf_data(d, wf_wide_with_zero, wf_keep_fisheries = c(8, 9))

  expect_equal(d_normal$n_wf, d_extras$n_wf)
})

# ---- 4. Variance adjustment ------------------------------------------------------

test_that("wf_var_adjust halves wf_n when adjust = 2 for all fisheries", {
  d_base   <- make_wf_data()
  adjust   <- rep(2, d_base$n_fishery)
  d_scaled <- make_wf_data(wf_var_adjust = adjust)

  expect_equal(d_scaled$wf_n, d_base$wf_n / 2, tolerance = 1e-10)
})

test_that("wf_var_adjust is stored on data and matches the input", {
  d      <- make_wf_data()
  adjust <- seq(0.5, length.out = d$n_fishery, by = 0.1)
  d2     <- make_wf_data(wf_var_adjust = adjust)
  expect_equal(d2$wf_var_adjust, adjust)
})

test_that("wf_var_adjust = 1 (default) is a no-op", {
  d_default  <- make_wf_data()
  d_explicit <- make_wf_data(wf_var_adjust = rep(1, d_default$n_fishery))

  expect_equal(d_default$wf_n,        d_explicit$wf_n)
  expect_equal(d_default$wf_obs_flat, d_explicit$wf_obs_flat)
  expect_equal(d_default$wf_obs_ints, d_explicit$wf_obs_ints)
  expect_equal(d_default$wf_obs_prop, d_explicit$wf_obs_prop)
})

# ---- 5. Rebinning matrix dimensions ---------------------------------------------

test_that("wf_rebin_matrix has dimensions [n_wt x n_len]", {
  d <- make_wf_data()
  expect_equal(dim(d$wf_rebin_matrix), c(d$n_wt, d$n_len))
})

# ---- 6. Rebinning matrix properties ---------------------------------------------

test_that("wf_rebin_matrix has non-negative entries and column sums <= 1", {
  d <- make_wf_data()
  expect_true(all(d$wf_rebin_matrix >= 0))
  expect_true(all(colSums(d$wf_rebin_matrix) <= 1 + 1e-10))
})

# ---- 7. Flattened vector lengths ------------------------------------------------

test_that("wf_obs_flat, wf_obs_prop, wf_obs_ints all have the same length", {
  d <- make_wf_data()
  expect_equal(length(d$wf_obs_flat), length(d$wf_obs_prop))
  expect_equal(length(d$wf_obs_flat), length(d$wf_obs_ints))
})

# ---- 8. Bin alignment error when columns don't match expected structure ---------

test_that("prep_wf_data errors when weight bin columns don't match data scalars", {
  data(wcpo_bet_data, package = "opal", envir = environment())
  data(wcpo_bet_wf,   package = "opal", envir = environment())
  d <- wcpo_bet_data
  # Set wrong wt_bin_start so bins won't match
  d$wt_bin_start <- 5
  d$wt_bin_width <- 1
  d$n_wt         <- 200L

  wf_wide <- wcpo_bet_wf |>
    tidyr::pivot_wider(
      id_cols     = c(fishery, year, month, ts),
      names_from  = bin,
      values_from = value,
      values_fill = 0
    ) |>
    dplyr::arrange(fishery, ts)

  expect_error(prep_wf_data(d, wf_wide, wf_keep_fisheries = c(8, 9)))
})

# ---- 9. Weight bin derivation consistency ---------------------------------------

test_that("wt_lower, wt_upper, wt_mid, wt_bin_edges are consistent with scalars", {
  d <- make_wf_data()

  # wt_lower: starts at wt_bin_start, step = wt_bin_width, length = n_wt
  expect_equal(d$wt_lower[1],   d$wt_bin_start)
  expect_equal(length(d$wt_lower), d$n_wt)
  expect_equal(d$wt_lower[2] - d$wt_lower[1], d$wt_bin_width)

  # wt_upper = wt_lower + wt_bin_width
  expect_equal(d$wt_upper, d$wt_lower + d$wt_bin_width)

  # wt_mid = midpoints
  expect_equal(d$wt_mid, d$wt_lower + d$wt_bin_width / 2)

  # wt_bin_edges has length n_wt + 1
  expect_equal(length(d$wt_bin_edges), d$n_wt + 1)
  expect_equal(d$wt_bin_edges[1], d$wt_lower[1])
  expect_equal(d$wt_bin_edges[d$n_wt + 1], d$wt_upper[d$n_wt])
})
