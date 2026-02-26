# Tests for prep_lf_data()

# Shared fixture ----------------------------------------------------------------
# Mirrors the setup in vignettes/bet.Rmd: load wcpo_bet_data + wcpo_bet_lf,
# pivot to wide format, and call prep_lf_data() with lf_keep_fisheries = c(8, 9).

make_lf_data <- function(lf_keep_fisheries = c(8, 9), ...) {
  data(wcpo_bet_data,  package = "opal", envir = environment())
  data(wcpo_bet_lf,    package = "opal", envir = environment())
  d <- wcpo_bet_data

  lf_wide <- wcpo_bet_lf |>
    tidyr::pivot_wider(
      id_cols     = c(fishery, year, month, ts),
      names_from  = bin,
      values_from = value,
      values_fill = 0
    ) |>
    dplyr::arrange(fishery, ts)

  prep_lf_data(d, lf_wide, lf_keep_fisheries = lf_keep_fisheries, ...)
}

# ---- 1. Output fields exist and have correct types / dimensions ---------------

test_that("prep_lf_data attaches expected fields to data", {
  d <- make_lf_data()

  expect_true(is.numeric(d$lf_obs_flat))
  expect_true(is.numeric(d$lf_obs_prop))
  expect_true(is.integer(d$lf_obs_ints))
  expect_true(is.numeric(d$lf_n))
  expect_true(is.integer(d$lf_fishery))
  expect_true(is.integer(d$lf_year))
  expect_true(is.integer(d$lf_season))
  expect_true(is.integer(d$lf_n_f))
  expect_true(is.integer(d$lf_fishery_f))
  expect_true(is.matrix(d$lf_obs_in))

  # lf_obs_flat / lf_obs_prop / lf_obs_ints must be the same length
  expect_equal(length(d$lf_obs_flat), length(d$lf_obs_prop))
  expect_equal(length(d$lf_obs_flat), length(d$lf_obs_ints))

  # Row counts must be consistent
  expect_equal(nrow(d$lf_obs_in), d$n_lf)
  expect_equal(length(d$lf_n),    d$n_lf)
  expect_equal(length(d$lf_fishery), d$n_lf)
  expect_equal(length(d$lf_year),    d$n_lf)
  expect_equal(length(d$lf_season),  d$n_lf)
  expect_equal(sum(d$lf_n_f),        d$n_lf)

  # Only the two requested fisheries
  expect_setequal(d$lf_fishery_f, c(8L, 9L))

  # lf_var_adjust default: rep(1, n_fishery)
  expect_equal(d$lf_var_adjust, rep(1, d$n_fishery))
})

# ---- 2. Golden-output regression: lf_obs_flat[1:200] -------------------------

test_that("lf_obs_flat[1:200] matches known values", {
  d <- make_lf_data()

  expected_flat <- c(
     0.00000000,  0.00000000,  0.00000000,  0.00000000,  0.00000000,
     0.00000000,  0.00000000,  0.00000000,  0.00000000,  0.00000000,
     0.00000000,  1.01861462,  1.76922569,  0.00000000,  0.00000000,
     0.75061107,  2.52402971,  6.55766744,  9.19502522, 14.26895822,
    26.29484431, 28.55549289, 23.34038495, 25.30035357, 13.19208554,
    13.53962883,  3.76982713,  1.50541508,  3.33759227,  0.00000000,
     0.00000000,  0.00000000,  0.00000000,  0.00000000,  0.00000000,
     0.00000000,  0.00000000,  0.00000000,  0.00000000,  0.00000000,
     0.00000000,  0.00000000,  0.00000000,  0.00000000,  0.00000000,
     0.00000000,  0.00000000,  0.00000000,  0.00000000,  0.00000000,
     0.00000000,  0.00000000,  0.00000000,  0.00000000,  0.00000000,
     0.00000000,  0.00000000,  0.00000000,  0.00000000,  0.00000000,
     0.00000000,  0.00000000,  0.00000000,  0.00000000,  0.00000000,
     0.00000000,  0.00000000,  0.00000000,  0.00000000,  0.00000000,
     0.00000000,  0.00000000,  0.00000000,  0.00000000,  0.00000000,
     0.00000000,  0.00000000,  0.00000000,  0.00000000,  0.00000000,
     0.00000000,  0.00000000,  0.00000000,  0.00000000,  0.00000000,
     0.00000000,  0.00000000,  0.00000000,  0.00000000,  0.00000000,
     0.00000000,  0.00000000,  0.00000000,  0.00000000,  0.00000000,
     0.00000000,  0.00000000,  0.00000000,  0.00000000,  0.00000000,
     0.00000000,  0.00000000,  0.00000000,  0.00000000,  0.00000000,
     0.00000000,  0.00000000,  0.00000000,  0.00000000,  0.77526480,
     0.07417891,  5.31169023, 25.87370641, 21.46456085, 61.48681472,
    73.40762967, 64.04619473, 81.25245244, 81.52538940, 81.14549726,
    67.39284448, 70.75469280, 61.61939057,  6.79805970,  3.83511902,
     0.47938988,  0.95877975,  0.00000000,  0.00000000,  0.00000000,
     0.00000000,  0.00000000,  0.00000000,  0.00000000,  0.00000000,
     0.00000000,  0.00000000,  0.00000000,  0.00000000,  0.00000000,
     0.00000000,  0.00000000,  0.00000000,  0.00000000,  0.00000000,
     0.00000000,  0.00000000,  0.00000000,  0.00000000,  0.00000000,
     0.00000000,  0.00000000,  0.00000000,  0.00000000,  0.00000000,
     0.00000000,  0.00000000,  0.00000000,  0.00000000,  0.00000000,
     0.00000000,  0.00000000,  0.00000000,  0.00000000,  0.00000000,
     0.00000000,  0.00000000,  0.00000000,  0.00000000,  0.00000000,
     0.00000000,  0.00000000,  0.00000000,  0.00000000,  0.00000000,
     0.00000000,  0.00000000,  0.00000000,  0.00000000,  0.00000000,
     0.00000000,  0.00000000,  0.00000000,  0.00000000,  0.00000000,
     0.00000000,  0.00000000,  0.00000000,  0.00000000,  0.00000000,
     0.00000000,  0.00000000,  0.00000000,  0.00000000,  0.00000000,
     0.00000000,  0.00000000,  0.00000000,  0.00000000,  0.00000000
  )

  expect_equal(d$lf_obs_flat[1:200], expected_flat, tolerance = 1e-5)
})

# ---- 3. Golden-output regression: lf_obs_prop[1:200] -------------------------

test_that("lf_obs_prop[1:200] matches known values", {
  d <- make_lf_data()

  expected_prop <- c(
    9.999991e-09, 9.999991e-09, 9.999991e-09, 9.999991e-09, 9.999991e-09,
    9.999991e-09, 9.999991e-09, 9.999991e-09, 9.999991e-09, 9.999991e-09,
    9.999991e-09, 5.823330e-03, 1.011450e-02, 9.999991e-09, 9.999991e-09,
    4.291180e-03, 1.442964e-02, 3.748955e-02, 5.256707e-02, 8.157424e-02,
    1.503250e-01, 1.632490e-01, 1.334347e-01, 1.446396e-01, 7.541787e-02,
    7.740474e-02, 2.155174e-02, 8.606320e-03, 1.908070e-02, 9.999991e-09,
    9.999991e-09, 9.999991e-09, 9.999991e-09, 9.999991e-09, 9.999991e-09,
    9.999991e-09, 9.999991e-09, 9.999991e-09, 9.999991e-09, 9.999991e-09,
    9.999991e-09, 9.999991e-09, 9.999991e-09, 9.999991e-09, 9.999991e-09,
    9.999991e-09, 9.999991e-09, 9.999991e-09, 9.999991e-09, 9.999991e-09,
    9.999991e-09, 9.999991e-09, 9.999991e-09, 9.999991e-09, 9.999991e-09,
    9.999991e-09, 9.999991e-09, 9.999991e-09, 9.999991e-09, 9.999991e-09,
    9.999991e-09, 9.999991e-09, 9.999991e-09, 9.999991e-09, 9.999991e-09,
    9.999991e-09, 9.999991e-09, 9.999991e-09, 9.999991e-09, 9.999991e-09,
    9.999991e-09, 9.999991e-09, 9.999991e-09, 9.999991e-09, 9.999991e-09,
    9.999991e-09, 9.999991e-09, 9.999991e-09, 9.999991e-09, 9.999991e-09,
    9.999991e-09, 9.999991e-09, 9.999991e-09, 9.999991e-09, 9.999991e-09,
    9.999991e-09, 9.999991e-09, 9.999991e-09, 9.999991e-09, 9.999991e-09,
    9.999991e-09, 9.999991e-09, 9.999991e-09, 9.999991e-09, 9.999991e-09,
    9.999991e-09, 9.999991e-09, 9.999991e-09, 9.999991e-09, 9.999991e-09,
    9.999991e-09, 9.999991e-09, 9.999991e-09, 9.999991e-09, 9.999991e-09,
    9.999991e-09, 9.999991e-09, 9.999991e-09, 9.999991e-09, 1.094704e-03,
    1.047525e-04, 7.500254e-03, 3.653435e-02, 3.030852e-02, 8.682098e-02,
    1.036535e-01, 9.043489e-02, 1.147306e-01, 1.151160e-01, 1.145796e-01,
    9.516045e-02, 9.990747e-02, 8.700818e-02, 9.599046e-03, 5.415297e-03,
    6.769209e-04, 1.353832e-03, 9.999991e-09, 9.999991e-09, 9.999991e-09,
    9.999991e-09, 9.999991e-09, 9.999991e-09, 9.999991e-09, 9.999991e-09,
    9.999991e-09, 9.999991e-09, 9.999991e-09, 9.999991e-09, 9.999991e-09,
    9.999991e-09, 9.999991e-09, 9.999991e-09, 9.999991e-09, 9.999991e-09,
    9.999991e-09, 9.999991e-09, 9.999991e-09, 9.999991e-09, 9.999991e-09,
    9.999991e-09, 9.999991e-09, 9.999991e-09, 9.999991e-09, 9.999991e-09,
    9.999991e-09, 9.999991e-09, 9.999991e-09, 9.999991e-09, 9.999991e-09,
    9.999991e-09, 9.999991e-09, 9.999991e-09, 9.999991e-09, 9.999991e-09,
    9.999991e-09, 9.999991e-09, 9.999991e-09, 9.999991e-09, 9.999991e-09,
    9.999991e-09, 9.999991e-09, 9.999991e-09, 9.999991e-09, 9.999991e-09,
    9.999991e-09, 9.999991e-09, 9.999991e-09, 9.999991e-09, 9.999991e-09,
    9.999991e-09, 9.999991e-09, 9.999991e-09, 9.999991e-09, 9.999991e-09,
    9.999991e-09, 9.999991e-09, 9.999991e-09, 9.999991e-09, 9.999991e-09,
    9.999991e-09, 9.999991e-09, 9.999991e-09, 9.999991e-09, 9.999991e-09,
    9.999991e-09, 9.999991e-09, 9.999991e-09, 9.999991e-09, 9.999991e-09
  )

  expect_equal(d$lf_obs_prop[1:200], expected_prop, tolerance = 1e-6)
})

# ---- 4. Golden-output regression: lf_obs_ints[1:200] -------------------------

test_that("lf_obs_ints[1:200] matches known values", {
  d <- make_lf_data()

  expected_ints <- as.integer(c(
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  2,  0,  0,  1,  3,  7,
     9, 14, 26, 29, 23, 25, 13, 14,  4,  2,  3,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  1,  0,  5, 26, 21, 61, 73, 64, 81, 82, 81, 67, 71, 62,  7,  4,  0,
     1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0
  ))

  expect_equal(d$lf_obs_ints[1:200], expected_ints)
})

# ---- 5. lf_obs_prop sums to ~1 per observation row ---------------------------

test_that("lf_obs_prop sums to approximately 1 for each observation", {
  d <- make_lf_data()
  nbins <- d$lf_nbins
  n_obs <- length(d$lf_obs_prop) / nbins
  mat   <- matrix(d$lf_obs_prop, nrow = n_obs, byrow = TRUE)
  row_sums <- rowSums(mat)
  expect_true(all(abs(row_sums - 1) < 1e-6),
              info = paste("Row sums range:", range(row_sums)))
})

# ---- 6. lf_var_adjust scales lf_n by 1/adjust --------------------------------

test_that("lf_var_adjust halves lf_n when adjust = 2 for all fisheries", {
  d_base   <- make_lf_data()
  adjust   <- rep(2, d_base$n_fishery)
  d_scaled <- make_lf_data(lf_var_adjust = adjust)

  expect_equal(d_scaled$lf_n, d_base$lf_n / 2, tolerance = 1e-10)
})

test_that("lf_var_adjust per-fishery: only targeted fishery lf_n is scaled", {
  d_base <- make_lf_data()

  # Build a n_fishery-length adjust vector with 4x only for fishery 8
  adjust <- rep(1, d_base$n_fishery)
  adjust[8] <- 4
  d_scaled <- make_lf_data(lf_var_adjust = adjust)

  # Rows belonging to fishery 8 should have lf_n scaled by 1/4
  f8_rows <- d_base$lf_fishery == 8L
  expect_equal(d_scaled$lf_n[f8_rows],  d_base$lf_n[f8_rows]  / 4, tolerance = 1e-10)
  # Rows belonging to fishery 9 should be unchanged
  expect_equal(d_scaled$lf_n[!f8_rows], d_base$lf_n[!f8_rows],     tolerance = 1e-10)
})

test_that("lf_var_adjust = 1 (default) is a no-op", {
  d_default  <- make_lf_data()
  d_explicit <- make_lf_data(lf_var_adjust = rep(1, d_default$n_fishery))

  expect_equal(d_default$lf_n,        d_explicit$lf_n)
  expect_equal(d_default$lf_obs_flat, d_explicit$lf_obs_flat)
  expect_equal(d_default$lf_obs_ints, d_explicit$lf_obs_ints)
  expect_equal(d_default$lf_obs_prop, d_explicit$lf_obs_prop)
})

# ---- 7. Zero-sample-size rows are silently dropped ---------------------------

test_that("rows with zero total sample size are dropped before processing", {
  data(wcpo_bet_data, package = "opal", envir = environment())
  data(wcpo_bet_lf,   package = "opal", envir = environment())
  d <- wcpo_bet_data

  lf_wide <- wcpo_bet_lf |>
    tidyr::pivot_wider(
      id_cols     = c(fishery, year, month, ts),
      names_from  = bin,
      values_from = value,
      values_fill = 0
    ) |>
    dplyr::arrange(fishery, ts) |>
    dplyr::filter(fishery %in% c(8, 9))

  # Inject a zero row for fishery 8 at an unused ts
  zero_row            <- lf_wide[1, ]
  zero_row$ts         <- 9999L
  bin_cols            <- setdiff(names(zero_row), c("fishery", "year", "month", "ts"))
  zero_row[, bin_cols] <- 0
  lf_wide_with_zero   <- rbind(lf_wide, zero_row)

  d_normal <- prep_lf_data(d, lf_wide,           lf_keep_fisheries = c(8, 9))
  d_extras <- prep_lf_data(d, lf_wide_with_zero, lf_keep_fisheries = c(8, 9))

  expect_equal(d_normal$n_lf, d_extras$n_lf)
})

# ---- 8. lf_var_adjust is stored on the returned data object ------------------

# ---- 9. lf_keep_fisheries ----------------------------------------------------

test_that("lf_keep_fisheries retains only the requested fisheries", {
  d <- make_lf_data(lf_keep_fisheries = c(8, 9))
  expect_setequal(d$lf_fishery_f, c(8L, 9L))
  expect_true(all(d$lf_fishery %in% c(8L, 9L)))
})

test_that("lf_keep_fisheries = single fishery returns only that fishery", {
  d <- make_lf_data(lf_keep_fisheries = 9)
  expect_setequal(d$lf_fishery_f, 9L)
  expect_true(all(d$lf_fishery == 9L))
})

test_that("lf_keep_fisheries reduces n_lf relative to keeping both fisheries", {
  d_both <- make_lf_data(lf_keep_fisheries = c(8, 9))
  d_one  <- make_lf_data(lf_keep_fisheries = 9)
  expect_lt(d_one$n_lf, d_both$n_lf)
  # one-fishery n_lf equals the count for fishery 9 in the two-fishery result
  expect_equal(d_one$n_lf, sum(d_both$lf_fishery == 9L))
})

test_that("lf_keep_fisheries = NULL keeps all fisheries present in lf_wide", {
  data(wcpo_bet_data, package = "opal", envir = environment())
  data(wcpo_bet_lf,   package = "opal", envir = environment())
  d <- wcpo_bet_data

  lf_wide <- wcpo_bet_lf |>
    tidyr::pivot_wider(
      id_cols     = c(fishery, year, month, ts),
      names_from  = bin,
      values_from = value,
      values_fill = 0
    ) |>
    dplyr::arrange(fishery, ts)

  d_all <- prep_lf_data(d, lf_wide, lf_keep_fisheries = NULL)
  expected_fisheries <- sort(unique(wcpo_bet_lf$fishery))
  expect_setequal(d_all$lf_fishery_f, expected_fisheries)
})

# ---- 10. lf_var_adjust is stored on the returned data object ------------------

test_that("lf_var_adjust is stored on data and matches the input", {
  d      <- make_lf_data()
  adjust <- seq(0.5, length.out = d$n_fishery, by = 0.1)
  d2     <- make_lf_data(lf_var_adjust = adjust)
  expect_equal(d2$lf_var_adjust, adjust)
})
