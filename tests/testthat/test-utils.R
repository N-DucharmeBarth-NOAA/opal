# Tests for resolve_bio_vector
library(RTMB)

# Define project paths
proj_dir <- this.path::this.proj()
r_dir <- file.path(proj_dir, "code", "rtmb", "R")

# Source helper functions
suppressMessages({
  source(file.path(r_dir, "selectivity.R"))
  source(file.path(r_dir, "utils.R"))
})

# Shared test fixtures ----

make_pla <- function(n_len = 10, n_age = 5) {
  # Construct a simple PLA: each column is a unit vector peaking at a bin
  # proportional to the age index (columns sum to 1)
  pla <- matrix(0, nrow = n_len, ncol = n_age)
  for (a in seq_len(n_age)) {
    bin <- min(a * 2, n_len)  # place mass near the middle for each age
    pla[bin, a] <- 1
  }
  pla
}

# Test: age-basis vector passed through unchanged ----

test_that("resolve_bio_vector returns age-basis vector unchanged", {
  n_age <- 5
  n_len <- 10
  pla <- make_pla(n_len, n_age)
  vec <- c(0.1, 0.3, 0.5, 0.4, 0.2)

  result <- resolve_bio_vector(vec, n_age, n_len, pla, "test_vec")
  expect_equal(result, vec)
})

# Test: length-basis vector converted to age via PLA ----

test_that("resolve_bio_vector converts length-basis vector to age-basis", {
  n_age <- 40
  n_len <- 95
  # Build a simple diagonal-ish PLA
  len_lower <- seq(10, by = 2, length.out = n_len)
  len_upper <- len_lower + 2
  mu_a <- 30 + (180 - 30) * (1 - exp(-0.2 * (0:(n_age - 1))))
  sd_a <- 0.1 * mu_a
  pla <- get_pla(len_lower, len_upper, mu_a, sd_a)

  # Create a simple length-basis maturity vector (logistic)
  len_mid <- len_lower + 1
  mat_l <- 1 / (1 + exp(-0.1 * (len_mid - 100)))

  result <- resolve_bio_vector(mat_l, n_age, n_len, pla, "maturity")

  # Result should be length n_age
  expect_equal(length(result), n_age)
  # Should equal manual PLA multiplication
  expected <- as.vector(t(pla) %*% mat_l)
  expect_equal(result, expected)
  # Values should remain in [0, 1] for a maturity-at-length input
  expect_true(all(result >= 0 & result <= 1))
})

# Test: informative error when length matches neither n_age nor n_len ----

test_that("resolve_bio_vector errors when length matches neither n_age nor n_len", {
  n_age <- 40
  n_len <- 95
  pla <- matrix(1 / n_len, nrow = n_len, ncol = n_age)  # uniform PLA
  vec_bad <- rep(0.3, 20)  # length 20 matches neither 40 nor 95

  expect_error(
    resolve_bio_vector(vec_bad, n_age, n_len, pla, "bad_vec"),
    regexp = "bad_vec has length 20, which matches neither n_age \\(40\\) nor n_len \\(95\\)"
  )
})

# Test: warning when n_age == n_len ----

test_that("resolve_bio_vector warns when n_age == n_len", {
  n_age <- 10
  n_len <- 10  # equal to n_age
  pla <- diag(n_age)  # identity PLA (trivially: each length bin = one age)
  vec <- rep(0.5, n_age)

  expect_warning(
    resolve_bio_vector(vec, n_age, n_len, pla, "ambiguous_vec"),
    regexp = "ambiguous_vec has length 10 which matches both n_age and n_len"
  )
})

# Test: age-basis vector returned when n_age == n_len (warning fired, no error) ----

test_that("resolve_bio_vector returns vector unchanged when n_age == n_len", {
  n_age <- 10
  n_len <- 10
  pla <- diag(n_age)
  vec <- seq(0.1, 1.0, length.out = n_age)

  result <- suppressWarnings(
    resolve_bio_vector(vec, n_age, n_len, pla, "vec")
  )
  expect_equal(result, vec)
})

# Test: n_len vector is correctly converted even when n_age < n_len ----

test_that("resolve_bio_vector uses PLA for length-basis when n_age < n_len", {
  n_age <- 5
  n_len <- 10
  pla <- make_pla(n_len, n_age)
  vec_l <- rep(1, n_len)  # constant length vector → should give 1 at each age

  result <- resolve_bio_vector(vec_l, n_age, n_len, pla, "ones")
  expected <- as.vector(t(pla) %*% vec_l)
  expect_equal(result, expected)
  expect_equal(length(result), n_age)
})
