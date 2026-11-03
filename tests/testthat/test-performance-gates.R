test_that("benchmark tolerance rejects invalid input", {
  old <- Sys.getenv("OPAL_BENCHMARK_TOL", unset = NA_character_)
  on.exit(if (is.na(old)) Sys.unsetenv("OPAL_BENCHMARK_TOL") else
    Sys.setenv(OPAL_BENCHMARK_TOL = old))
  for (value in c("garbage", "-0.1", "Inf", "NA")) {
    Sys.setenv(OPAL_BENCHMARK_TOL = value)
    expect_error(benchmark_tolerance(), "finite non-negative")
  }
  Sys.setenv(OPAL_BENCHMARK_TOL = "0.1")
  expect_equal(benchmark_tolerance(), 0.1)
})

test_that("paired timing confidence intervals gate confident regressions", {
  expect_equal(unname(benchmark_ratio_ci(rep(1, 10), rep(1, 10))),
               c(1, 1, 1))
  expect_gt(benchmark_ratio_ci(rep(1, 10), rep(1.15, 10))[["lower"]], 1.10)
  expect_lt(benchmark_ratio_ci(rep(1, 10), rep(0.8, 10))[["upper"]], 1)
  # Pair resampling must preserve shared runner drift.
  expect_equal(unname(benchmark_ratio_ci(1:10, 1.15 * (1:10))),
               rep(1.15, 3))
  expect_error(benchmark_ratio_ci(1:4, 1:4), "five paired")
  expect_error(benchmark_ratio_ci(rep(0, 10), rep(1, 10)), "positive finite")
})

test_that("bootstrap comparisons preserve the caller's random seed", {
  set.seed(163)
  seed <- .Random.seed
  invisible(benchmark_ratio_ci(rep(1, 10), rep(1, 10)))
  expect_identical(.Random.seed, seed)
})

test_that("Cachegrind parser selects instructions rather than cache misses", {
  path <- tempfile()
  on.exit(unlink(path))
  writeLines(c("events: Dr Ir Dw", "summary: 5 123456 7"), path)
  expect_equal(benchmark_instruction_count(path), 123456)
  writeLines("events: Ir", path)
  expect_error(benchmark_instruction_count(path), "Missing")
  writeLines(c("events: Dr", "summary: 5"), path)
  expect_error(benchmark_instruction_count(path), "Invalid")
})
