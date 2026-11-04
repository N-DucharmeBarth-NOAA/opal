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

test_that("PR timing smoke covers only Opakapaka and IO YFT objectives and gradients", {
  old <- Sys.getenv(c("OPAL_BENCHMARK_CASES", "OPAL_BENCHMARK_FULL_FITS"), unset = NA_character_)
  on.exit({
    for (name in names(old)) {
      if (is.na(old[[name]])) Sys.unsetenv(name) else
        do.call(Sys.setenv, setNames(list(old[[name]]), name))
    }
  })
  Sys.unsetenv(c("OPAL_BENCHMARK_CASES", "OPAL_BENCHMARK_FULL_FITS"))
  smoke <- benchmark_timing_cases("smoke")
  expect_setequal(smoke, c("opaka_fn", "opaka_gr",
                          "io_yft_fn", "io_yft_gr"))
  expect_length(benchmark_timing_cases("full"), 17L)
  expect_true(all(smoke %in% benchmark_timing_cases("full")))
  expect_false(any(grepl("_(fit|fit_evals|sdreport|trajectory)$", smoke)))
  expect_true(all(c("bet_fn", "bet_gr") %in% benchmark_timing_cases("full")))
  expect_error(benchmark_timing_cases("invalid"), "smoke or full")
  expect_error(benchmark_timing_cases(c("smoke", "full")), "smoke or full")
  expect_error(benchmark_timing_cases(NA_character_), "smoke or full")
  Sys.setenv(OPAL_BENCHMARK_CASES = "bet_fn,bet_gr")
  expect_error(benchmark_timing_cases("smoke"), "No timing cases")
  expect_setequal(benchmark_timing_cases("full"), c("bet_fn", "bet_gr"))
  Sys.setenv(OPAL_BENCHMARK_CASES = "bet_fit")
  expect_error(benchmark_timing_cases("smoke"), "No timing cases")
  expect_identical(benchmark_timing_cases("full"), "bet_fit")
})

test_that("A/B fit counts gate only the Opakapaka fit on PRs", {
  old <- Sys.getenv(c("OPAL_BENCHMARK_CASES", "OPAL_BENCHMARK_FULL_FITS"), unset = NA_character_)
  on.exit({
    for (name in names(old)) {
      if (is.na(old[[name]])) Sys.unsetenv(name) else
        do.call(Sys.setenv, setNames(list(old[[name]]), name))
    }
  })
  Sys.unsetenv("OPAL_BENCHMARK_CASES")
  Sys.setenv(OPAL_BENCHMARK_FULL_FITS = "false")
  expect_identical(benchmark_count_cases("smoke"), "opaka_fit")
  expect_identical(benchmark_count_cases("full"), c("opaka_fit", "bet_fit"))
  Sys.setenv(OPAL_BENCHMARK_FULL_FITS = "true")
  expect_identical(benchmark_count_cases("full"),
                   c("opaka_fit", "bet_fit", "io_yft_fit_evals"))
  expect_identical(benchmark_gated_count_cases(), "opaka_fit")
  Sys.setenv(OPAL_BENCHMARK_CASES = "bet_fn")
  expect_length(benchmark_count_cases("smoke"), 0L)
  expect_error(benchmark_count_cases("invalid"), "smoke or full")
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

test_that("15 percent memory and fit-count regressions fail without waiving the RSS floor", {
  baseline <- data.frame(case = "example",
                         metric = c("mem_alloc_bytes", "fn_evals", "gr_evals",
                                    "iterations", "peak_rss_bytes"),
                         value = c(100, 100, 100, 100, 40 * 1024^2),
                         units = c("bytes", "count", "count", "count", "bytes"))
  observed <- baseline
  observed$value <- baseline$value * 1.15
  # Same-runner A/B fit counts gate; committed-baseline counts are recorded.
  expect_false(any(compare_to_baseline(observed, baseline, tolerance = 0.10,
                                       same_runner = TRUE)$pass))
  committed <- compare_to_baseline(observed, baseline, tolerance = 0.10)
  expect_identical(committed$pass, committed$metric %in% c("fn_evals", "gr_evals", "iterations"))
  peak <- baseline[baseline$metric == "peak_rss_bytes", ]
  observed <- peak
  observed$value <- peak$value + 5 * 1024^2 - 1
  expect_true(compare_to_baseline(observed, peak, tolerance = 0.10)$pass)
  observed$value <- peak$value + 5 * 1024^2
  expect_false(compare_to_baseline(observed, peak, tolerance = 0.10)$pass)
})

test_that("A/B can reuse an idle reference without reusing measured subprocesses", {
  idle_calls <- 0L
  measured_calls <- 0L
  testthat::local_mocked_bindings(r = function(func, args, ...) {
    if (args[[5L]]) {
      idle_calls <<- idle_calls + 1L
      return(list(hwm = 100, metadata = list(), load_mode = "installed"))
    }
    measured_calls <<- measured_calls + 1L
    data.frame(case = args[[2L]], metric = "peak_rss_bytes", value = 200,
                units = "bytes")
  }, .package = "callr")
  first <- run_benchmarks("opaka_fn")
  second <- run_benchmarks("bet_fn", idle = attr(first, "benchmark_idle"))
  expect_identical(idle_calls, 1L)
  expect_identical(measured_calls, 2L)
  expect_equal(first$value, 100)
  expect_equal(second$value, 100)
  expect_identical(attr(second, "benchmark_metadata")$cases, "bet_fn")
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
