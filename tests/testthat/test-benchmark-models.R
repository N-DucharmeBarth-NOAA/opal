# Baseline calibration: 2026-11-03, Ubuntu 22.04, R 4.5.2, RTMB 1.9,
# one OMP/OpenBLAS/MKL thread. Generated explicitly with
# dev/update-benchmark-baseline.R 20 and dev/run-benchmark-instructions.R.
# Hosted memory rows: run 41724716265 attempt 2, artifact 13275492072.
# IO YFT full-fit and instruction rows remain sandbox candidates; metadata
# marks the pending hosted metrics. Complete hosted calibration/noise validation
# before making the checks required.
test_that("benchmark comparison gates memory and structure, not single-run time", {
  old <- Sys.getenv("OPAL_BENCHMARK_TOL", unset = NA_character_)
  on.exit(if (is.na(old)) Sys.unsetenv("OPAL_BENCHMARK_TOL") else
    Sys.setenv(OPAL_BENCHMARK_TOL = old))
  Sys.setenv(OPAL_BENCHMARK_TOL = "0.10")
  x <- data.frame(case = "example", metric = c("median_seconds", "mem_alloc_bytes",
                                             "structural_bytes", "peak_rss_bytes"),
                  value = c(1, 100, 100, 0), units = c("seconds", rep("bytes", 3)))
  observed <- x
  observed$value <- c(100, 109, 111, 4 * 1024^2)
  comparison <- compare_to_baseline(observed, x)
  expect_true(comparison$pass[comparison$metric == "median_seconds"])
  expect_true(comparison$pass[comparison$metric == "mem_alloc_bytes"])
  expect_false(comparison$pass[comparison$metric == "structural_bytes"])
  expect_true(comparison$pass[comparison$metric == "peak_rss_bytes"])
  peak <- x[x$metric == "peak_rss_bytes", ]
  peak$value <- 20 * 1024^2
  raised <- peak
  raised$value <- 24 * 1024^2
  expect_true(compare_to_baseline(raised, peak)$pass)
  raised$value <- 26 * 1024^2
  expect_false(compare_to_baseline(raised, peak)$pass)
  peak$value <- 0
  raised$value <- 5.1 * 1024^2
  expect_false(compare_to_baseline(raised, peak)$pass)
  boundary <- x[x$metric == "mem_alloc_bytes", ]
  raised <- boundary
  raised$value <- 110
  expect_true(compare_to_baseline(raised, boundary)$pass)
  boundary$metric <- raised$metric <- "fn_evals"
  raised$value <- 111
  # Committed-baseline fit counts depend on the runner CPU, so they are recorded;
  # same-runner A/B comparisons gate them.
  expect_true(compare_to_baseline(raised, boundary)$pass)
  expect_false(compare_to_baseline(raised, boundary, recorded = character())$pass)
  improved <- x
  improved$value[2] <- 80
  expect_message(compare_to_baseline(improved, x), "refresh.*baseline")
  missing <- compare_to_baseline(observed, x[-2, ])
  expect_false(missing$pass[missing$metric == "mem_alloc_bytes"])
  expect_true(is.na(missing$baseline[missing$metric == "mem_alloc_bytes"]))
  expect_error(compare_to_baseline(rbind(x, x), x), "Duplicate")
  expect_error(compare_to_baseline(x, x, tolerance = -1), "Invalid tolerance")
  invalid <- x
  invalid$value[2] <- NA_real_
  expect_false(compare_to_baseline(invalid, x)$pass[2])
  invalid <- x
  invalid$units[2] <- "seconds"
  expect_false(compare_to_baseline(invalid, x)$pass[2])
})

test_that("benchmark operations cover the required assessment workloads", {
  old <- Sys.getenv("OPAL_BENCHMARK_CASES", unset = NA_character_)
  full <- Sys.getenv("OPAL_BENCHMARK_FULL_FITS", unset = NA_character_)
  on.exit(if (is.na(old)) Sys.unsetenv("OPAL_BENCHMARK_CASES") else
    Sys.setenv(OPAL_BENCHMARK_CASES = old))
  on.exit(if (is.na(full)) Sys.unsetenv("OPAL_BENCHMARK_FULL_FITS") else
    Sys.setenv(OPAL_BENCHMARK_FULL_FITS = full), add = TRUE)
  Sys.setenv(OPAL_BENCHMARK_FULL_FITS = "true")
  Sys.unsetenv("OPAL_BENCHMARK_CASES")
  expect_length(benchmark_cases(), 18L)
  expect_true(all(c("bet_trajectory", "io_yft_fit_evals", "opal_model_plain") %in%
                    benchmark_cases()))
  Sys.setenv(OPAL_BENCHMARK_FULL_FITS = "false")
  expect_false("io_yft_fit_evals" %in% benchmark_cases())
  expect_true(all(c("opaka_fit", "bet_fit") %in% benchmark_cases()))
  Sys.setenv(OPAL_BENCHMARK_CASES = "opaka_fn,io_yft_gr")
  expect_identical(benchmark_cases(), c("opaka_fn", "io_yft_gr"))
  Sys.setenv(OPAL_BENCHMARK_CASES = "invalid")
  expect_error(benchmark_cases(), "Unknown benchmark")
})

test_that("benchmark artifacts retain measured values and provenance", {
  out <- tempfile("benchmark-artifacts-")
  on.exit(unlink(out, recursive = TRUE))
  x <- data.frame(case = "example", metric = "mem_alloc_bytes",
                  value = 123, units = "bytes")
  plain <- x
  metadata <- list(timestamp_utc = "unit test", git_sha = "unit test",
                   R = R.version.string, platform = R.version$platform,
                   packages = character(), source_md5 = character(),
                   scientific_version = "unit test")
  attr(x, "benchmark_metadata") <- metadata
  write_benchmark_results(x, out)
  expect_equal(utils::read.csv(file.path(out, "results.csv")), plain)
  stored <- readRDS(file.path(out, "results.rds"))
  expect_identical(stored$results, x)
  expect_true(all(c("git_sha", "R", "platform", "packages", "source_md5",
                    "scientific_version") %in% names(stored$metadata)))
  expect_true(file.exists(file.path(out, "summary.md")))
  expect_false(file.exists(file.path(out, "baseline.rds")))
  comparison <- compare_to_baseline(x, x)
  attr(comparison, "benchmark_metadata") <- metadata
  write_benchmark_results(comparison, out)
  expect_match(paste(readLines(file.path(out, "summary.md")), collapse = "\n"),
               "baseline.*observed.*ratio.*pass")
  expect_error(run_benchmarks("opaka_fn", iterations = 0), "positive integer")
  expect_error(run_benchmarks(c("opaka_fn", "opaka_fn")), "unique")
})

test_that("two-pass benchmark fitting counts both optimizer passes", {
  calls <- c(fn = 0L, gr = 0L)
  obj <- list(par = c(log_B0 = 2))
  obj$fn <- function(x) {
    calls["fn"] <<- calls["fn"] + 1L
    sum((x - 1)^2)
  }
  obj$gr <- function(x) {
    calls["gr"] <<- calls["gr"] + 1L
    2 * (x - 1)
  }
  fit <- benchmark_fit(obj, list(log_B0 = 2), record = TRUE)
  expect_length(fit$passes, 2L)
  expect_equal(unname(fit$counts[c("fn_evals", "gr_evals")]), unname(calls))
  expect_equal(fit$counts[["iterations"]], sum(vapply(fit$passes, `[[`, 0L, "iterations")))
  expect_equal(length(fit$events), sum(calls))
  expect_equal(fit$opt$par, c(log_B0 = 1), tolerance = 1e-8)
  for (event in fit$events) {
    expected <- if (event$kind == "fn") sum((event$par - 1)^2) else 2 * (event$par - 1)
    expect_equal(event$value, expected)
  }
})

test_that("trajectory sampling preserves exact ordered reference evaluations", {
  events <- lapply(seq_len(20L), function(i) {
    list(kind = if (i %% 2L) "fn" else "gr", par = i, value = i)
  })
  trajectory <- list(events = events)
  sampled <- benchmark_sample_trajectory(trajectory, per_kind = 3L)
  expect_length(sampled$events, 6L)
  expect_identical(sampled$reference_events, events)
  expect_identical(sampled$events, events[sampled$replay_indices])
  expect_identical(sampled$replay_indices, sort(sampled$replay_indices))
  expect_identical(benchmark_sample_trajectory(sampled, 3L), sampled)
  expect_error(benchmark_sample_trajectory(trajectory, 0), "positive integer")
})

test_that("trajectory inputs cannot silently change and references are never synthesized", {
  path <- tempfile("benchmark-trajectories-")
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(path, recursive = TRUE))
  inputs <- list(data = list(), parameters = list(x = 0), map = list())
  expect_error(benchmark_read_trajectory("example", path, inputs), "Missing benchmark trajectory")
  trajectory <- list(inputs = inputs, events = list(list(kind = "fn", par = 0, value = 0)))
  saveRDS(trajectory, file.path(path, "example.rds"))
  expect_identical(benchmark_read_trajectory("example", path, inputs), trajectory)
  expect_error(benchmark_read_trajectory("example", path, list()), "inputs changed")
  trajectory$events[[1]]$par <- NA_real_
  saveRDS(trajectory, file.path(path, "example.rds"))
  expect_error(benchmark_read_trajectory("example", path, inputs), "Invalid benchmark trajectory")
})

test_that("the committed baseline covers every workload and regression metric", {
  # test_path() also resolves inside an R CMD check tarball tree.
  baseline <- readRDS(test_path("_benchmark", "baseline.rds"))
  results <- baseline$results
  cases <- c(
    paste0("opaka_", c("makeadfun", "fn", "gr", "fit", "report", "sdreport")),
    paste0("bet_", c("makeadfun", "fn", "gr", "trajectory", "fit", "sdreport")),
    paste0("io_yft_", c("makeadfun", "fn", "gr", "trajectory", "fit_evals")),
    "opal_model_plain"
  )
  expect_setequal(unique(results$case), cases)
  expect_false(anyDuplicated(results[c("case", "metric")]) > 0L)
  expect_true(all(is.finite(results$value) & results$value >= 0))
  for (case in cases) {
    expect_true(all(c("median_seconds", "mem_alloc_bytes", "structural_bytes",
                      "peak_rss_bytes") %in% results$metric[results$case == case]),
                info = case)
  }
  for (case in c("opaka_fit", "bet_fit", "io_yft_fit_evals")) {
    expect_true(all(c("fn_evals", "gr_evals", "iterations") %in%
                      results$metric[results$case == case]), info = case)
  }
  expect_setequal(results$case[results$metric == "instructions"],
                  c("opaka_fn", "opaka_gr", "io_yft_fn", "io_yft_gr"))
  expect_true(all(c("git_sha", "R", "platform", "packages", "source_md5",
                    "instruction_measurement") %in% names(baseline$metadata)))
  expect_no_error(benchmark_check_baseline(results))
})

test_that("baseline merges require complete instruction artifacts and coverage", {
  baseline <- readRDS(test_path("_benchmark", "baseline.rds"))
  instructions <- list(
    results = baseline$results[baseline$results$metric == "instructions", ],
    metadata = list(git_sha = "instructions"))
  merged <- benchmark_merge_instructions(baseline$results, list(git_sha = "memory"),
                                         instructions, list(md5 = "abc"))
  expect_no_error(benchmark_check_baseline(merged$results))
  expect_identical(merged$metadata$instruction_measurement,
                   list(md5 = "abc", metadata = list(git_sha = "instructions")))
  partial <- instructions
  partial$results <- partial$results[partial$results$case != "io_yft_gr", ]
  expect_error(benchmark_merge_instructions(baseline$results, list(), partial),
               "exactly the cases")
  invalid <- instructions
  invalid$results$value[1] <- 0
  expect_error(benchmark_merge_instructions(baseline$results, list(), invalid),
               "Invalid measured instruction")
  incomplete <- baseline$results[!(baseline$results$case == "io_yft_fit_evals" &
                                     baseline$results$metric == "fn_evals"), ]
  expect_error(benchmark_check_baseline(incomplete), "io_yft_fit_evals fn_evals")
  expect_error(benchmark_check_baseline(
    baseline$results[baseline$results$case != "bet_fit", ]), "bet_fit")
  extra <- rbind(baseline$results, data.frame(case = "unknown", metric = "median_seconds",
                                              value = 1, units = "seconds"))
  expect_error(benchmark_check_baseline(extra), "Unknown baseline case")
})

test_that("A/B provenance records each library's revision, not the harness checkout", {
  record <- function(sha) list(metadata = list(git_sha = sha, git_dirty = character(),
                                               source_md5 = c(a = "1"),
                                               packages = c(opal = "0.1.3")))
  metadata <- benchmark_ab_library_metadata(
    list(base = record("harness"), head = record("harness")),
    c(base = "base-sha", head = "head-sha"))
  expect_identical(names(metadata), c("base", "head"))
  expect_identical(metadata$base$git_sha, "base-sha")
  expect_identical(metadata$head$git_sha, "head-sha")
  expect_identical(metadata$base$harness_checkout$git_sha, "harness")
  expect_null(metadata$base$source_md5)
  expect_identical(metadata$head$packages, c(opal = "0.1.3"))
})

test_that("memory-only baseline refreshes retain instruction provenance", {
  instruction_measurement <- list(artifact = "instructions.rds", md5 = "abc",
                                  metadata = list(git_sha = "instructions"))
  refreshed <- benchmark_retain_instruction_metadata(
    list(git_sha = "memory"), list(instruction_measurement = instruction_measurement), "")
  expect_identical(refreshed$instruction_measurement, instruction_measurement)
  replaced <- benchmark_retain_instruction_metadata(
    list(git_sha = "memory"), list(instruction_measurement = instruction_measurement),
    "new-instructions.rds")
  expect_null(replaced$instruction_measurement)
})

test_that("model performance stays within the explicit measured baseline", {
  skip_if_not_benchmark()
  skip_if(Sys.getenv("OPAL_BENCHMARK_LAYER") %in% c("instructions", "timing-ab"),
          "This layer is handled by its explicit development script.")
  root <- benchmark_root()
  library <- Sys.getenv("OPAL_BENCHMARK_LIBRARY", "")
  if (!nzchar(library)) library <- NULL
  iterations <- as.integer(Sys.getenv("OPAL_BENCHMARK_ITERATIONS", "20"))
  results <- run_benchmarks(library = library, iterations = iterations)
  out <- Sys.getenv("OPAL_BENCHMARK_OUT", file.path(root, "benchmark-out"))
  write_benchmark_results(results, out)
  # Update mode produces reviewable artifacts, never overwrites the baseline.
  if (identical(Sys.getenv("OPAL_BENCHMARK_UPDATE"), "true")) {
    saveRDS(list(results = results, metadata = attr(results, "benchmark_metadata")),
            file.path(out, "baseline.rds"), version = 3)
    expect_true(file.exists(file.path(out, "baseline.rds")))
    return(invisible(NULL))
  }
  path <- file.path(root, "tests/testthat/_benchmark/baseline.rds")
  expect_true(file.exists(path), info = "Explicitly run dev/update-benchmark-baseline.R first.")
  if (!file.exists(path)) return(invisible(NULL))
  comparison <- compare_to_baseline(results, readRDS(path))
  utils::write.csv(results, file.path(out, "measurements.csv"), row.names = FALSE)
  attr(comparison, "benchmark_metadata") <- attr(results, "benchmark_metadata")
  attr(comparison, "idle_hwm_bytes") <- attr(results, "idle_hwm_bytes")
  write_benchmark_results(comparison, out)
  utils::write.csv(comparison, file.path(out, "comparison.csv"), row.names = FALSE)
  failures <- comparison[!comparison$pass, ]
  details <- sprintf("%s %s: %.6g vs baseline %.6g (%+.1f%%; limit +%.1f%%)",
                     failures$case, failures$metric, failures$observed,
                     failures$baseline, 100 * (failures$ratio - 1),
                     100 * benchmark_tolerance())
  expect_true(all(comparison$pass), info = paste(details, collapse = "\n"))
})
