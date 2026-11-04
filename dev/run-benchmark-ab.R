# The head checkout supplies the harness for both installed package versions.
benchmark_ab_main <- function() {
  root <- normalizePath(Sys.getenv("OPAL_BENCHMARK_ROOT", "."))
  source(file.path(root, "tests/testthat/helper-benchmark.R"))
  source(file.path(root, "tests/testthat/helper-performance.R"))
  libraries <- c(base = Sys.getenv("OPAL_BENCHMARK_BASE_LIBRARY"),
                 head = Sys.getenv("OPAL_BENCHMARK_HEAD_LIBRARY"))
  if (any(!dir.exists(libraries))) stop("Both installed A/B libraries are required")
  rounds <- as.integer(Sys.getenv("OPAL_BENCHMARK_ROUNDS", "10"))
  if (!is.finite(rounds) || rounds < 5L) stop("At least five A/B rounds are required")
  # A single long fit collects deterministic counts, not a timing CI.
  timing_cases <- benchmark_timing_cases()
  out <- Sys.getenv("OPAL_BENCHMARK_OUT", file.path(root, "benchmark-out"))
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  measurements <- list()
  # Each library needs one fresh idle reference, as in the multi-case runner.
  idle <- setNames(vector("list", length(libraries)), names(libraries))
  for (case in timing_cases) {
    for (round in seq_len(rounds)) {
      # Alternate which revision runs first so within-pair drift cancels.
      order <- if (round %% 2L) names(libraries) else rev(names(libraries))
      for (position in seq_along(order)) {
        version <- order[[position]]
        iterations <- if (grepl("_fit$", case)) 1L else 20L
        value <- run_benchmarks(cases = case, library = libraries[[version]],
                                iterations = iterations, idle = idle[[version]])
        idle[[version]] <- attr(value, "benchmark_idle")
        value$round <- round
        value$version <- version
        value$position <- position
        measurements[[length(measurements) + 1L]] <- value
        # Preserve partial results even when a later benchmark fails.
        write.csv(do.call(rbind, measurements), file.path(out, "rounds.csv"),
                  row.names = FALSE)
      }
    }
  }
  raw <- do.call(rbind, measurements)
  tolerance <- benchmark_tolerance()
  results <- list()
  for (case in unique(raw$case)) {
    values <- raw[raw$case == case, ]
    for (metric in unique(values$metric)) {
      x <- values[values$metric == metric, ]
      base <- x$value[x$version == "base"]
      head <- x$value[x$version == "head"]
      if (metric == "median_seconds") {
        ci <- benchmark_ratio_ci(base, head)
        pass <- ci[["lower"]] <= 1 + tolerance
        if (ci[["ratio"]] > 1 + tolerance && pass) {
          warning(case, ": timing estimate exceeds tolerance but CI overlaps limit")
        }
        if (ci[["ratio"]] < 1 - tolerance) {
          message(case, ": timing improved by more than tolerance; refresh baseline")
        }
        row <- data.frame(case = case, metric = metric,
                          baseline = median(base), observed = median(head),
                          ratio = ci[["ratio"]], pass = pass,
                          lower = ci[["lower"]], upper = ci[["upper"]])
      } else {
        observed <- data.frame(case = case, metric = metric,
                               value = median(head), units = x$units[1])
        baseline <- observed
        baseline$value <- median(base)
        row <- compare_to_baseline(observed, baseline, tolerance)
        row <- row[, c("case", "metric", "baseline", "observed", "ratio", "pass")]
        row$lower <- NA_real_
        row$upper <- NA_real_
      }
      results[[length(results) + 1L]] <- row
    }
  }
  results <- do.call(rbind, results)
  metadata <- idle$head$metadata
  shas <- Sys.getenv(c(base = "OPAL_BENCHMARK_BASE_SHA", head = "OPAL_BENCHMARK_HEAD_SHA"),
                     unset = NA_character_)
  names(shas) <- names(libraries)
  metadata$ab_library_metadata <- benchmark_ab_library_metadata(idle, shas)
  metadata$comparison <- Sys.getenv("OPAL_BENCHMARK_COMPARISON", unset = NA_character_)
  metadata$pr_head_sha <- Sys.getenv("OPAL_BENCHMARK_PR_HEAD_SHA", unset = NA_character_)
  metadata$timing_cases <- timing_cases
  metadata$rounds <- rounds
  metadata$round_order <- "BASE first in odd rounds, HEAD first in even rounds"
  attr(results, "benchmark_metadata") <- metadata
  write_benchmark_results(results, out)
  saveRDS(list(metadata = metadata, libraries = libraries,
               commits = Sys.getenv(c("OPAL_BENCHMARK_BASE_SHA",
                                      "OPAL_BENCHMARK_HEAD_SHA"))),
          file.path(out, "provenance.rds"))
  failures <- results[!results$pass, ]
  if (nrow(failures)) {
    details <- sprintf("%s %s: %.6g vs base %.6g (%+.1f%%; limit +%.1f%%)",
                       failures$case, failures$metric, failures$observed,
                       failures$baseline, 100 * (failures$ratio - 1),
                       100 * tolerance)
    # Set only by the workflow's permission-checked maintainer-label step.
    count_failure <- failures$metric %in% c("fn_evals", "gr_evals", "iterations",
                                           "instructions")
    if (identical(Sys.getenv("OPAL_BENCHMARK_ACCEPTED"), "true") &&
        !any(count_failure)) {
      warning(paste(c("Maintainer accepted performance regression:", details),
                    collapse = "\n"))
    } else {
      stop(paste(details, collapse = "\n"))
    }
  }
  invisible(results)
}

if (sys.nframe() == 0L) benchmark_ab_main()
