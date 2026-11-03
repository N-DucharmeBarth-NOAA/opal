# Cachegrind totals include setup; subtract an otherwise identical zero-replay
# process so package loading and tape construction cannot hide replay regressions.
benchmark_instructions_worker <- function(case, iterations) {
  root <- normalizePath(Sys.getenv("OPAL_BENCHMARK_ROOT", "."))
  library <- Sys.getenv("OPAL_BENCHMARK_LIBRARY")
  if (!nzchar(library)) stop("OPAL_BENCHMARK_LIBRARY is required")
  dependencies <- strsplit(Sys.getenv("OPAL_BENCHMARK_DEPENDENCY_LIBRARIES"),
                           .Platform$path.sep, fixed = TRUE)[[1]]
  .libPaths(c(library, dependencies, .libPaths()))
  library(opal)
  source(file.path(root, "tests/testthat/helper-benchmark.R"))
  operation <- benchmark_operation(case, file.path(root, "tests/testthat/_benchmark/trajectories"))
  operation()
  gc()
  for (i in seq_len(iterations)) invisible(operation())
}

benchmark_instructions_main <- function() {
  if (!nzchar(Sys.which("valgrind"))) stop("Valgrind is required")
  root <- normalizePath(Sys.getenv("OPAL_BENCHMARK_ROOT", "."))
  source(file.path(root, "tests/testthat/helper-skip.R"))
  source(file.path(root, "tests/testthat/helper-benchmark.R"))
  source(file.path(root, "tests/testthat/helper-performance.R"))
  Sys.setenv(OPAL_BENCHMARK_DEPENDENCY_LIBRARIES =
               paste(.libPaths(), collapse = .Platform$path.sep))
  out <- Sys.getenv("OPAL_BENCHMARK_OUT", file.path(root, "benchmark-out"))
  dir.create(out, showWarnings = FALSE, recursive = TRUE)
  cases <- c("opaka_fn", "opaka_gr", "io_yft_fn", "io_yft_gr")
  selection <- Sys.getenv("OPAL_BENCHMARK_CASES")
  if (nzchar(selection)) {
    selected <- strsplit(selection, ",", fixed = TRUE)[[1]]
    if (any(!selected %in% cases)) stop("Unknown instruction benchmark case")
    cases <- intersect(cases, selected)
  }
  results <- lapply(cases, function(case) {
    counts <- vapply(c(0L, 3L), function(iterations) {
      prefix <- file.path(out, paste0(case, "-", iterations))
      status <- system2("valgrind", c(
        "--tool=cachegrind", "--cache-sim=no", "--branch-sim=no",
        shQuote(paste0("--cachegrind-out-file=", prefix, ".cachegrind")),
        shQuote(file.path(R.home("bin"), "exec", "R")),
        "--vanilla", "--slave", "-f",
        shQuote(file.path(root, "dev/run-benchmark-instructions.R")),
        "--args", "--worker", shQuote(case), iterations
      ), stdout = paste0(prefix, ".stdout"), stderr = paste0(prefix, ".stderr"),
      env = paste0("R_HOME=", shQuote(R.home())))
      if (status != 0L) stop("Cachegrind failed: see ", prefix, ".stderr")
      benchmark_instruction_count(paste0(prefix, ".cachegrind"))
    }, numeric(1))
    count <- (counts[2] - counts[1]) / 3
    if (!is.finite(count) || count <= 0) stop("Non-positive net instruction count")
    data.frame(case = case, metric = "instructions", value = count,
               units = "instructions/call")
  })
  observed <- do.call(rbind, results)
  metadata <- benchmark_metadata()
  metadata$cachegrind <- list(
    version = system2("valgrind", "--version", stdout = TRUE),
    calls = 3L, method = "Three warmed calls minus identical zero-call process"
  )
  saveRDS(list(results = observed, metadata = metadata),
          file.path(out, "instructions.rds"))
  if (identical(Sys.getenv("OPAL_BENCHMARK_UPDATE"), "true")) {
    baseline <- observed
  } else {
    baseline <- readRDS(file.path(root, "tests/testthat/_benchmark/baseline.rds"))$results
  }
  comparison <- compare_to_baseline(observed, baseline)
  write_benchmark_results(comparison, out)
  failures <- comparison[!comparison$pass, ]
  if (nrow(failures)) {
    stop(paste(sprintf(
      "%s instructions: %.6g vs baseline %.6g (%+.1f%%; limit +%.1f%%)",
      failures$case, failures$observed, failures$baseline,
      100 * (failures$ratio - 1), 100 * benchmark_tolerance()
    ), collapse = "\n"))
  }
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) && args[1] == "--worker") {
    benchmark_instructions_worker(args[2], as.integer(args[3]))
  } else {
    benchmark_instructions_main()
  }
}
