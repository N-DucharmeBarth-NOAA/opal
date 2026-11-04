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
  source(file.path(root, "tests/testthat/helper-benchmark.R"))
  source(file.path(root, "tests/testthat/helper-performance.R"))
  library_path <- Sys.getenv("OPAL_BENCHMARK_LIBRARY")
  if (!nzchar(library_path)) stop("OPAL_BENCHMARK_LIBRARY is required")
  dependency_libraries <- .libPaths()
  Sys.setenv(OPAL_BENCHMARK_DEPENDENCY_LIBRARIES =
               paste(dependency_libraries, collapse = .Platform$path.sep))
  .libPaths(c(library_path, dependency_libraries, .libPaths()))
  benchmark_load(library_path)
  out <- Sys.getenv("OPAL_BENCHMARK_OUT", file.path(root, "benchmark-out"))
  dir.create(out, showWarnings = FALSE, recursive = TRUE)
  cases <- benchmark_instruction_cases()
  selection <- Sys.getenv("OPAL_BENCHMARK_CASES")
  if (nzchar(selection)) {
    selected <- strsplit(selection, ",", fixed = TRUE)[[1]]
    if (any(!selected %in% cases)) stop("Unknown instruction benchmark case")
    cases <- intersect(cases, selected)
  }
  cachegrind <- function(label, worker_args) {
    prefix <- file.path(out, label)
    status <- system2("valgrind", c(
      "--tool=cachegrind", "--cache-sim=no", "--branch-sim=no",
      shQuote(paste0("--cachegrind-out-file=", prefix, ".cachegrind")),
      shQuote(file.path(R.home("bin"), "exec", "R")),
      "--vanilla", "--slave", "-f",
      shQuote(file.path(root, "dev/run-benchmark-instructions.R")),
      "--args", worker_args
    ), stdout = paste0(prefix, ".stdout"), stderr = paste0(prefix, ".stderr"),
    env = paste0("R_HOME=", shQuote(R.home())))
    if (status != 0L) stop("Cachegrind failed: see ", prefix, ".stderr")
    benchmark_instruction_count(paste0(prefix, ".cachegrind"))
  }
  # TEMPORARY validation-only (#164): size the injected sum() from its measured
  # instructions per element and the committed per-call baseline count.
  slowdown <- benchmark_inject_value("OPAL_BENCHMARK_INJECT_SLOWDOWN")
  injection_sizes <- NULL
  if (slowdown) {
    probe_n <- 1e6
    probe <- vapply(c(0L, 3L), function(iterations) {
      cachegrind(paste0("inject-probe-", iterations),
                 c("--inject-probe", format(probe_n, scientific = FALSE), iterations))
    }, numeric(1))
    per_element <- (probe[2] - probe[1]) / (3 * probe_n)
    committed <- readRDS(file.path(root, "tests/testthat/_benchmark/baseline.rds"))$results
    committed <- committed[committed$metric == "instructions", ]
    injection_sizes <- vapply(cases, function(case) {
      round(slowdown * committed$value[committed$case == case] / per_element)
    }, numeric(1))
  }
  results <- lapply(cases, function(case) {
    if (!is.null(injection_sizes)) {
      Sys.setenv(OPAL_BENCHMARK_INJECT_SIZE = format(injection_sizes[[case]], scientific = FALSE))
      on.exit(Sys.unsetenv("OPAL_BENCHMARK_INJECT_SIZE"), add = TRUE)
    }
    counts <- vapply(c(0L, 3L), function(iterations) {
      cachegrind(paste0(case, "-", iterations), c("--worker", shQuote(case), iterations))
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
  if (!is.null(injection_sizes)) {
    metadata$injection <- list(slowdown = slowdown, per_element = per_element,
                               sizes = injection_sizes)
  }
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
  } else if (length(args) && args[1] == "--inject-probe") {
    # TEMPORARY validation-only (#164): instructions of the injected sum().
    probe <- rep(1, as.numeric(args[2]))
    for (i in seq_len(as.integer(args[3]))) sum(probe)
  } else {
    benchmark_instructions_main()
  }
}
