# Deliberate baseline maintenance, never run automatically by benchmark tests.
root <- Sys.getenv("OPAL_BENCHMARK_ROOT", getwd())
Sys.setenv(OPAL_BENCHMARK_ROOT = normalizePath(root, mustWork = TRUE))
source(file.path(root, "tests/testthat/helper-benchmark.R"))
args <- commandArgs(trailingOnly = TRUE)
iterations <- if (length(args)) as.integer(args[[1L]]) else 20L
library <- Sys.getenv("OPAL_BENCHMARK_LIBRARY", "")
if (!nzchar(library)) library <- NULL
results <- run_benchmarks(library = library, iterations = iterations)
metadata <- attr(results, "benchmark_metadata")
metadata$iterations <- iterations
metadata$idle_hwm_bytes <- attr(results, "idle_hwm_bytes")
metadata$library <- library
path <- file.path(root, "tests/testthat/_benchmark/baseline.rds")
dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
if (file.exists(path)) {
  previous <- readRDS(path)
  key <- function(x) paste(x$case, x$metric, sep = "\r")
  retained <- previous$results[!key(previous$results) %in% key(results), ]
  results <- rbind(retained, results)
  history <- previous$metadata$measurement_runs
  if (is.null(history)) history <- list(previous$metadata)
  metadata$measurement_runs <- c(history, list(metadata))
}
instructions_path <- if (length(args) > 1L) args[[2L]] else
  Sys.getenv("OPAL_BENCHMARK_INSTRUCTIONS", "")
if (nzchar(instructions_path)) {
  # Optional second argument (or environment variable): measured instructions.rds
  # from dev/run-benchmark-instructions.R, never fabricated calibration values.
  instruction <- readRDS(instructions_path)
  required <- c("case", "metric", "value", "units")
  x <- instruction$results
  if (!is.data.frame(x) || !all(required %in% names(x)) ||
      any(x$metric != "instructions") || any(!is.finite(x$value)) ||
      any(x$value <= 0) || anyDuplicated(x[c("case", "metric")])) {
    stop("Invalid measured instruction artifact: ", instructions_path)
  }
  key <- function(x) paste(x$case, x$metric, sep = "\r")
  results <- rbind(results[!key(results) %in% key(x), ], x[required])
  metadata$instruction_measurement <- list(
    artifact = normalizePath(instructions_path, mustWork = TRUE),
    md5 = tools::md5sum(instructions_path), metadata = instruction$metadata
  )
}
saveRDS(list(results = results, metadata = metadata), path, version = 3)
message("Measured baseline written to ", path, "; review provenance before committing.")
