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
saveRDS(list(results = results, metadata = metadata), path, version = 3)
message("Measured baseline written to ", path, "; review provenance before committing.")
