# Combine reviewed hosted artifacts into the committed baseline. Never run by tests.
# Usage: Rscript dev/combine-benchmark-baseline.R memory/baseline.rds \
#   instructions/instructions.rds --hosted <run/artifact URL> [<URL> ...]
root <- Sys.getenv("OPAL_BENCHMARK_ROOT", getwd())
Sys.setenv(OPAL_BENCHMARK_ROOT = normalizePath(root, mustWork = TRUE))
source(file.path(root, "tests/testthat/helper-benchmark.R"))
args <- commandArgs(trailingOnly = TRUE)
flag <- match("--hosted", args)
if (is.na(flag) || flag != 3L || length(args) < 4L) {
  stop("Supply the memory baseline.rds, instructions.rds and --hosted with ",
       "the reviewed run/artifact URLs.")
}
memory_path <- args[[1L]]
instructions_path <- args[[2L]]
memory <- readRDS(memory_path)
if (!is.list(memory) || !is.data.frame(memory$results) || is.null(memory$metadata)) {
  stop("Invalid memory baseline artifact: ", memory_path)
}
if (any(memory$results$metric == "instructions")) {
  stop("The memory artifact must not contain instruction rows: ", memory_path)
}
metadata <- memory$metadata
metadata$memory_measurement <- list(artifact = basename(memory_path),
                                    md5 = unname(tools::md5sum(memory_path)))
merged <- benchmark_merge_instructions(
  memory$results, metadata, readRDS(instructions_path),
  list(artifact = basename(instructions_path),
       md5 = unname(tools::md5sum(instructions_path))))
benchmark_check_baseline(merged$results)
merged$metadata$hosted_artifacts <- args[-seq_len(3L)]
path <- file.path(root, "tests/testthat/_benchmark/baseline.rds")
saveRDS(merged, path, version = 3)
message("Combined baseline written to ", path, "; review provenance before committing.")
