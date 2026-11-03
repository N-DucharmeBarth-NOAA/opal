# Explicit regeneration only; never sourced by tests.
# Keep complete reference fits, but replay at most 50 fn and 50 gr calls spread
# across the path. --reuse-fits verifies and samples existing fits without nlminb.
root <- Sys.getenv("OPAL_BENCHMARK_ROOT", getwd())
Sys.setenv(OPAL_BENCHMARK_ROOT = normalizePath(root, mustWork = TRUE))
source(file.path(root, "tests/testthat/helper-benchmark.R"))
scratch <- tempdir()
models <- commandArgs(trailingOnly = TRUE)
reuse <- "--reuse-fits" %in% models
models <- setdiff(models, "--reuse-fits")
if (!length(models)) models <- c("opaka", "bet", "io_yft")
if (any(!models %in% c("opaka", "bet", "io_yft"))) stop("Unknown trajectory model.")
out <- file.path(root, "tests/testthat/_benchmark/trajectories")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
for (model in models) {
  message(if (reuse) "Validating and sampling stored reference fit: " else
    "Recording two-pass reference fit: ", model)
  trajectory <- callr::r(function(root, model, reuse) {
    source(file.path(root, "tests/testthat/helper-benchmark.R"))
    benchmark_load()
    inputs <- benchmark_inputs(model)
    if (reuse) {
      trajectory <- benchmark_read_trajectory(model,
        file.path(root, "tests/testthat/_benchmark/trajectories"), inputs)
    } else {
      obj <- benchmark_objective(inputs)
      fit <- benchmark_fit(obj, inputs$parameters, record = TRUE)
      trajectory <- c(list(model = model, inputs = inputs, metadata = benchmark_metadata()), fit)
    }
    trajectory <- benchmark_sample_trajectory(trajectory)
    trajectory$replay_metadata <- benchmark_metadata()
    benchmark_validate_trajectory(trajectory)
    trajectory
  }, args = list(root = root, model = model, reuse = reuse), libpath = .libPaths(),
  env = c(OPAL_BENCHMARK_ROOT = root, TMPDIR = scratch, OMP_NUM_THREADS = "1",
          OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1"),
  system_profile = FALSE, user_profile = FALSE, show = TRUE)
  saveRDS(trajectory, file.path(out, paste0(model, ".rds")), version = 3)
  message(model, ": ", length(trajectory$events), " calls; ",
          trajectory$counts[["iterations"]], " iterations")
}
