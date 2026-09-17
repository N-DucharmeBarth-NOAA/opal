#!/usr/bin/env Rscript

root <- normalizePath(".")
pkgload::load_all(root, quiet = TRUE)
source(file.path(root, "tests", "testthat", "helper-opaka.R"))
source(file.path(root, "tests", "testthat", "helper-selftest.R"))

n_sim <- as.integer(Sys.getenv("OPAL_SELFTEST_NSIM", "100"))
output_directory <- Sys.getenv(
  "OPAL_SELFTEST_OUT", file.path(tempdir(), "opal-selftest-calibration")
)
dir.create(output_directory, showWarnings = FALSE, recursive = TRUE)

for (start in c("truth", "default")) {
  started <- Sys.time()
  results <- run_selftest(n_sim, start = start)
  elapsed_seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  summary <- summarize_selftest(results)
  summary$start <- start
  summary$elapsed_seconds <- elapsed_seconds

  saveRDS(results, file.path(output_directory, sprintf("selftest-opaka-%s.rds", start)))
  utils::write.csv(
    summary,
    file.path(output_directory, sprintf("selftest-opaka-%s-summary.csv", start)),
    row.names = FALSE
  )
  print(summary)
}