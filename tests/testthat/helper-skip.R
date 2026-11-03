skip_if_not_selftest <- function() {
  testthat::skip_if_not(
    identical(Sys.getenv("OPAL_SELFTEST"), "true"),
    "Set OPAL_SELFTEST=true to run simulation self-tests"
  )
}

reference_tolerance <- function() {
  as.numeric(Sys.getenv("OPAL_REFERENCE_TOL", "1e-6"))
}

skip_if_not_benchmark <- function() {
  testthat::skip_if_not(
    identical(Sys.getenv("OPAL_BENCHMARK"), "true"),
    "Set OPAL_BENCHMARK=true to run performance benchmarks"
  )
}

benchmark_tolerance <- function() {
  value <- suppressWarnings(as.numeric(Sys.getenv("OPAL_BENCHMARK_TOL", "0.10")))
  if (length(value) != 1L || !is.finite(value) || value < 0) {
    stop("OPAL_BENCHMARK_TOL must be a finite non-negative number")
  }
  value
}