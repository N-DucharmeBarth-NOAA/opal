skip_if_not_selftest <- function() {
  testthat::skip_if_not(
    identical(Sys.getenv("OPAL_SELFTEST"), "true"),
    "Set OPAL_SELFTEST=true to run simulation self-tests"
  )
}

reference_tolerance <- function() {
  as.numeric(Sys.getenv("OPAL_REFERENCE_TOL", "1e-6"))
}