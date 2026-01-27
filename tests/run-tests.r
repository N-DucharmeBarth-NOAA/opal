# Run all tests for biomass extraction functions
# Usage: source this file from R console in the project directory

library(testthat)

# Set project directory (assumes this script is in tests/)
test_dir("tests/testthat", reporter = "progress")
