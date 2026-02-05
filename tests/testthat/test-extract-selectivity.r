# Tests for selectivity extraction functions

library(testthat)
library(data.table)
library(magrittr)
library(r4ss)
library(FLR4MFCL)

# Define project paths
proj_dir = this.path::this.proj()
dir_model = file.path(proj_dir, "model-files")
dir_ss3 = file.path(dir_model, "ss3")
dir_mfcl = file.path(dir_model, "mfcl")
dir_helper_fns_ss3 = file.path(proj_dir, "code", "ss3", "helper-fns")
dir_helper_fns_mfcl = file.path(proj_dir, "code", "mfcl", "helper-fns")

# Source helper functions
sapply(file.path(dir_helper_fns_ss3, list.files(dir_helper_fns_ss3)), source)
sapply(file.path(dir_helper_fns_mfcl, list.files(dir_helper_fns_mfcl)), source)
context("extract_ss3_selectivity")

# Test data paths
ss3_dir = file.path(dir_ss3, "01-bet-base")
mfcl_rep = file.path(dir_mfcl, "v11", "plot-10.par.rep")
mfcl_par = file.path(dir_mfcl, "v11", "10.par")
mfcl_dir = file.path(dir_mfcl, "v11")

# Compute result once and reuse across tests to avoid repeatedly parsing Report.sso
selex_ss3 = extract_ss3_selectivity(ss3_dir, "01-bet-base", verbose = FALSE)

test_that("SS3 selectivity extractor returns data.table with correct structure", {
  expect_s3_class(selex_ss3, "data.table")
  expect_named(selex_ss3, c("id", "Fleet", "Fleet_name", "Yr", "Sex", "variable", "value"))
})

test_that("SS3 selectivity extractor returns correct data types", {
  expect_is(selex_ss3$id, "character")
  expect_is(selex_ss3$Fleet, "integer")
  expect_is(selex_ss3$Fleet_name, "character")
  expect_is(selex_ss3$Yr, "integer")
  expect_is(selex_ss3$Sex, "integer")
  expect_is(selex_ss3$variable, "numeric")
  expect_is(selex_ss3$value, "numeric")
})

test_that("SS3 selectivity values are between 0 and 1", {
  expect_true(all(selex_ss3$value >= 0 & selex_ss3$value <= 1), 
              info = paste("Found selectivity values outside [0,1]:", 
                           paste(selex_ss3[value < 0 | value > 1]$value, collapse = ", ")))
})

test_that("SS3 selectivity length values are positive", {
  expect_true(all(selex_ss3$variable > 0),
              info = paste("Found non-positive length values"))
})

test_that("SS3 selectivity fleet numbers are positive integers", {
  expect_true(all(selex_ss3$Fleet > 0))
  expect_is(selex_ss3$Fleet, "integer")
})

test_that("SS3 selectivity extractor writes CSV file", {
  expect_true(file.exists(file.path(ss3_dir, "selex_l.csv")))
})

test_that("SS3 selectivity CSV file is readable and matches output", {
  csv_data = fread(file.path(ss3_dir, "selex_l.csv"))
  expect_equal(nrow(csv_data), nrow(selex_ss3))
  expect_equal(names(csv_data), names(selex_ss3))
})

test_that("SS3 selectivity has no NA values in required columns", {
  required_cols = c("id", "Fleet", "Fleet_name", "Yr", "Sex", "variable", "value")
  expect_false(anyNA(selex_ss3[, ..required_cols]), 
               info = paste("Found NA values in columns:", 
                            paste(colnames(selex_ss3)[colSums(is.na(selex_ss3)) > 0], collapse = ", ")))
})

test_that("SS3 model ID matches input", {
  expect_equal(unique(selex_ss3$id), "01-bet-base")
})

test_that("SS3 selectivity extractor with write_csv=FALSE does not write file", {
  # Use a temporary directory to test
  temp_dir = file.path(ss3_dir, "temp_test")
  dir.create(temp_dir, showWarnings = FALSE)
  
  # Copy Report.sso to temp dir for testing
  file.copy(file.path(ss3_dir, "Report.sso"), temp_dir, overwrite = TRUE)
  
  # Extract without writing CSV
  selex_no_csv = extract_ss3_selectivity(temp_dir, "test-no-csv", write_csv = FALSE, verbose = FALSE)
  
  # Check that CSV was not created
  expect_false(file.exists(file.path(temp_dir, "selex_l.csv")))
  
  # Check that data.table was still returned
  expect_s3_class(selex_no_csv, "data.table")
  expect_true(nrow(selex_no_csv) > 0)
  
  # Cleanup
  unlink(temp_dir, recursive = TRUE)
})

context("extract_mfcl_selectivity")

# Compute MFCL result once
selex_mfcl = extract_mfcl_selectivity(
  rep_file = mfcl_rep,
  par_file = mfcl_par,
  model_id = "v11",
  output_dir = mfcl_dir,
  verbose = FALSE
)

test_that("MFCL selectivity extractor returns data.table with correct structure", {
  expect_s3_class(selex_mfcl, "data.table")
  expect_named(selex_mfcl, c("id", "Fleet", "Fleet_name", "Yr", "Sex", "variable", "value"))
})

test_that("MFCL selectivity extractor returns correct data types", {
  expect_is(selex_mfcl$id, "character")
  expect_is(selex_mfcl$Fleet, "integer")
  expect_is(selex_mfcl$Fleet_name, "character")
  expect_is(selex_mfcl$Yr, "integer")
  expect_is(selex_mfcl$Sex, "integer")
  expect_is(selex_mfcl$variable, "numeric")
  expect_is(selex_mfcl$value, "numeric")
})

test_that("MFCL selectivity values are between 0 and 1", {
  expect_true(all(selex_mfcl$value >= 0 & selex_mfcl$value <= 1), 
              info = paste("Found selectivity values outside [0,1]:", 
                           paste(selex_mfcl[value < 0 | value > 1]$value, collapse = ", ")))
})

test_that("MFCL selectivity length values are positive", {
  expect_true(all(selex_mfcl$variable > 0),
              info = paste("Found non-positive length values"))
})

test_that("MFCL selectivity fleet numbers are positive integers", {
  expect_true(all(selex_mfcl$Fleet > 0))
  expect_is(selex_mfcl$Fleet, "integer")
})

test_that("MFCL selectivity extractor writes CSV file", {
  expect_true(file.exists(file.path(mfcl_dir, "selex_l.csv")))
})

test_that("MFCL selectivity CSV file is readable and matches output", {
  csv_data = fread(file.path(mfcl_dir, "selex_l.csv"))
  expect_equal(nrow(csv_data), nrow(selex_mfcl))
  expect_equal(names(csv_data), names(selex_mfcl))
})

test_that("MFCL selectivity has no NA values in required columns", {
  required_cols = c("id", "Fleet", "Fleet_name", "Yr", "Sex", "variable", "value")
  expect_false(anyNA(selex_mfcl[, ..required_cols]), 
               info = paste("Found NA values in columns:", 
                            paste(colnames(selex_mfcl)[colSums(is.na(selex_mfcl)) > 0], collapse = ", ")))
})

test_that("MFCL model ID matches input", {
  expect_equal(unique(selex_mfcl$id), "v11")
})

test_that("MFCL selectivity Sex is 0 (aggregated)", {
  expect_true(all(selex_mfcl$Sex == 0))
})

test_that("MFCL selectivity extractor with write_csv=FALSE does not write file", {
  # Extract without writing CSV
  selex_no_csv = extract_mfcl_selectivity(
    rep_file = mfcl_rep,
    par_file = mfcl_par,
    model_id = "test-no-csv",
    write_csv = FALSE,
    verbose = FALSE
  )
  
  # Check that data.table was still returned
  expect_s3_class(selex_no_csv, "data.table")
  expect_true(nrow(selex_no_csv) > 0)
  expect_equal(names(selex_no_csv), c("id", "Fleet", "Fleet_name", "Yr", "Sex", "variable", "value"))
})

test_that("MFCL selectivity extractor requires output_dir when write_csv=TRUE", {
  # Should error when output_dir is NULL and write_csv is TRUE
  expect_error(
    extract_mfcl_selectivity(
      rep_file = mfcl_rep,
      par_file = mfcl_par,
      model_id = "test-error",
      output_dir = NULL,
      write_csv = TRUE,
      verbose = FALSE
    ),
    "output_dir must be specified when write_csv = TRUE"
  )
})

context("extractor_compatibility")

test_that("SS3 and MFCL outputs have identical column names", {
  expect_equal(names(selex_ss3), names(selex_mfcl))
})

test_that("SS3 and MFCL outputs have identical column types", {
  expect_equal(sapply(selex_ss3, class), sapply(selex_mfcl, class))
})

test_that("SS3 and MFCL outputs have same column order", {
  expect_equal(names(selex_ss3), c("id", "Fleet", "Fleet_name", "Yr", "Sex", "variable", "value"))
  expect_equal(names(selex_mfcl), c("id", "Fleet", "Fleet_name", "Yr", "Sex", "variable", "value"))
})

test_that("Both extractors produce data in expected value ranges", {
  # Both should have selectivity 0-1
  expect_true(all(selex_ss3$value >= 0 & selex_ss3$value <= 1))
  expect_true(all(selex_mfcl$value >= 0 & selex_mfcl$value <= 1))
  
  # Both should have positive lengths
  expect_true(all(selex_ss3$variable > 0))
  expect_true(all(selex_mfcl$variable > 0))
})

