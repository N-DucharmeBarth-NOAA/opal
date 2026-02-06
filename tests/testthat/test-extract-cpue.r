# Tests for extract_ss3_cpue and extract_mfcl_cpue functions

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
context("extract_ss3_cpue")

# Test SS3 extractor with CSV saving (default behavior)
ss3_dir = file.path(dir_ss3, "01-bet-base")
cpue_ss3 = extract_ss3_cpue(ss3_dir, "01-bet-base", verbose = FALSE)

test_that("SS3 CPUE extractor produces correct structure", {
  expect_s3_class(cpue_ss3, "data.table")
  expect_true(all(c("id", "Fleet", "Fleet_name", "Time", "Obs", "Exp", "SE", "Dev", "Use") %in% names(cpue_ss3)))
})

test_that("SS3 CPUE extractor produces correct data types", {
  expect_type(cpue_ss3$id, "character")
  expect_type(cpue_ss3$Fleet, "double")
  expect_type(cpue_ss3$Fleet_name, "character")
  expect_type(cpue_ss3$Time, "double")
  expect_type(cpue_ss3$Obs, "double")
  expect_type(cpue_ss3$Exp, "double")
  expect_type(cpue_ss3$SE, "double")
  expect_type(cpue_ss3$Dev, "double")
  expect_type(cpue_ss3$Use, "integer")
})

test_that("SS3 CPUE extractor produces non-negative SE", {
  expect_true(all(cpue_ss3$SE >= 0), 
              info = paste("Found SE values < 0"))
})

test_that("SS3 CPUE extractor Use flag is binary", {
  expect_true(all(cpue_ss3$Use %in% c(0, 1)), 
              info = paste("Found Use values not in {0, 1}"))
})

test_that("SS3 CPUE extractor writes CSV file by default", {
  expect_true(file.exists(file.path(ss3_dir, "cpue.csv")))
})

test_that("SS3 CPUE extractor can skip CSV saving", {
  # Remove CSV if it exists
  csv_path = file.path(ss3_dir, "cpue_test.csv")
  if(file.exists(csv_path)) file.remove(csv_path)
  
  # Extract without saving CSV
  cpue_no_csv = extract_ss3_cpue(ss3_dir, "01-bet-base", save_csv = FALSE, verbose = FALSE)
  
  # Verify data.table is still returned
  expect_s3_class(cpue_no_csv, "data.table")
  expect_equal(nrow(cpue_no_csv), nrow(cpue_ss3))
})

test_that("SS3 CPUE extractor model ID matches", {
  expect_equal(unique(cpue_ss3$id), "01-bet-base")
})

context("extract_mfcl_cpue")

# Test MFCL extractor
mfcl_frq = file.path(dir_mfcl, "v11", "bet.frq")
mfcl_rep = file.path(dir_mfcl, "v11", "plot-10.par.rep")
mfcl_dir = file.path(dir_mfcl, "v11")
cpue_mfcl = extract_mfcl_cpue(mfcl_frq, mfcl_rep, "mfcl-v11", output_dir = mfcl_dir, verbose = FALSE)

test_that("MFCL CPUE extractor produces correct structure", {
  expect_s3_class(cpue_mfcl, "data.table")
  expect_true(all(c("id", "Fleet", "Fleet_name", "Time", "Obs", "Exp", "SE", "Dev", "Use") %in% names(cpue_mfcl)))
})

test_that("MFCL CPUE extractor produces correct data types", {
  expect_type(cpue_mfcl$id, "character")
  expect_type(cpue_mfcl$Fleet, "double")
  expect_type(cpue_mfcl$Fleet_name, "character")
  expect_type(cpue_mfcl$Time, "double")
  expect_type(cpue_mfcl$Obs, "double")
  expect_type(cpue_mfcl$Exp, "double")
  expect_type(cpue_mfcl$SE, "double")
  expect_type(cpue_mfcl$Dev, "double")
  expect_type(cpue_mfcl$Use, "integer")
})

test_that("MFCL CPUE extractor produces positive values after exp() conversion", {
  expect_true(all(cpue_mfcl$Obs > 0), 
              info = paste("Found Obs values <= 0 - log space conversion may have failed"))
  expect_true(all(cpue_mfcl$Exp > 0), 
              info = paste("Found Exp values <= 0 - log space conversion may have failed"))
})

test_that("MFCL CPUE extractor produces non-negative SE", {
  expect_true(all(cpue_mfcl$SE >= 0), 
              info = paste("Found SE values < 0"))
})

test_that("MFCL CPUE extractor calculates Dev correctly", {
  expect_equal(cpue_mfcl$Dev, cpue_mfcl$Obs - cpue_mfcl$Exp, 
               tolerance = 1e-10)
})

test_that("MFCL CPUE extractor Use flag is binary", {
  expect_true(all(cpue_mfcl$Use %in% c(0, 1)), 
              info = paste("Found Use values not in {0, 1}"))
})

test_that("MFCL CPUE extractor writes CSV file by default", {
  expect_true(file.exists(file.path(mfcl_dir, "cpue.csv")))
})

test_that("MFCL CPUE extractor can skip CSV saving", {
  # Extract without saving CSV
  cpue_no_csv = extract_mfcl_cpue(mfcl_frq, mfcl_rep, "mfcl-v11", 
                                   save_csv = FALSE, verbose = FALSE)
  
  # Verify data.table is still returned
  expect_s3_class(cpue_no_csv, "data.table")
  expect_equal(nrow(cpue_no_csv), nrow(cpue_mfcl))
})

test_that("MFCL CPUE extractor requires output_dir when save_csv = TRUE", {
  expect_error(
    extract_mfcl_cpue(mfcl_frq, mfcl_rep, "mfcl-v11", save_csv = TRUE, verbose = FALSE),
    "output_dir must be provided when save_csv = TRUE"
  )
})

test_that("MFCL CPUE extractor model ID matches", {
  expect_equal(unique(cpue_mfcl$id), "mfcl-v11")
})

context("SS3 and MFCL compatibility")

test_that("SS3 and MFCL outputs have identical column names", {
  expect_equal(names(cpue_ss3), names(cpue_mfcl))
})

test_that("SS3 and MFCL outputs have compatible column types", {
  # Note: Fleet may be integer in SS3 and numeric in MFCL, both are acceptable
  ss3_types = sapply(cpue_ss3, class)
  mfcl_types = sapply(cpue_mfcl, class)
  
  # Check key columns have compatible types (numeric or integer are compatible)
  expect_true(ss3_types["id"] == mfcl_types["id"])
  expect_true(ss3_types["Fleet_name"] == mfcl_types["Fleet_name"])
  expect_true(all(c(ss3_types[c("Time", "Obs", "Exp", "SE", "Dev")] == "numeric")))
  expect_true(all(c(mfcl_types[c("Time", "Obs", "Exp", "SE", "Dev")] == "numeric")))
})


context("Error handling")

test_that("SS3 extractor stops if Report.sso not found", {
  expect_error(extract_ss3_cpue(file.path(dir_ss3, "nonexistent-model"), "test"))
})

test_that("MFCL extractor stops if rep file not found", {
  expect_error(extract_mfcl_cpue(mfcl_frq, "nonexistent.rep", "test", output_dir = mfcl_dir))
})

test_that("MFCL extractor stops if frq file not found", {
  expect_error(extract_mfcl_cpue("nonexistent.frq", mfcl_rep, "test", output_dir = mfcl_dir))
})
