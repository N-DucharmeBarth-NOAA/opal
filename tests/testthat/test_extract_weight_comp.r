# Tests for extract_ss3_weight_comp and extract_mfcl_weight_comp functions

library(testthat)
library(data.table)
library(magrittr)
library(r4ss)

# Define project paths
proj_dir = this.path::this.proj()
dir_model = file.path(proj_dir, "model-files")
dir_ss3 = file.path(dir_model, "ss3")
dir_mfcl = file.path(dir_model, "mfcl")
dir_helper_fns = file.path(proj_dir, "code", "ss3", "helper-fns")

# Source helper functions
source(file.path(dir_helper_fns, "extract_ss3_weight_comp.r"))
source(file.path(proj_dir, "code", "mfcl", "helper-fns", "extract_mfcl_weight_comp.r"))
source(file.path(dir_helper_fns, "rebin_composition.r"))

context("extract_weight_comp")

# Test SS3 weight composition extractor
test_that("SS3 weight comp extractor produces correct structure", {
  ss3_model_dir = file.path(dir_ss3, "01-bet-base")
  
  # Skip if model directory doesn't exist
  skip_if_not(file.exists(ss3_model_dir), "SS3 model directory not found")
  
  wt_comp = extract_ss3_weight_comp(ss3_model_dir, "01-bet-base", verbose = FALSE)
  
  required_cols = c("id", "Fleet", "Fleet_name", "Used", "Kind", "Sex",
                   "Bin", "Obs", "Exp", "Dev", "effN", "Nsamp_in", "Nsamp_adj")
  
  expect_s3_class(wt_comp, "data.table")
  expect_true(all(required_cols %in% names(wt_comp)))
  expect_true(file.exists(file.path(ss3_model_dir, "comp_size.csv")))
})

test_that("SS3 weight comp has correct column types", {
  ss3_model_dir = file.path(dir_ss3, "01-bet-base")
  skip_if_not(file.exists(ss3_model_dir), "SS3 model directory not found")
  
  wt_comp = extract_ss3_weight_comp(ss3_model_dir, "01-bet-base", verbose = FALSE)
  
  expect_is(wt_comp$id, "character")
  expect_is(wt_comp$Fleet, "integer")
  expect_is(wt_comp$Fleet_name, "character")
  expect_is(wt_comp$Used, "integer")
  expect_is(wt_comp$Kind, "integer")
  expect_is(wt_comp$Sex, "integer")
  expect_is(wt_comp$Bin, "numeric")
  expect_is(wt_comp$Obs, "numeric")
  expect_is(wt_comp$Exp, "numeric")
  expect_is(wt_comp$Dev, "numeric")
  expect_is(wt_comp$effN, "numeric")
  expect_is(wt_comp$Nsamp_in, "numeric")
  expect_is(wt_comp$Nsamp_adj, "numeric")
})

test_that("SS3 weight comp deviations calculated correctly", {
  ss3_model_dir = file.path(dir_ss3, "01-bet-base")
  skip_if_not(file.exists(ss3_model_dir), "SS3 model directory not found")
  
  wt_comp = extract_ss3_weight_comp(ss3_model_dir, "01-bet-base", verbose = FALSE)
  
  # Check that Dev = Obs - Exp for all rows
  expect_equal(wt_comp$Dev, wt_comp$Obs - wt_comp$Exp, tolerance = 1e-10)
})

test_that("SS3 weight comp model ID is correct", {
  ss3_model_dir = file.path(dir_ss3, "01-bet-base")
  skip_if_not(file.exists(ss3_model_dir), "SS3 model directory not found")
  
  wt_comp = extract_ss3_weight_comp(ss3_model_dir, "01-bet-base", verbose = FALSE)
  
  expect_equal(unique(wt_comp$id), "01-bet-base")
})

test_that("MFCL handles missing weight.fit gracefully", {
  mfcl_dir = file.path(dir_mfcl, "v11")
  skip_if_not(file.exists(mfcl_dir), "MFCL directory not found")
  
  # Test with non-existent file
  result = extract_mfcl_weight_comp(
    file.path(mfcl_dir, "nonexistent.fit"),
    file.path(mfcl_dir, "bet.frq"),
    "test",
    output_dir = mfcl_dir,
    verbose = FALSE
  )
  
  # Should return empty data.table with correct structure
  expect_s3_class(result, "data.table")
  expect_equal(nrow(result), 0)
  
  required_cols = c("id", "Fleet", "Fleet_name", "Used", "Kind", "Sex",
                   "Bin", "Obs", "Exp", "Dev", "effN", "Nsamp_in", "Nsamp_adj")
  expect_true(all(required_cols %in% names(result)))
})

test_that("MFCL returns correct structure for empty data", {
  mfcl_dir = file.path(dir_mfcl, "v11")
  skip_if_not(file.exists(mfcl_dir), "MFCL directory not found")
  
  result = extract_mfcl_weight_comp(
    file.path(mfcl_dir, "weight.fit"),
    file.path(mfcl_dir, "bet.frq"),
    "mfcl-v11",
    output_dir = mfcl_dir,
    verbose = FALSE
  )
  
  expect_s3_class(result, "data.table")
  
  required_cols = c("id", "Fleet", "Fleet_name", "Used", "Kind", "Sex",
                   "Bin", "Obs", "Exp", "Dev", "effN", "Nsamp_in", "Nsamp_adj")
  expect_true(all(required_cols %in% names(result)))
})

test_that("Bin harmonization preserves total", {
  ss3_model_dir = file.path(dir_ss3, "01-bet-base")
  skip_if_not(file.exists(ss3_model_dir), "SS3 model directory not found")
  
  # Get original data
  wt_orig = extract_ss3_weight_comp(ss3_model_dir, "01-bet-base",
                                    harmonize_bins = FALSE, verbose = FALSE)
  
  # Skip test if no data
  skip_if(nrow(wt_orig) == 0, "No weight composition data available")
  
  # Define target bins based on data range
  bin_range = range(wt_orig$Bin)
  target_bins = seq(floor(bin_range[1]), ceiling(bin_range[2]) + 5, by = 5)
  
  # Get rebinned data
  wt_rebin = extract_ss3_weight_comp(ss3_model_dir, "01-bet-base",
                                     harmonize_bins = TRUE,
                                     target_bins = target_bins,
                                     verbose = FALSE)
  
  # Check that total observations preserved within each fleet
  for(fleet in unique(wt_orig$Fleet)) {
    total_orig = wt_orig[Fleet == fleet, sum(Obs)]
    total_rebin = wt_rebin[Fleet == fleet, sum(Obs)]
    expect_equal(total_orig, total_rebin, tolerance = 0.001,
                 info = paste("Fleet", fleet, "totals don't match"))
  }
})

test_that("SS3 and MFCL outputs have identical structure", {
  ss3_model_dir = file.path(dir_ss3, "01-bet-base")
  mfcl_dir = file.path(dir_mfcl, "v11")
  
  skip_if_not(file.exists(ss3_model_dir), "SS3 model directory not found")
  skip_if_not(file.exists(mfcl_dir), "MFCL directory not found")
  
  wt_ss3 = extract_ss3_weight_comp(ss3_model_dir, "01-bet-base", verbose = FALSE)
  wt_mfcl = extract_mfcl_weight_comp(
    file.path(mfcl_dir, "weight.fit"),
    file.path(mfcl_dir, "bet.frq"),
    "mfcl-v11",
    output_dir = mfcl_dir,
    verbose = FALSE
  )
  
  # Both should have same column names
  expect_equal(names(wt_ss3), names(wt_mfcl))
  
  # Both should have same column types
  expect_equal(sapply(wt_ss3, class), sapply(wt_mfcl, class))
})

test_that("SS3 weight comp proportions are reasonable", {
  ss3_model_dir = file.path(dir_ss3, "01-bet-base")
  skip_if_not(file.exists(ss3_model_dir), "SS3 model directory not found")
  
  wt_comp = extract_ss3_weight_comp(ss3_model_dir, "01-bet-base", verbose = FALSE)
  
  # Skip if no data
  skip_if(nrow(wt_comp) == 0, "No weight composition data available")
  
  # All proportions should be non-negative
  expect_true(all(wt_comp$Obs >= 0))
  expect_true(all(wt_comp$Exp >= 0))
  
  # Sample sizes should be positive
  expect_true(all(wt_comp$Nsamp_in > 0 | is.na(wt_comp$Nsamp_in)))
})

test_that("SS3 weight comp bins are sorted", {
  ss3_model_dir = file.path(dir_ss3, "01-bet-base")
  skip_if_not(file.exists(ss3_model_dir), "SS3 model directory not found")
  
  wt_comp = extract_ss3_weight_comp(ss3_model_dir, "01-bet-base", verbose = FALSE)
  
  # Skip if no data
  skip_if(nrow(wt_comp) == 0, "No weight composition data available")
  
  # Within each fleet, bins should be sorted
  for(fleet in unique(wt_comp$Fleet)) {
    fleet_bins = wt_comp[Fleet == fleet, Bin]
    expect_equal(fleet_bins, sort(fleet_bins))
  }
})

test_that("Harmonized bins match target bins", {
  ss3_model_dir = file.path(dir_ss3, "01-bet-base")
  skip_if_not(file.exists(ss3_model_dir), "SS3 model directory not found")
  
  # Get original data to determine range
  wt_orig = extract_ss3_weight_comp(ss3_model_dir, "01-bet-base",
                                    harmonize_bins = FALSE, verbose = FALSE)
  
  # Skip if no data
  skip_if(nrow(wt_orig) == 0, "No weight composition data available")
  
  # Define specific target bins
  target_bins = seq(0, 140, by = 2)
  
  wt_rebin = extract_ss3_weight_comp(ss3_model_dir, "01-bet-base",
                                     harmonize_bins = TRUE,
                                     target_bins = target_bins,
                                     verbose = FALSE)
  
  # Bins should be from target_bins (minus the last edge)
  expected_bins = target_bins[-length(target_bins)]
  unique_bins = sort(unique(wt_rebin$Bin))
  
  # Check that rebinned bins are a subset of expected bins
  expect_true(all(unique_bins %in% expected_bins))
})

test_that("CSV files are created in correct locations", {
  ss3_model_dir = file.path(dir_ss3, "01-bet-base")
  skip_if_not(file.exists(ss3_model_dir), "SS3 model directory not found")
  
  # Remove existing file if present
  csv_file = file.path(ss3_model_dir, "comp_size.csv")
  if(file.exists(csv_file)) file.remove(csv_file)
  
  # Extract and verify file is created
  wt_comp = extract_ss3_weight_comp(ss3_model_dir, "01-bet-base", verbose = FALSE)
  expect_true(file.exists(csv_file))
  
  # Verify CSV can be read back
  csv_data = fread(csv_file)
  expect_equal(nrow(csv_data), nrow(wt_comp))
})
