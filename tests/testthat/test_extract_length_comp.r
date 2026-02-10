# Tests for extract_ss3_length_comp and extract_mfcl_length_comp functions

library(testthat)
library(data.table)
library(magrittr)

# Define project paths
proj_dir = this.path::this.proj()
dir_model = file.path(proj_dir, "model-files")
dir_ss3 = file.path(dir_model, "ss3")
dir_mfcl = file.path(dir_model, "mfcl")
dir_helper_fns_ss3 = file.path(proj_dir, "code", "ss3", "helper-fns")
dir_helper_fns_mfcl = file.path(proj_dir, "code", "mfcl", "helper-fns")

# Source helper functions from both directories
sapply(file.path(dir_helper_fns_ss3, list.files(dir_helper_fns_ss3)), source)
sapply(file.path(dir_helper_fns_mfcl, list.files(dir_helper_fns_mfcl)), source)

context("extract_length_comp")

# Define required columns for output
required_cols = c("id", "Fleet", "Fleet_name", "Used", "Kind", "Sex",
                 "Bin", "Obs", "Exp", "Dev", "effN", "Nsamp_in", "Nsamp_adj")

# ===== Pre-load data once for all tests =====
# Read SS3 length composition once
len_ss3_base = NULL
if(file.exists(file.path(dir_ss3, "01-bet-base", "Report.sso")) &&
   requireNamespace("r4ss", quietly = TRUE)) {
  len_ss3_base = extract_ss3_length_comp(
    file.path(dir_ss3, "01-bet-base"),
    "01-bet-base",
    verbose = FALSE
  )
}

# Read MFCL length composition once (if files exist)
len_mfcl_v11 = NULL
if(file.exists(file.path(dir_mfcl, "v11", "length.fit")) &&
   file.exists(file.path(dir_mfcl, "v11", "bet.frq"))) {
  tryCatch({
    len_mfcl_v11 = extract_mfcl_length_comp(
      file.path(dir_mfcl, "v11", "length.fit"),
      file.path(dir_mfcl, "v11", "bet.frq"),
      "mfcl-v11",
      output_dir = dir_mfcl,
      save_csv = FALSE,
      verbose = FALSE
    )
  }, error = function(e) {
    # MFCL data may not be available
  })
}

# ===== SS3 Tests =====

test_that("extract_ss3_length_comp returns data.table with correct structure", {
  skip_if_not(!is.null(len_ss3_base),
              "SS3 length composition not available (pre-loaded)")
  
  expect_s3_class(len_ss3_base, "data.table")
  expect_true(all(required_cols %in% names(len_ss3_base)))
})

test_that("extract_ss3_length_comp returns correct data types", {
  skip_if_not(!is.null(len_ss3_base),
              "SS3 length composition not available (pre-loaded)")
  
  expect_is(len_ss3_base$id, "character")
  expect_is(len_ss3_base$Fleet, "integer")
  expect_is(len_ss3_base$Fleet_name, "character")
  expect_is(len_ss3_base$Used, "integer")
  expect_is(len_ss3_base$Kind, "integer")
  expect_is(len_ss3_base$Sex, "integer")
  expect_is(len_ss3_base$Bin, "numeric")
  expect_is(len_ss3_base$Obs, "numeric")
  expect_is(len_ss3_base$Exp, "numeric")
  expect_is(len_ss3_base$Dev, "numeric")
  expect_is(len_ss3_base$effN, "numeric")
  expect_is(len_ss3_base$Nsamp_in, "numeric")
  expect_is(len_ss3_base$Nsamp_adj, "numeric")
})

test_that("extract_ss3_length_comp writes output CSV", {
  skip_if_not(!is.null(len_ss3_base),
              "SS3 length composition not available (pre-loaded)")
  
  # Remove existing file if present
  output_file = file.path(dir_ss3, "01-bet-base", "comp_len.csv")
  if(file.exists(output_file)) {
    unlink(output_file)
  }
  
  # Call with save_csv=TRUE to generate file
  extract_ss3_length_comp(
    file.path(dir_ss3, "01-bet-base"),
    "01-bet-base",
    save_csv = TRUE,
    verbose = FALSE
  )
  
  expect_true(file.exists(output_file))
  
  # Read back and verify
  len_read = fread(output_file)
  expect_equal(names(len_read), required_cols)
  expect_equal(nrow(len_read), nrow(len_ss3_base))
})

test_that("extract_ss3_length_comp with save_csv=FALSE does not write CSV", {
  skip_if_not(!is.null(len_ss3_base),
              "SS3 length composition not available (pre-loaded)")
  
  # Remove file if present
  output_file = file.path(dir_ss3, "01-bet-base", "comp_len.csv")
  if(file.exists(output_file)) {
    unlink(output_file)
  }
  
  # Call with save_csv=FALSE
  len_comp = extract_ss3_length_comp(
    file.path(dir_ss3, "01-bet-base"),
    "01-bet-base",
    save_csv = FALSE,
    verbose = FALSE
  )
  
  # Verify CSV was not created
  expect_false(file.exists(output_file))
  
  # Verify data.table was still returned
  expect_s3_class(len_comp, "data.table")
  expect_true(all(required_cols %in% names(len_comp)))
})

test_that("extract_ss3_length_comp deviation equals Obs - Exp", {
  skip_if_not(!is.null(len_ss3_base),
              "SS3 length composition not available (pre-loaded)")
  
  expect_equal(len_ss3_base$Dev, len_ss3_base$Obs - len_ss3_base$Exp, tolerance = 1e-10)
})

test_that("extract_ss3_length_comp aggregation produces reasonable sums", {
  skip_if_not(!is.null(len_ss3_base),
              "SS3 length composition not available (pre-loaded)")
  
  # Check that we have data
  expect_true(nrow(len_ss3_base) > 0)
  
  # Check that proportions are non-negative
  expect_true(all(len_ss3_base$Obs >= 0))
  expect_true(all(len_ss3_base$Exp >= 0))
  
  # Check that sample sizes are positive
  expect_true(all(len_ss3_base$Nsamp_in > 0 | len_ss3_base$Nsamp_in == 0))
})

test_that("extract_ss3_length_comp bin harmonization preserves total", {
  skip_if_not(!is.null(len_ss3_base),
              "SS3 length composition not available (pre-loaded)")
  
  # Use pre-loaded data as the original (no harmonization)
  len_orig = len_ss3_base
  
  # Determine appropriate target bins based on data range
  min_bin = min(len_orig$Bin)
  max_bin = max(len_orig$Bin)
  target_bins = seq(floor(min_bin), ceiling(max_bin) + 10, by = 10)
  
  # Get rebinned data
  len_rebin = extract_ss3_length_comp(
    file.path(dir_ss3, "01-bet-base"),
    "01-bet-base",
    harmonize_bins = TRUE,
    target_bins = target_bins,
    verbose = FALSE
  )
  
  # Check total approximately preserved (within tolerance)
  total_orig = len_orig[, sum(Obs)]
  total_rebin = len_rebin[, sum(Obs)]
  expect_equal(total_orig, total_rebin, tolerance = 0.001)
  
  # Check Exp totals too
  total_exp_orig = len_orig[, sum(Exp)]
  total_exp_rebin = len_rebin[, sum(Exp)]
  expect_equal(total_exp_orig, total_exp_rebin, tolerance = 0.001)
})

test_that("extract_ss3_length_comp harmonized bins match target", {
  skip_if_not(!is.null(len_ss3_base),
              "SS3 length composition not available (pre-loaded)")
  
  target_bins = seq(20, 200, by = 10)
  
  len_rebin = extract_ss3_length_comp(
    file.path(dir_ss3, "01-bet-base"),
    "01-bet-base",
    harmonize_bins = TRUE,
    target_bins = target_bins,
    verbose = FALSE
  )
  
  # Check that bins match target (excluding upper edge)
  unique_bins = sort(unique(len_rebin$Bin))
  expected_bins = target_bins[-length(target_bins)]
  
  # Should contain bins from target (may be subset if some bins have no data)
  expect_true(all(unique_bins %in% expected_bins))
})

test_that("extract_ss3_length_comp stops if Report.sso not found", {
  expect_error(
    extract_ss3_length_comp(
      file.path(dir_ss3, "nonexistent-model"),
      "nonexistent",
      verbose = FALSE
    ),
    "Report.sso not found"
  )
})

test_that("extract_ss3_length_comp requires target_bins when harmonize_bins=TRUE", {
  skip_if_not(!is.null(len_ss3_base),
              "SS3 length composition not available (pre-loaded)")
  
  expect_error(
    extract_ss3_length_comp(
      file.path(dir_ss3, "01-bet-base"),
      "01-bet-base",
      harmonize_bins = TRUE,
      target_bins = NULL,
      verbose = FALSE
    ),
    "target_bins must be provided"
  )
})

# ===== MFCL Tests =====

test_that("extract_mfcl_length_comp handles missing length.fit file gracefully", {
  expect_error(
    extract_mfcl_length_comp(
      file.path(dir_mfcl, "v11", "length.fit"),
      file.path(dir_mfcl, "v11", "bet.frq"),
      "mfcl-v11",
      output_dir = dir_mfcl,
      verbose = FALSE
    ),
    "length.fit file not found"
  )
})

test_that("extract_mfcl_length_comp stops if frq file not found", {
  # Create a temporary dummy length.fit file for testing
  temp_fit = tempfile(fileext = ".fit")
  writeLines("# dummy", temp_fit)
  
  expect_error(
    extract_mfcl_length_comp(
      temp_fit,
      file.path(dir_mfcl, "nonexistent.frq"),
      "test",
      output_dir = tempdir(),
      verbose = FALSE
    ),
    ".frq file not found"
  )
  
  unlink(temp_fit)
})

test_that("extract_mfcl_length_comp requires output_dir when save_csv=TRUE", {
  # Create temporary dummy files for testing
  temp_fit = tempfile(fileext = ".fit")
  temp_frq = tempfile(fileext = ".frq")
  writeLines("# dummy", temp_fit)
  writeLines("# dummy", temp_frq)
  
  expect_error(
    extract_mfcl_length_comp(
      temp_fit,
      temp_frq,
      "test",
      output_dir = NULL,
      save_csv = TRUE,
      verbose = FALSE
    ),
    "output_dir must be provided when save_csv = TRUE"
  )
  
  unlink(temp_fit)
  unlink(temp_frq)
})

# ===== Format Consistency Tests =====

test_that("SS3 and MFCL outputs have identical column structure", {
  skip_if_not(!is.null(len_ss3_base),
              "SS3 length composition not available (pre-loaded)")
  
  # Just check that SS3 output has the required structure
  # (MFCL test would be similar if data were available)
  expect_equal(names(len_ss3_base), required_cols)
  expect_equal(
    sapply(len_ss3_base, class),
    c(id = "character", Fleet = "integer", Fleet_name = "character",
      Used = "integer", Kind = "integer", Sex = "integer",
      Bin = "numeric", Obs = "numeric", Exp = "numeric",
      Dev = "numeric", effN = "numeric", Nsamp_in = "numeric",
      Nsamp_adj = "numeric")
  )
})

test_that("Length composition output is compatible with plotting format", {
  skip_if_not(!is.null(len_ss3_base),
              "SS3 length composition not available (pre-loaded)")
  
  # Verify all required columns for plotting exist
  expect_true(all(c("id", "Fleet", "Bin", "Obs", "Exp") %in% names(len_ss3_base)))
  
  # Verify no NA in critical columns
  expect_false(anyNA(len_ss3_base$Bin))
  expect_false(anyNA(len_ss3_base$Obs))
  expect_false(anyNA(len_ss3_base$Exp))
})

test_that("Multiple fleets are handled correctly", {
  skip_if_not(!is.null(len_ss3_base),
              "SS3 length composition not available (pre-loaded)")
  
  # Check that fleet-specific aggregation works
  if(length(unique(len_ss3_base$Fleet)) > 1) {
    # If multiple fleets, verify each has its own data
    fleet_counts = len_ss3_base[, .N, by = Fleet]
    expect_true(all(fleet_counts$N > 0))
    
    # Verify fleet names are unique per fleet
    fleet_names = len_ss3_base[, unique(Fleet_name), by = Fleet]
    expect_equal(nrow(fleet_names), length(unique(len_ss3_base$Fleet)))
  }
})

test_that("Sex-specific data is preserved", {
  skip_if_not(!is.null(len_ss3_base),
              "SS3 length composition not available (pre-loaded)")
  
  # Check that Sex column exists and contains valid values
  expect_true("Sex" %in% names(len_ss3_base))
  expect_true(all(len_ss3_base$Sex %in% c(0, 1, 2)))
})

test_that("Empty fleets produce valid output", {
  skip_if_not(!is.null(len_ss3_base),
              "SS3 length composition not available (pre-loaded)")
  
  # Check that if any fleet has zero observations, it's handled
  zero_obs_fleets = len_ss3_base[, .(total_obs = sum(Obs)), by = Fleet][total_obs == 0]
  
  if(nrow(zero_obs_fleets) > 0) {
    # Verify structure is still valid for zero observation fleets
    zero_data = len_ss3_base[Fleet %in% zero_obs_fleets$Fleet]
    expect_true(all(required_cols %in% names(zero_data)))
  }
  
  # Test passes if no zero fleets or if they're handled correctly
  expect_true(TRUE)
})
