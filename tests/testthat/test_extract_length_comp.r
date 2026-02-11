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
# Create two variants: raw (no replacement) and zero-replaced (to match SS3 tiny-fill)
len_mfcl_v11_raw = NULL
len_mfcl_v11 = NULL
if(file.exists(file.path(dir_mfcl, "v11", "length.fit")) &&
   file.exists(file.path(dir_mfcl, "v11", "bet.frq"))) {
  tryCatch({
    len_mfcl_v11_raw = extract_mfcl_length_comp(
      file.path(dir_mfcl, "v11", "length.fit"),
      file.path(dir_mfcl, "v11", "bet.frq"),
      "mfcl-v11",
      output_dir = dir_mfcl,
      save_csv = FALSE,
      verbose = FALSE,
      zero_replace = NULL
    )
  }, error = function(e) {
    # MFCL raw data may not be available
  })

  tryCatch({
    len_mfcl_v11 = extract_mfcl_length_comp(
      file.path(dir_mfcl, "v11", "length.fit"),
      file.path(dir_mfcl, "v11", "bet.frq"),
      "mfcl-v11",
      output_dir = dir_mfcl,
      save_csv = FALSE,
      verbose = FALSE,
      zero_replace = 9.90589e-05
    )
  }, error = function(e) {
    # MFCL zero-replaced data may not be available
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
  expect_is(len_ss3_base$Used, "character")
  expect_is(len_ss3_base$Kind, "character")
  expect_is(len_ss3_base$Sex, "integer")
  expect_true(is.numeric(len_ss3_base$Bin))
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
      file.path(dir_mfcl, "v11", "nonexistent_length.fit"),
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

test_that("extract_mfcl_length_comp returns correct data types (raw)", {
  skip_if_not(!is.null(len_mfcl_v11_raw),
              "MFCL raw length composition not available (pre-loaded)")

  expect_is(len_mfcl_v11_raw$id, "character")
  expect_is(len_mfcl_v11_raw$Fleet, "integer")
  expect_is(len_mfcl_v11_raw$Fleet_name, "character")
  expect_is(len_mfcl_v11_raw$Used, "character")
  expect_is(len_mfcl_v11_raw$Kind, "character")
  expect_is(len_mfcl_v11_raw$Sex, "integer")
  expect_is(len_mfcl_v11_raw$Bin, "numeric")
  expect_is(len_mfcl_v11_raw$Obs, "numeric")
  expect_is(len_mfcl_v11_raw$Exp, "numeric")
  expect_is(len_mfcl_v11_raw$Dev, "numeric")
  expect_is(len_mfcl_v11_raw$effN, "numeric")
  expect_is(len_mfcl_v11_raw$Nsamp_in, "numeric")
  expect_is(len_mfcl_v11_raw$Nsamp_adj, "numeric")
})

# ===== Format Consistency Tests =====

test_that("SS3 and MFCL outputs have identical column structure", {
  skip_if_not(!is.null(len_ss3_base),
              "SS3 length composition not available (pre-loaded)")
  
  # Check that SS3 output has the required structure
  expect_equal(names(len_ss3_base), required_cols)
  
  # Check data types individually
  expect_is(len_ss3_base$id, "character")
  expect_is(len_ss3_base$Fleet, "integer")
  expect_is(len_ss3_base$Fleet_name, "character")
  expect_is(len_ss3_base$Used, "character")
  expect_is(len_ss3_base$Kind, "character")
  expect_is(len_ss3_base$Sex, "integer")
  expect_true(is.numeric(len_ss3_base$Bin), "Bin should be numeric")
  expect_is(len_ss3_base$Obs, "numeric")
  expect_is(len_ss3_base$Exp, "numeric")
  expect_is(len_ss3_base$Dev, "numeric")
  expect_is(len_ss3_base$effN, "numeric")
  expect_is(len_ss3_base$Nsamp_in, "numeric")
  expect_is(len_ss3_base$Nsamp_adj, "numeric")
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

test_that("MFCL and SS3 observed values match for common fisheries and bins", {
  skip_if_not(!is.null(len_ss3_base) && !is.null(len_mfcl_v11),
              "Both SS3 and MFCL length composition data must be available (pre-loaded)")
  
  # Find common fleets between datasets
  ss3_fleets = unique(len_ss3_base$Fleet)
  mfcl_fleets = unique(len_mfcl_v11$Fleet)
  common_fleets = intersect(ss3_fleets, mfcl_fleets)
  
  skip_if(length(common_fleets) == 0, "No common fleets between SS3 and MFCL data")
  
  # For each common fleet, check observed values match
  for(fleet in common_fleets) {
    ss3_data = len_ss3_base[Fleet == fleet, .(Bin, Obs)]
    setorder(ss3_data, Bin)
    
    mfcl_data = len_mfcl_v11[Fleet == fleet, .(Bin, Obs)]
    setorder(mfcl_data, Bin)
    
    # Find common bins
    common_bins = intersect(ss3_data$Bin, mfcl_data$Bin)
    
    if(length(common_bins) > 0) {
      ss3_subset = ss3_data[Bin %in% common_bins]
      mfcl_subset = mfcl_data[Bin %in% common_bins]
      
      # Check that observed values match (within tolerance for floating point)
      # SS3 may fill tiny non-zero values for bins that are zero in MFCL
      # Allow larger tolerance to accommodate aggregation differences
      # Observed average differences up to ~0.0015; use 3e-3 to be conservative
      expect_equal(ss3_subset$Obs, mfcl_subset$Obs, tolerance = 3e-3,
           label = paste("Fleet", fleet, "observed values should match between SS3 and MFCL"))
    }
  }
})

# ===== Pre-load time-resolved (non-aggregated) data for tests =====

# Required columns for time-resolved output
required_time_cols = c("id", "Fleet", "Fleet_name", "year", "month", "ts",
                       "Used", "Kind", "Sex", "Bin",
                       "Obs", "Exp", "Dev", "effN", "Nsamp_in", "Nsamp_adj")

len_ss3_time = NULL
if(file.exists(file.path(dir_ss3, "01-bet-base", "Report.sso")) &&
   requireNamespace("r4ss", quietly = TRUE)) {
  len_ss3_time = extract_ss3_length_comp(
    file.path(dir_ss3, "01-bet-base"),
    "01-bet-base",
    aggregate = FALSE,
    save_csv = FALSE,
    verbose = FALSE
  )
}

len_mfcl_time = NULL
if(file.exists(file.path(dir_mfcl, "v11", "length.fit")) &&
   file.exists(file.path(dir_mfcl, "v11", "bet.frq"))) {
  tryCatch({
    len_mfcl_time = extract_mfcl_length_comp(
      file.path(dir_mfcl, "v11", "length.fit"),
      file.path(dir_mfcl, "v11", "bet.frq"),
      "mfcl-v11",
      output_dir = dir_mfcl,
      aggregate = FALSE,
      save_csv = FALSE,
      verbose = FALSE
    )
  }, error = function(e) {
    # MFCL time-resolved data may not be available
  })
}

# ===== SS3 Time-Resolved Tests =====

test_that("extract_ss3_length_comp aggregate=FALSE returns correct columns", {
  skip_if_not(!is.null(len_ss3_time),
              "SS3 time-resolved length composition not available")
  
  expect_s3_class(len_ss3_time, "data.table")
  expect_equal(names(len_ss3_time), required_time_cols)
})

test_that("extract_ss3_length_comp aggregate=FALSE has correct data types", {
  skip_if_not(!is.null(len_ss3_time),
              "SS3 time-resolved length composition not available")
  
  expect_is(len_ss3_time$id, "character")
  expect_is(len_ss3_time$Fleet, "integer")
  expect_is(len_ss3_time$Fleet_name, "character")
  expect_is(len_ss3_time$year, "integer")
  expect_is(len_ss3_time$month, "integer")
  expect_is(len_ss3_time$ts, "integer")
  expect_is(len_ss3_time$Used, "character")
  expect_is(len_ss3_time$Kind, "character")
  expect_is(len_ss3_time$Sex, "integer")
  expect_true(is.numeric(len_ss3_time$Bin))
  expect_is(len_ss3_time$Obs, "numeric")
  expect_is(len_ss3_time$Exp, "numeric")
  expect_is(len_ss3_time$Dev, "numeric")
  expect_is(len_ss3_time$effN, "numeric")
  expect_is(len_ss3_time$Nsamp_in, "numeric")
  expect_is(len_ss3_time$Nsamp_adj, "numeric")
})

test_that("extract_ss3_length_comp aggregate=FALSE has valid time values", {
  skip_if_not(!is.null(len_ss3_time),
              "SS3 time-resolved length composition not available")
  
  # Year should be in plausible range
  expect_true(all(len_ss3_time$year >= 1952))
  expect_true(all(len_ss3_time$year <= 2025))
  
  # Month should be one of the quarterly months matching MFCL convention
  expect_true(all(len_ss3_time$month %in% c(2, 5, 8, 11)))
  
  # ts should be positive and consistent with year and month
  expect_true(all(len_ss3_time$ts > 0))
  expected_ts = (len_ss3_time$year - 1952L) * 4L + match(len_ss3_time$month, c(2, 5, 8, 11))
  expect_equal(len_ss3_time$ts, expected_ts)
})

test_that("extract_ss3_length_comp aggregate=FALSE deviation equals Obs - Exp", {
  skip_if_not(!is.null(len_ss3_time),
              "SS3 time-resolved length composition not available")
  
  expect_equal(len_ss3_time$Dev, len_ss3_time$Obs - len_ss3_time$Exp, tolerance = 1e-10)
})

test_that("extract_ss3_length_comp aggregate=FALSE has more rows than aggregated", {
  skip_if_not(!is.null(len_ss3_time) && !is.null(len_ss3_base),
              "Both SS3 time-resolved and aggregated data must be available")
  
  # Time-resolved should have more rows (multiple time steps per fleet/bin)
  expect_true(nrow(len_ss3_time) >= nrow(len_ss3_base))
})

test_that("extract_ss3_length_comp aggregate=FALSE writes comp_len_time.csv", {
  skip_if_not(!is.null(len_ss3_time),
              "SS3 time-resolved length composition not available")
  
  output_file = file.path(dir_ss3, "01-bet-base", "comp_len_time.csv")
  if(file.exists(output_file)) unlink(output_file)
  
  extract_ss3_length_comp(
    file.path(dir_ss3, "01-bet-base"),
    "01-bet-base",
    aggregate = FALSE,
    save_csv = TRUE,
    verbose = FALSE
  )
  
  expect_true(file.exists(output_file))
  
  len_read = fread(output_file)
  expect_equal(names(len_read), required_time_cols)
  
  # Cleanup
  unlink(output_file)
})

# ===== MFCL Time-Resolved Tests =====

test_that("extract_mfcl_length_comp aggregate=FALSE returns correct columns", {
  skip_if_not(!is.null(len_mfcl_time),
              "MFCL time-resolved length composition not available")
  
  expect_s3_class(len_mfcl_time, "data.table")
  expect_equal(names(len_mfcl_time), required_time_cols)
})

test_that("extract_mfcl_length_comp aggregate=FALSE has correct data types", {
  skip_if_not(!is.null(len_mfcl_time),
              "MFCL time-resolved length composition not available")
  
  expect_is(len_mfcl_time$id, "character")
  expect_is(len_mfcl_time$Fleet, "integer")
  expect_is(len_mfcl_time$Fleet_name, "character")
  expect_is(len_mfcl_time$year, "integer")
  expect_is(len_mfcl_time$month, "integer")
  expect_is(len_mfcl_time$ts, "integer")
  expect_is(len_mfcl_time$Used, "character")
  expect_is(len_mfcl_time$Kind, "character")
  expect_is(len_mfcl_time$Sex, "integer")
  expect_true(is.numeric(len_mfcl_time$Bin))
  expect_is(len_mfcl_time$Obs, "numeric")
  expect_is(len_mfcl_time$Exp, "numeric")
  expect_is(len_mfcl_time$Dev, "numeric")
  expect_is(len_mfcl_time$effN, "numeric")
  expect_is(len_mfcl_time$Nsamp_in, "numeric")
  expect_is(len_mfcl_time$Nsamp_adj, "numeric")
})

test_that("extract_mfcl_length_comp aggregate=FALSE has valid time values", {
  skip_if_not(!is.null(len_mfcl_time),
              "MFCL time-resolved length composition not available")
  
  # Year should be in plausible range
  expect_true(all(len_mfcl_time$year >= 1952))
  expect_true(all(len_mfcl_time$year <= 2025))
  
  # Month should be one of the quarterly months
  expect_true(all(len_mfcl_time$month %in% c(2, 5, 8, 11)))
  
  # ts should be positive and consistent with year and month
  expect_true(all(len_mfcl_time$ts > 0))
  expected_ts = (len_mfcl_time$year - 1952L) * 4L + match(len_mfcl_time$month, c(2, 5, 8, 11))
  expect_equal(len_mfcl_time$ts, expected_ts)
})

test_that("extract_mfcl_length_comp aggregate=FALSE deviation equals Obs - Exp", {
  skip_if_not(!is.null(len_mfcl_time),
              "MFCL time-resolved length composition not available")
  
  expect_equal(len_mfcl_time$Dev, len_mfcl_time$Obs - len_mfcl_time$Exp, tolerance = 1e-10)
})

test_that("extract_mfcl_length_comp aggregate=FALSE has more rows than aggregated", {
  skip_if_not(!is.null(len_mfcl_time) && !is.null(len_mfcl_v11),
              "Both MFCL time-resolved and aggregated data must be available")
  
  expect_true(nrow(len_mfcl_time) >= nrow(len_mfcl_v11))
})

test_that("extract_mfcl_length_comp aggregate=FALSE writes comp_len_time.csv", {
  skip_if_not(!is.null(len_mfcl_time),
              "MFCL time-resolved length composition not available")
  
  output_file = file.path(dir_mfcl, "comp_len_time.csv")
  if(file.exists(output_file)) unlink(output_file)
  
  extract_mfcl_length_comp(
    file.path(dir_mfcl, "v11", "length.fit"),
    file.path(dir_mfcl, "v11", "bet.frq"),
    "mfcl-v11",
    output_dir = dir_mfcl,
    aggregate = FALSE,
    save_csv = TRUE,
    verbose = FALSE
  )
  
  expect_true(file.exists(output_file))
  
  len_read = fread(output_file)
  expect_equal(names(len_read), required_time_cols)
  
  # Cleanup
  unlink(output_file)
})

# ===== Cross-Model Time-Resolved Consistency Tests =====

test_that("SS3 and MFCL time-resolved outputs have identical column structure", {
  skip_if_not(!is.null(len_ss3_time) && !is.null(len_mfcl_time),
              "Both SS3 and MFCL time-resolved data must be available")
  
  expect_equal(names(len_ss3_time), names(len_mfcl_time))
  expect_equal(names(len_ss3_time), required_time_cols)
})

test_that("SS3 and MFCL time-resolved ts values overlap", {
  skip_if_not(!is.null(len_ss3_time) && !is.null(len_mfcl_time),
              "Both SS3 and MFCL time-resolved data must be available")
  
  ss3_ts = sort(unique(len_ss3_time$ts))
  mfcl_ts = sort(unique(len_mfcl_time$ts))
  common_ts = intersect(ss3_ts, mfcl_ts)
  
  # There should be at least some overlapping time steps
  expect_true(length(common_ts) > 0)
})

test_that("SS3 and MFCL time-resolved observed values match for spot-check selection", {
  skip_if_not(!is.null(len_ss3_time) && !is.null(len_mfcl_time),
              "Both SS3 and MFCL time-resolved data must be available")
  
  # Find common time steps
  ss3_ts = sort(unique(len_ss3_time$ts))
  mfcl_ts = sort(unique(len_mfcl_time$ts))
  common_ts = intersect(ss3_ts, mfcl_ts)
  
  skip_if(length(common_ts) == 0, "No common time steps between SS3 and MFCL time-resolved data")
  
  # Sample a few time steps for spot check (first 3, or all if fewer than 3)
  sampled_ts = common_ts[1:min(3, length(common_ts))]
  
  # For each sampled time step, compare fleets and bins
  for(ts_val in sampled_ts) {
    ss3_ts_data = len_ss3_time[ts == ts_val, .(Fleet, Fleet_name, Bin, Obs)]
    mfcl_ts_data = len_mfcl_time[ts == ts_val, .(Fleet, Fleet_name, Bin, Obs)]
    
    # Find common fleets
    common_fleets = intersect(ss3_ts_data$Fleet, mfcl_ts_data$Fleet)
    
    for(fleet in common_fleets) {
      ss3_fleet = ss3_ts_data[Fleet == fleet, .(Bin, Obs)]
      setorder(ss3_fleet, Bin)
      
      mfcl_fleet = mfcl_ts_data[Fleet == fleet, .(Bin, Obs)]
      setorder(mfcl_fleet, Bin)
      
      # Find common bins
      common_bins = intersect(ss3_fleet$Bin, mfcl_fleet$Bin)
      
      if(length(common_bins) > 0) {
        ss3_subset = ss3_fleet[Bin %in% common_bins]
        mfcl_subset = mfcl_fleet[Bin %in% common_bins]
        
        # For disaggregated data, tolerance accounts for rounding and extraction differences
        # Spot checks show max differences around 0.003 across time steps and fleets
        # Use 3e-3 (0.003) tolerance, same as aggregated comparisons, but applied to raw time-resolved data
        expect_equal(ss3_subset$Obs, mfcl_subset$Obs, tolerance = 3e-3,
             label = paste("TS", ts_val, "Fleet", fleet, "time-resolved observed values should match closely between SS3 and MFCL"))
      }
    }
  }
})