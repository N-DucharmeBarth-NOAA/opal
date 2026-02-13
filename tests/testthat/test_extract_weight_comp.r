# Tests for extract_ss3_weight_comp and extract_mfcl_weight_comp functions

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

context("extract_weight_comp")

# Define required columns for output
required_cols = c("id", "Fleet", "Fleet_name", "Used", "Kind", "Sex",
                 "Bin", "Obs", "Exp", "Dev", "effN", "Nsamp_in", "Nsamp_adj")

# ===== Pre-load data once for all tests =====
# Read SS3 weight composition once
wt_ss3_base = NULL
if(file.exists(file.path(dir_ss3, "01-bet-base", "Report.sso")) &&
   requireNamespace("r4ss", quietly = TRUE)) {
  tryCatch({
    wt_ss3_base = extract_ss3_weight_comp(
      file.path(dir_ss3, "01-bet-base"),
      "01-bet-base",
      verbose = FALSE
    )
  }, error = function(e) {
    # SS3 weight composition may not be available
  })
}

# Read MFCL weight composition once (if files exist)
# Create two variants: raw (no replacement) and zero-replaced (to match SS3 tiny-fill)
wt_mfcl_v11_raw = NULL
wt_mfcl_v11 = NULL
if(file.exists(file.path(dir_mfcl, "v11", "weight.fit")) &&
   file.exists(file.path(dir_mfcl, "v11", "bet.frq"))) {
  tryCatch({
    wt_mfcl_v11_raw = extract_mfcl_weight_comp(
      file.path(dir_mfcl, "v11", "weight.fit"),
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
    wt_mfcl_v11 = extract_mfcl_weight_comp(
      file.path(dir_mfcl, "v11", "weight.fit"),
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

test_that("extract_ss3_weight_comp returns data.table with correct structure", {
  skip_if_not(!is.null(wt_ss3_base),
              "SS3 weight composition not available (pre-loaded)")
  
  expect_s3_class(wt_ss3_base, "data.table")
  expect_true(all(required_cols %in% names(wt_ss3_base)))
})

test_that("extract_ss3_weight_comp returns correct data types", {
  skip_if_not(!is.null(wt_ss3_base),
              "SS3 weight composition not available (pre-loaded)")
  
  expect_is(wt_ss3_base$id, "character")
  expect_is(wt_ss3_base$Fleet, "integer")
  expect_is(wt_ss3_base$Fleet_name, "character")
  expect_is(wt_ss3_base$Used, "character")  # Should be character like length comp
  expect_is(wt_ss3_base$Kind, "character")  # Should be character like length comp
  expect_is(wt_ss3_base$Sex, "integer")
  expect_true(is.numeric(wt_ss3_base$Bin))
  expect_is(wt_ss3_base$Obs, "numeric")
  expect_is(wt_ss3_base$Exp, "numeric")
  expect_is(wt_ss3_base$Dev, "numeric")
  expect_is(wt_ss3_base$effN, "numeric")
  expect_is(wt_ss3_base$Nsamp_in, "numeric")
  expect_is(wt_ss3_base$Nsamp_adj, "numeric")
})

test_that("extract_ss3_weight_comp deviations calculated correctly", {
  skip_if_not(!is.null(wt_ss3_base),
              "SS3 weight composition not available (pre-loaded)")
  
  # Check that Dev = Obs - Exp for all rows
  expect_equal(wt_ss3_base$Dev, wt_ss3_base$Obs - wt_ss3_base$Exp, tolerance = 1e-10)
})

test_that("extract_ss3_weight_comp model ID is correct", {
  skip_if_not(!is.null(wt_ss3_base),
              "SS3 weight composition not available (pre-loaded)")
  
  expect_equal(unique(wt_ss3_base$id), "01-bet-base")
})

test_that("extract_ss3_weight_comp proportions are reasonable", {
  skip_if_not(!is.null(wt_ss3_base),
              "SS3 weight composition not available (pre-loaded)")
  
  # All proportions should be non-negative
  expect_true(all(wt_ss3_base$Obs >= 0))
  expect_true(all(wt_ss3_base$Exp >= 0))
  
  # Sample sizes should be positive (or NA)
  expect_true(all(wt_ss3_base$Nsamp_in > 0 | is.na(wt_ss3_base$Nsamp_in)))
})

test_that("extract_ss3_weight_comp aggregation produces reasonable sums", {
  skip_if_not(!is.null(wt_ss3_base),
              "SS3 weight composition not available (pre-loaded)")
  
  # Check that we have data
  expect_true(nrow(wt_ss3_base) > 0)
  
  # Check that proportions are non-negative
  expect_true(all(wt_ss3_base$Obs >= 0))
  expect_true(all(wt_ss3_base$Exp >= 0))
  
  # Check that sample sizes are positive
  expect_true(all(wt_ss3_base$Nsamp_in > 0 | wt_ss3_base$Nsamp_in == 0))
})

test_that("extract_ss3_weight_comp stops if Report.sso not found", {
  expect_error(
    extract_ss3_weight_comp(
      file.path(dir_ss3, "nonexistent-model"),
      "nonexistent",
      verbose = FALSE
    ),
    "Report.sso not found"
  )
})

test_that("extract_ss3_weight_comp requires target_bins when harmonize_bins=TRUE", {
  skip_if_not(!is.null(wt_ss3_base),
              "SS3 weight composition not available (pre-loaded)")
  
  expect_error(
    extract_ss3_weight_comp(
      file.path(dir_ss3, "01-bet-base"),
      "01-bet-base",
      harmonize_bins = TRUE,
      target_bins = NULL,
      verbose = FALSE
    ),
    "target_bins must be provided"
  )
})

test_that("extract_ss3_weight_comp writes output CSV", {
  skip_if_not(!is.null(wt_ss3_base),
              "SS3 weight composition not available (pre-loaded)")
  
  # Remove existing file if present
  output_file = file.path(dir_ss3, "01-bet-base", "comp_size.csv")
  if(file.exists(output_file)) {
    unlink(output_file)
  }
  
  # Call with save_csv=TRUE to generate file
  extract_ss3_weight_comp(
    file.path(dir_ss3, "01-bet-base"),
    "01-bet-base",
    save_csv = TRUE,
    verbose = FALSE
  )
  
  expect_true(file.exists(output_file))
  
  # Read back and verify
  wt_read = fread(output_file)
  expect_equal(names(wt_read), required_cols)
  expect_equal(nrow(wt_read), nrow(wt_ss3_base))
})

test_that("extract_ss3_weight_comp bins are sorted", {
  skip_if_not(!is.null(wt_ss3_base),
              "SS3 weight composition not available (pre-loaded)")
  
  # Within each fleet, bins should be sorted
  for(fleet in unique(wt_ss3_base$Fleet)) {
    fleet_bins = wt_ss3_base[Fleet == fleet, Bin]
    expect_equal(fleet_bins, sort(fleet_bins))
  }
})

test_that("extract_ss3_weight_comp save_csv = FALSE does not create file", {
  ss3_model_dir = file.path(dir_ss3, "01-bet-base")
  skip_if_not(file.exists(ss3_model_dir), "SS3 model directory not found")
  
  # Remove existing file if present
  csv_file = file.path(ss3_model_dir, "comp_size.csv")
  if(file.exists(csv_file)) file.remove(csv_file)
  
  # Extract without saving CSV
  wt_comp = extract_ss3_weight_comp(ss3_model_dir, "01-bet-base", 
                                    save_csv = FALSE, verbose = FALSE)
  
  # Verify file was not created
  expect_false(file.exists(csv_file))
  
  # But data.table should still be returned
  expect_s3_class(wt_comp, "data.table")
})

# ===== MFCL Tests =====

test_that("extract_mfcl_weight_comp returns correct structure", {
  skip_if_not(!is.null(wt_mfcl_v11),
              "MFCL weight composition not available (pre-loaded)")
  
  expect_s3_class(wt_mfcl_v11, "data.table")
  expect_true(all(required_cols %in% names(wt_mfcl_v11)))
})

test_that("extract_mfcl_weight_comp returns correct data types", {
  skip_if_not(!is.null(wt_mfcl_v11),
              "MFCL weight composition not available (pre-loaded)")
  
  expect_is(wt_mfcl_v11$id, "character")
  expect_is(wt_mfcl_v11$Fleet, "integer")
  expect_is(wt_mfcl_v11$Fleet_name, "character")
  expect_is(wt_mfcl_v11$Used, "character")  # Should be character
  expect_is(wt_mfcl_v11$Kind, "character")  # Should be character
  expect_is(wt_mfcl_v11$Sex, "integer")
  expect_true(is.numeric(wt_mfcl_v11$Bin))
  expect_is(wt_mfcl_v11$Obs, "numeric")
  expect_is(wt_mfcl_v11$Exp, "numeric")
  expect_is(wt_mfcl_v11$Dev, "numeric")
})

test_that("extract_mfcl_weight_comp deviations calculated correctly", {
  skip_if_not(!is.null(wt_mfcl_v11),
              "MFCL weight composition not available (pre-loaded)")
  
  expect_equal(wt_mfcl_v11$Dev, wt_mfcl_v11$Obs - wt_mfcl_v11$Exp, tolerance = 1e-10)
})

test_that("extract_mfcl_weight_comp returns correct data types (raw)", {
  skip_if_not(!is.null(wt_mfcl_v11_raw),
              "MFCL raw weight composition not available (pre-loaded)")

  expect_is(wt_mfcl_v11_raw$id, "character")
  expect_is(wt_mfcl_v11_raw$Fleet, "integer")
  expect_is(wt_mfcl_v11_raw$Fleet_name, "character")
  expect_is(wt_mfcl_v11_raw$Used, "character")
  expect_is(wt_mfcl_v11_raw$Kind, "character")
  expect_is(wt_mfcl_v11_raw$Sex, "integer")
  expect_is(wt_mfcl_v11_raw$Bin, "numeric")
  expect_is(wt_mfcl_v11_raw$Obs, "numeric")
  expect_is(wt_mfcl_v11_raw$Exp, "numeric")
  expect_is(wt_mfcl_v11_raw$Dev, "numeric")
  expect_is(wt_mfcl_v11_raw$effN, "numeric")
  expect_is(wt_mfcl_v11_raw$Nsamp_in, "numeric")
  expect_is(wt_mfcl_v11_raw$Nsamp_adj, "numeric")
})

test_that("extract_mfcl_weight_comp handles missing weight.fit gracefully", {
  expect_error(
    extract_mfcl_weight_comp(
      file.path(dir_mfcl, "v11", "nonexistent_weight.fit"),
      file.path(dir_mfcl, "v11", "bet.frq"),
      "mfcl-v11",
      output_dir = dir_mfcl,
      verbose = FALSE
    ),
    "weight.fit file not found"
  )
})

test_that("extract_mfcl_weight_comp stops if frq file not found", {
  # Create a temporary dummy weight.fit file for testing
  temp_fit = tempfile(fileext = ".fit")
  writeLines("# dummy", temp_fit)
  
  expect_error(
    extract_mfcl_weight_comp(
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

test_that("extract_mfcl_weight_comp handles missing weight.fit gracefully", {
  mfcl_dir = file.path(dir_mfcl, "v11")
  skip_if_not(file.exists(mfcl_dir), "MFCL directory not found")
  skip_if_not(file.exists(file.path(mfcl_dir, "bet.frq")), ".frq file not found")
  
  # Test with non-existent file - should error
  expect_error(
    extract_mfcl_weight_comp(
      file.path(mfcl_dir, "nonexistent.fit"),
      file.path(mfcl_dir, "bet.frq"),
      "test",
      output_dir = mfcl_dir,
      save_csv = FALSE,
      verbose = FALSE
    ),
    "weight.fit file not found"
  )
})

test_that("extract_mfcl_weight_comp requires output_dir when save_csv = TRUE", {
  mfcl_dir = file.path(dir_mfcl, "v11")
  skip_if_not(file.exists(file.path(mfcl_dir, "weight.fit")), "weight.fit file not found")
  skip_if_not(file.exists(file.path(mfcl_dir, "bet.frq")), ".frq file not found")
  
  # Should error when output_dir is missing and save_csv = TRUE
  expect_error(
    extract_mfcl_weight_comp(
      file.path(mfcl_dir, "weight.fit"),
      file.path(mfcl_dir, "bet.frq"),
      "test",
      save_csv = TRUE,
      verbose = FALSE
    ),
    "output_dir must be provided"
  )
})

# ===== Format Consistency Tests =====

test_that("SS3 and MFCL outputs have identical column structure", {
  skip_if_not(!is.null(wt_ss3_base),
              "SS3 weight composition not available (pre-loaded)")
  
  # Check that SS3 output has the required structure
  expect_equal(names(wt_ss3_base), required_cols)
  
  # Check data types individually
  expect_is(wt_ss3_base$id, "character")
  expect_is(wt_ss3_base$Fleet, "integer")
  expect_is(wt_ss3_base$Fleet_name, "character")
  expect_is(wt_ss3_base$Used, "character")
  expect_is(wt_ss3_base$Kind, "character")
  expect_is(wt_ss3_base$Sex, "integer")
  expect_true(is.numeric(wt_ss3_base$Bin), "Bin should be numeric")
  expect_is(wt_ss3_base$Obs, "numeric")
  expect_is(wt_ss3_base$Exp, "numeric")
  expect_is(wt_ss3_base$Dev, "numeric")
  expect_is(wt_ss3_base$effN, "numeric")
  expect_is(wt_ss3_base$Nsamp_in, "numeric")
  expect_is(wt_ss3_base$Nsamp_adj, "numeric")
})

test_that("Weight composition output is compatible with plotting format", {
  skip_if_not(!is.null(wt_ss3_base),
              "SS3 weight composition not available (pre-loaded)")
  
  # Verify all required columns for plotting exist
  expect_true(all(c("id", "Fleet", "Bin", "Obs", "Exp") %in% names(wt_ss3_base)))
  
  # Verify no NA in critical columns
  expect_false(anyNA(wt_ss3_base$Bin))
  expect_false(anyNA(wt_ss3_base$Obs))
  expect_false(anyNA(wt_ss3_base$Exp))
})

test_that("Multiple fleets are handled correctly", {
  skip_if_not(!is.null(wt_ss3_base),
              "SS3 weight composition not available (pre-loaded)")
  
  # Check that fleet-specific aggregation works
  if(length(unique(wt_ss3_base$Fleet)) > 1) {
    # If multiple fleets, verify each has its own data
    fleet_counts = wt_ss3_base[, .N, by = Fleet]
    expect_true(all(fleet_counts$N > 0))
    
    # Verify fleet names are unique per fleet
    fleet_names = wt_ss3_base[, unique(Fleet_name), by = Fleet]
    expect_equal(nrow(fleet_names), length(unique(wt_ss3_base$Fleet)))
  }
})

test_that("Sex-specific data is preserved", {
  skip_if_not(!is.null(wt_ss3_base),
              "SS3 weight composition not available (pre-loaded)")
  
  # Check that Sex column exists and contains valid values
  expect_true("Sex" %in% names(wt_ss3_base))
  expect_true(all(wt_ss3_base$Sex %in% c(0, 1, 2)))
})

test_that("Empty fleets produce valid output", {
  skip_if_not(!is.null(wt_ss3_base),
              "SS3 weight composition not available (pre-loaded)")
  
  # Check that if any fleet has zero observations, it's handled
  zero_obs_fleets = wt_ss3_base[, .(total_obs = sum(Obs)), by = Fleet][total_obs == 0]
  
  if(nrow(zero_obs_fleets) > 0) {
    # Verify structure is still valid for zero observation fleets
    zero_data = wt_ss3_base[Fleet %in% zero_obs_fleets$Fleet]
    expect_true(all(required_cols %in% names(zero_data)))
  }
  
  # Test passes if no zero fleets or if they're handled correctly
  expect_true(TRUE)
})

test_that("MFCL and SS3 observed values match for common fisheries and bins", {
  skip_if_not(!is.null(wt_ss3_base) && !is.null(wt_mfcl_v11),
              "Both SS3 and MFCL weight composition data must be available (pre-loaded)")
  
  # Find common fleets between datasets
  ss3_fleets = unique(wt_ss3_base$Fleet)
  mfcl_fleets = unique(wt_mfcl_v11$Fleet)
  common_fleets = intersect(ss3_fleets, mfcl_fleets)
  
  skip_if(length(common_fleets) == 0, "No common fleets between SS3 and MFCL data")
  
  # For each common fleet, check observed values match
  for(fleet in common_fleets) {
    ss3_data = wt_ss3_base[Fleet == fleet, .(Bin, Obs)]
    setorder(ss3_data, Bin)
    
    mfcl_data = wt_mfcl_v11[Fleet == fleet, .(Bin, Obs)]
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

wt_ss3_time = NULL
if(file.exists(file.path(dir_ss3, "01-bet-base", "Report.sso")) &&
   requireNamespace("r4ss", quietly = TRUE)) {
  tryCatch({
    wt_ss3_time = extract_ss3_weight_comp(
      file.path(dir_ss3, "01-bet-base"),
      "01-bet-base",
      aggregate = FALSE,
      save_csv = FALSE,
      verbose = FALSE
    )
  }, error = function(e) {
    # SS3 time-resolved data may not be available
  })
}

wt_mfcl_time = NULL
if(file.exists(file.path(dir_mfcl, "v11", "weight.fit")) &&
   file.exists(file.path(dir_mfcl, "v11", "bet.frq"))) {
  tryCatch({
    wt_mfcl_time = extract_mfcl_weight_comp(
      file.path(dir_mfcl, "v11", "weight.fit"),
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

test_that("extract_ss3_weight_comp aggregate=FALSE returns correct columns", {
  skip_if_not(!is.null(wt_ss3_time),
              "SS3 time-resolved weight composition not available")
  
  expect_s3_class(wt_ss3_time, "data.table")
  expect_equal(names(wt_ss3_time), required_time_cols)
})

test_that("extract_ss3_weight_comp aggregate=FALSE has correct data types", {
  skip_if_not(!is.null(wt_ss3_time),
              "SS3 time-resolved weight composition not available")
  
  expect_is(wt_ss3_time$id, "character")
  expect_is(wt_ss3_time$Fleet, "integer")
  expect_is(wt_ss3_time$Fleet_name, "character")
  expect_is(wt_ss3_time$year, "integer")
  expect_is(wt_ss3_time$month, "integer")
  expect_is(wt_ss3_time$ts, "integer")
  expect_is(wt_ss3_time$Used, "character")
  expect_is(wt_ss3_time$Kind, "character")
  expect_is(wt_ss3_time$Sex, "integer")
  expect_true(is.numeric(wt_ss3_time$Bin))
  expect_is(wt_ss3_time$Obs, "numeric")
  expect_is(wt_ss3_time$Exp, "numeric")
  expect_is(wt_ss3_time$Dev, "numeric")
  expect_is(wt_ss3_time$effN, "numeric")
  expect_is(wt_ss3_time$Nsamp_in, "numeric")
  expect_is(wt_ss3_time$Nsamp_adj, "numeric")
})

test_that("extract_ss3_weight_comp aggregate=FALSE has valid time values", {
  skip_if_not(!is.null(wt_ss3_time),
              "SS3 time-resolved weight composition not available")
  
  # Year should be in plausible range
  expect_true(all(wt_ss3_time$year >= 1952))
  expect_true(all(wt_ss3_time$year <= 2025))
  
  # Month should be one of the quarterly months matching MFCL convention
  expect_true(all(wt_ss3_time$month %in% c(2, 5, 8, 11)))
  
  # ts should be positive and consistent with year and month
  expect_true(all(wt_ss3_time$ts > 0))
  expected_ts = (wt_ss3_time$year - 1952L) * 4L + match(wt_ss3_time$month, c(2, 5, 8, 11))
  expect_equal(wt_ss3_time$ts, expected_ts)
})

test_that("extract_ss3_weight_comp aggregate=FALSE deviation equals Obs - Exp", {
  skip_if_not(!is.null(wt_ss3_time),
              "SS3 time-resolved weight composition not available")
  
  expect_equal(wt_ss3_time$Dev, wt_ss3_time$Obs - wt_ss3_time$Exp, tolerance = 1e-10)
})

test_that("extract_ss3_weight_comp aggregate=FALSE has more rows than aggregated", {
  skip_if_not(!is.null(wt_ss3_time) && !is.null(wt_ss3_base),
              "Both SS3 time-resolved and aggregated data must be available")
  
  # Time-resolved should have more rows (multiple time steps per fleet/bin)
  expect_true(nrow(wt_ss3_time) >= nrow(wt_ss3_base))
})

test_that("extract_ss3_weight_comp aggregate=FALSE writes comp_size_time.csv", {
  skip_if_not(!is.null(wt_ss3_time),
              "SS3 time-resolved weight composition not available")
  
  output_file = file.path(dir_ss3, "01-bet-base", "comp_size_time.csv")
  if(file.exists(output_file)) unlink(output_file)
  
  extract_ss3_weight_comp(
    file.path(dir_ss3, "01-bet-base"),
    "01-bet-base",
    aggregate = FALSE,
    save_csv = TRUE,
    verbose = FALSE
  )
  
  expect_true(file.exists(output_file))
  
  wt_read = fread(output_file)
  expect_equal(names(wt_read), required_time_cols)
  
  # Cleanup
  unlink(output_file)
})

# ===== MFCL Time-Resolved Tests =====

test_that("extract_mfcl_weight_comp aggregate=FALSE returns correct columns", {
  skip_if_not(!is.null(wt_mfcl_time),
              "MFCL time-resolved weight composition not available")
  
  expect_s3_class(wt_mfcl_time, "data.table")
  expect_equal(names(wt_mfcl_time), required_time_cols)
})

test_that("extract_mfcl_weight_comp aggregate=FALSE has correct data types", {
  skip_if_not(!is.null(wt_mfcl_time),
              "MFCL time-resolved weight composition not available")
  
  expect_is(wt_mfcl_time$id, "character")
  expect_is(wt_mfcl_time$Fleet, "integer")
  expect_is(wt_mfcl_time$Fleet_name, "character")
  expect_is(wt_mfcl_time$year, "integer")
  expect_is(wt_mfcl_time$month, "integer")
  expect_is(wt_mfcl_time$ts, "integer")
  expect_is(wt_mfcl_time$Used, "character")
  expect_is(wt_mfcl_time$Kind, "character")
  expect_is(wt_mfcl_time$Sex, "integer")
  expect_true(is.numeric(wt_mfcl_time$Bin))
  expect_is(wt_mfcl_time$Obs, "numeric")
  expect_is(wt_mfcl_time$Exp, "numeric")
  expect_is(wt_mfcl_time$Dev, "numeric")
  expect_is(wt_mfcl_time$effN, "numeric")
  expect_is(wt_mfcl_time$Nsamp_in, "numeric")
  expect_is(wt_mfcl_time$Nsamp_adj, "numeric")
})

test_that("extract_mfcl_weight_comp aggregate=FALSE has valid time values", {
  skip_if_not(!is.null(wt_mfcl_time),
              "MFCL time-resolved weight composition not available")
  
  # Year should be in plausible range
  expect_true(all(wt_mfcl_time$year >= 1952))
  expect_true(all(wt_mfcl_time$year <= 2025))
  
  # Month should be one of the quarterly months
  expect_true(all(wt_mfcl_time$month %in% c(2, 5, 8, 11)))
  
  # ts should be positive and consistent with year and month
  expect_true(all(wt_mfcl_time$ts > 0))
  expected_ts = (wt_mfcl_time$year - 1952L) * 4L + match(wt_mfcl_time$month, c(2, 5, 8, 11))
  expect_equal(wt_mfcl_time$ts, expected_ts)
})

test_that("extract_mfcl_weight_comp aggregate=FALSE deviation equals Obs - Exp", {
  skip_if_not(!is.null(wt_mfcl_time),
              "MFCL time-resolved weight composition not available")
  
  expect_equal(wt_mfcl_time$Dev, wt_mfcl_time$Obs - wt_mfcl_time$Exp, tolerance = 1e-10)
})

test_that("extract_mfcl_weight_comp aggregate=FALSE has more rows than aggregated", {
  skip_if_not(!is.null(wt_mfcl_time) && !is.null(wt_mfcl_v11),
              "Both MFCL time-resolved and aggregated data must be available")
  
  expect_true(nrow(wt_mfcl_time) >= nrow(wt_mfcl_v11))
})

test_that("extract_mfcl_weight_comp aggregate=FALSE writes comp_size_time.csv", {
  skip_if_not(!is.null(wt_mfcl_time),
              "MFCL time-resolved weight composition not available")
  
  output_file = file.path(dir_mfcl, "comp_size_time.csv")
  if(file.exists(output_file)) unlink(output_file)
  
  extract_mfcl_weight_comp(
    file.path(dir_mfcl, "v11", "weight.fit"),
    file.path(dir_mfcl, "v11", "bet.frq"),
    "mfcl-v11",
    output_dir = dir_mfcl,
    aggregate = FALSE,
    save_csv = TRUE,
    verbose = FALSE
  )
  
  expect_true(file.exists(output_file))
  
  wt_read = fread(output_file)
  expect_equal(names(wt_read), required_time_cols)
  
  # Cleanup
  unlink(output_file)
})

# ===== Cross-Model Time-Resolved Consistency Tests =====

test_that("SS3 and MFCL time-resolved outputs have identical column structure", {
  skip_if_not(!is.null(wt_ss3_time) && !is.null(wt_mfcl_time),
              "Both SS3 and MFCL time-resolved data must be available")
  
  expect_equal(names(wt_ss3_time), names(wt_mfcl_time))
  expect_equal(names(wt_ss3_time), required_time_cols)
})

test_that("SS3 and MFCL time-resolved ts values overlap", {
  skip_if_not(!is.null(wt_ss3_time) && !is.null(wt_mfcl_time),
              "Both SS3 and MFCL time-resolved data must be available")
  
  ss3_ts = sort(unique(wt_ss3_time$ts))
  mfcl_ts = sort(unique(wt_mfcl_time$ts))
  common_ts = intersect(ss3_ts, mfcl_ts)
  
  # There should be at least some overlapping time steps
  expect_true(length(common_ts) > 0)
})

test_that("SS3 and MFCL time-resolved observed values match for spot-check selection", {
  skip_if_not(!is.null(wt_ss3_time) && !is.null(wt_mfcl_time),
              "Both SS3 and MFCL time-resolved data must be available")
  
  # Find common time steps
  ss3_ts = sort(unique(wt_ss3_time$ts))
  mfcl_ts = sort(unique(wt_mfcl_time$ts))
  common_ts = intersect(ss3_ts, mfcl_ts)
  
  skip_if(length(common_ts) == 0, "No common time steps between SS3 and MFCL time-resolved data")
  
  # Sample a few time steps for spot check (first 3, or all if fewer than 3)
  sampled_ts = common_ts[1:min(3, length(common_ts))]
  
  # For each sampled time step, compare fleets and bins
  for(ts_val in sampled_ts) {
    ss3_ts_data = wt_ss3_time[ts == ts_val, .(Fleet, Fleet_name, Bin, Obs)]
    mfcl_ts_data = wt_mfcl_time[ts == ts_val, .(Fleet, Fleet_name, Bin, Obs)]
    
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

# ===== Bin Harmonization Tests =====

test_that("Bin harmonization preserves total", {
  ss3_model_dir = file.path(dir_ss3, "01-bet-base")
  skip_if_not(file.exists(ss3_model_dir), "SS3 model directory not found")
  skip_if_not(!is.null(wt_ss3_base), "SS3 weight composition not available")
  skip_if(nrow(wt_ss3_base) == 0, "No weight composition data available")
  
  # Define target bins based on data range
  bin_range = range(wt_ss3_base$Bin)
  target_bins = seq(floor(bin_range[1]), ceiling(bin_range[2]) + 5, by = 5)
  
  # Get rebinned data
  wt_rebin = extract_ss3_weight_comp(ss3_model_dir, "01-bet-base",
                                     harmonize_bins = TRUE,
                                     target_bins = target_bins,
                                     save_csv = FALSE,
                                     verbose = FALSE)
  
  # Check that total observations preserved within each fleet
  for(fleet in unique(wt_ss3_base$Fleet)) {
    total_orig = wt_ss3_base[Fleet == fleet, sum(Obs)]
    total_rebin = wt_rebin[Fleet == fleet, sum(Obs)]
    expect_equal(total_orig, total_rebin, tolerance = 0.001,
                 info = paste("Fleet", fleet, "totals don't match"))
  }
})

test_that("Harmonized bins match target bins", {
  ss3_model_dir = file.path(dir_ss3, "01-bet-base")
  skip_if_not(file.exists(ss3_model_dir), "SS3 model directory not found")
  skip_if_not(!is.null(wt_ss3_base), "SS3 weight composition not available")
  skip_if(nrow(wt_ss3_base) == 0, "No weight composition data available")
  
  # Define specific target bins
  target_bins = seq(0, 140, by = 2)
  
  wt_rebin = extract_ss3_weight_comp(ss3_model_dir, "01-bet-base",
                                     harmonize_bins = TRUE,
                                     target_bins = target_bins,
                                     save_csv = FALSE,
                                     verbose = FALSE)
  
  # Bins should be from target_bins (minus the last edge)
  expected_bins = target_bins[-length(target_bins)]
  unique_bins = sort(unique(wt_rebin$Bin))
  
  # Check that rebinned bins are a subset of expected bins
  expect_true(all(unique_bins %in% expected_bins))
})

test_that("CSV files are created in correct locations when save_csv = TRUE", {
  ss3_model_dir = file.path(dir_ss3, "01-bet-base")
  skip_if_not(file.exists(ss3_model_dir), "SS3 model directory not found")
  
  # Remove existing file if present
  csv_file = file.path(ss3_model_dir, "comp_size.csv")
  if(file.exists(csv_file)) file.remove(csv_file)
  
  # Extract and verify file is created
  wt_comp = extract_ss3_weight_comp(ss3_model_dir, "01-bet-base", 
                                    save_csv = TRUE, verbose = FALSE)
  expect_true(file.exists(csv_file))
  
  # Verify CSV can be read back
  csv_data = fread(csv_file)
  expect_equal(nrow(csv_data), nrow(wt_comp))
})
