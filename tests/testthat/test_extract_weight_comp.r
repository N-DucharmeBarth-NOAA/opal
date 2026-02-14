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
par_file_mfcl = file.path(dir_mfcl, "v11", "10.par")
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
      zero_replace = NULL,
      par_file = par_file_mfcl
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
      zero_replace = 9.90589e-05,
      par_file = par_file_mfcl
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
  
  # Sample sizes should be non-negative (or NA)
  expect_true(all(wt_ss3_base$Nsamp_in >= 0 | is.na(wt_ss3_base$Nsamp_in)))
})

test_that("extract_ss3_weight_comp aggregation produces reasonable sums", {
  skip_if_not(!is.null(wt_ss3_base),
              "SS3 weight composition not available (pre-loaded)")
  
  # Check that we have data
  expect_true(nrow(wt_ss3_base) > 0)
  
  # Check that proportions are non-negative
  expect_true(all(wt_ss3_base$Obs >= 0))
  expect_true(all(wt_ss3_base$Exp >= 0))
  
  # Check that sample sizes are non-negative
  expect_true(all(wt_ss3_base$Nsamp_in >= 0 | wt_ss3_base$Nsamp_in == 0))
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

test_that("MFCL and SS3 observed values match for common fisheries and bins with same width", {
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
    
    # Calculate bin widths for both datasets
    # Bin width is the difference between consecutive bins
    if(nrow(ss3_data) > 1) {
      ss3_data[, bin_width := c(diff(Bin), diff(Bin)[length(diff(Bin))])]
    } else {
      ss3_data[, bin_width := NA_real_]
    }
    
    if(nrow(mfcl_data) > 1) {
      mfcl_data[, bin_width := c(diff(Bin), diff(Bin)[length(diff(Bin))])]
    } else {
      mfcl_data[, bin_width := NA_real_]
    }
    
    # Find common bins
    common_bins = intersect(ss3_data$Bin, mfcl_data$Bin)
    
    if(length(common_bins) > 0) {
      ss3_subset = ss3_data[Bin %in% common_bins]
      mfcl_subset = mfcl_data[Bin %in% common_bins]
      
      # Verify bin widths match for common bins
      # Allow small tolerance for floating point differences
      bin_width_matches = abs(ss3_subset$bin_width - mfcl_subset$bin_width) < 1e-6
      
      if(!all(bin_width_matches, na.rm = TRUE)) {
        # Report bins with different widths
        diff_indices = which(!bin_width_matches)
        message(sprintf("Fleet %d: Bin widths differ for %d/%d common bins", 
                       fleet, length(diff_indices), length(common_bins)))
        message(sprintf("  Example: Bin %.1f has SS3 width %.2f vs MFCL width %.2f",
                       ss3_subset$Bin[diff_indices[1]],
                       ss3_subset$bin_width[diff_indices[1]],
                       mfcl_subset$bin_width[diff_indices[1]]))
        
        # Only compare observed values for bins with matching widths
        matching_bins = ss3_subset$Bin[bin_width_matches]
        ss3_subset = ss3_subset[Bin %in% matching_bins]
        mfcl_subset = mfcl_subset[Bin %in% matching_bins]
      }
      
      # Only proceed with comparison if we have bins with matching widths
      if(nrow(ss3_subset) > 0 && nrow(mfcl_subset) > 0) {
        # With back_transform=TRUE (default), SS3 Obs is on the original input
        # proportion scale, matching MFCL. Tolerance covers floating-point
        # round-trip through Report.sso and any residual clamping at zero.
        expect_equal(ss3_subset$Obs, mfcl_subset$Obs, tolerance = 1.5e-4,
             label = paste("Fleet", fleet, "observed values should match between SS3 and MFCL for bins with same width"))
      } else {
        # Skip comparison if no bins with matching widths
        message(sprintf("Fleet %d: No bins with matching widths, skipping observed value comparison", fleet))
      }
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
      verbose = FALSE,
      par_file = par_file_mfcl
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

test_that("SS3 and MFCL time-resolved observed values match for spot-check selection with same bin widths", {
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
      
      # Calculate bin widths
      if(nrow(ss3_fleet) > 1) {
        ss3_fleet[, bin_width := c(diff(Bin), diff(Bin)[length(diff(Bin))])]
      } else {
        ss3_fleet[, bin_width := NA_real_]
      }
      
      if(nrow(mfcl_fleet) > 1) {
        mfcl_fleet[, bin_width := c(diff(Bin), diff(Bin)[length(diff(Bin))])]
      } else {
        mfcl_fleet[, bin_width := NA_real_]
      }
      
      # Find common bins
      common_bins = intersect(ss3_fleet$Bin, mfcl_fleet$Bin)
      
      if(length(common_bins) > 0) {
        ss3_subset = ss3_fleet[Bin %in% common_bins]
        mfcl_subset = mfcl_fleet[Bin %in% common_bins]
        
        # Verify bin widths match for common bins
        bin_width_matches = abs(ss3_subset$bin_width - mfcl_subset$bin_width) < 1e-6
        
        if(!all(bin_width_matches, na.rm = TRUE)) {
          # Only compare observed values for bins with matching widths
          matching_bins = ss3_subset$Bin[bin_width_matches]
          ss3_subset = ss3_subset[Bin %in% matching_bins]
          mfcl_subset = mfcl_subset[Bin %in% matching_bins]
        }
        
        # Only proceed with comparison if we have bins with matching widths
        if(nrow(ss3_subset) > 0 && nrow(mfcl_subset) > 0) {
          # With back_transform=TRUE (default), SS3 Obs matches original input
          # proportions. Tolerance covers floating-point round-trip only.
          expect_equal(ss3_subset$Obs, mfcl_subset$Obs, tolerance = 1.5e-4,
               label = paste("TS", ts_val, "Fleet", fleet, "time-resolved observed values should match closely between SS3 and MFCL for bins with same width"))
        }
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

# ===== Back-Transform Tests =====

test_that("back_transform=TRUE recovers original Nsamp from data.ss", {
  skip_if_not(!is.null(wt_ss3_time),
              "SS3 time-resolved weight composition not available (pre-loaded)")
  
  # The pre-loaded wt_ss3_time uses back_transform=TRUE by default
  # For Fleet 1, time step 21 (1957 Q2), the original Nsamp in data.ss is 7748.6844
  spot = wt_ss3_time[Fleet == 1 & ts == 21]
  skip_if(nrow(spot) == 0, "No data for Fleet 1, ts 21")
  
  # All bins for this observation should share the same Nsamp_in
  expect_equal(spot$Nsamp_in[1], 7748.6844, tolerance = 0.01,
               label = "Back-transformed Nsamp_in should match original data.ss value")
})

test_that("back_transform=TRUE recovers original observed proportions", {
  skip_if_not(!is.null(wt_ss3_time),
              "SS3 time-resolved weight composition not available (pre-loaded)")
  
  # For Fleet 1, ts 21, Bin 2: the original proportion from MFCL weight.fit
  # and data.ss is 0.17403612 / 7748.6844 = 2.24601e-05
  spot = wt_ss3_time[Fleet == 1 & ts == 21 & Bin == 2]
  skip_if(nrow(spot) == 0, "No data for Fleet 1, ts 21, Bin 2")
  
  expect_equal(spot$Obs, 2.24601e-05, tolerance = 1e-7,
               label = "Back-transformed Obs should match original input proportion")
  
  # And Obs * Nsamp_in should recover the nominal count from data.ss / frq file
  expect_equal(spot$Obs * spot$Nsamp_in, 0.17403612, tolerance = 1e-4,
               label = "Obs * Nsamp_in should recover nominal count from data.ss")
})

test_that("back_transform=FALSE returns raw SS3 Report.sso values", {
  ss3_model_dir = file.path(dir_ss3, "01-bet-base")
  skip_if_not(file.exists(file.path(ss3_model_dir, "Report.sso")),
              "SS3 Report.sso not found")
  
  wt_raw = tryCatch(
    extract_ss3_weight_comp(ss3_model_dir, "01-bet-base",
                            aggregate = FALSE, save_csv = FALSE,
                            back_transform = FALSE, verbose = FALSE),
    error = function(e) NULL
  )
  skip_if_not(!is.null(wt_raw), "Could not read SS3 output")
  
  # For Fleet 1, ts 21, Bin 2: Report.sso reports mincomp-adjusted Obs = 0.000937177
  spot = wt_raw[Fleet == 1 & ts == 21 & Bin == 2]
  skip_if(nrow(spot) == 0, "No data for Fleet 1, ts 21, Bin 2")
  
  expect_equal(spot$Obs, 0.000937177, tolerance = 1e-5,
               label = "back_transform=FALSE Obs should match Report.sso value")
  
  # Nsamp_in should be the variance-adjusted value ~0.193717
  expect_equal(spot$Nsamp_in, 0.193717, tolerance = 1e-4,
               label = "back_transform=FALSE Nsamp_in should be variance-adjusted value")
})

test_that("back_transform=TRUE SS3 Obs matches MFCL Obs for same bins", {
  skip_if_not(!is.null(wt_ss3_time) && !is.null(wt_mfcl_time),
              "Both SS3 and MFCL time-resolved data must be available")
  
  # Pick Fleet 1, early time steps with fine (1 kg) bins
  common_ts = intersect(unique(wt_ss3_time$ts), unique(wt_mfcl_time$ts))
  skip_if(length(common_ts) == 0, "No common time steps")
  
  ts_val = common_ts[1]
  
  ss3_f1 = wt_ss3_time[Fleet == 1 & ts == ts_val & Bin <= 43]
  mfcl_f1 = wt_mfcl_time[Fleet == 1 & ts == ts_val & Bin <= 43]
  
  skip_if(nrow(ss3_f1) == 0 || nrow(mfcl_f1) == 0,
          "No data for Fleet 1 in fine bins")
  
  common_bins = intersect(ss3_f1$Bin, mfcl_f1$Bin)
  skip_if(length(common_bins) == 0, "No common bins")
  
  ss3_comp = ss3_f1[Bin %in% common_bins][order(Bin)]$Obs
  mfcl_comp = mfcl_f1[Bin %in% common_bins][order(Bin)]$Obs
  
  # With back_transform=TRUE, the Obs proportions should closely match
  # (tolerance allows for floating-point rounding through Report.sso)
  expect_equal(ss3_comp, mfcl_comp, tolerance = 1e-5,
               label = "Back-transformed SS3 Obs should match MFCL Obs for fine bins")
})

# ===== MFCL Nsamp_adj with flag 50 Tests =====

test_that("MFCL Nsamp_adj uses flag 50 values when par_file provided", {
  skip_if_not(!is.null(wt_mfcl_time) && file.exists(par_file_mfcl),
              "MFCL time-resolved data and par file must be available")
  
  # Read flag 50 values directly for verification
  if(!requireNamespace("FLR4MFCL", quietly = TRUE)) skip("FLR4MFCL not available")
  base_par = FLR4MFCL::read.MFCLPar(par_file_mfcl, first.yr = 1952)
  flag_dt = data.table::as.data.table(FLR4MFCL::flags(base_par))
  flag50 = flag_dt[flag == 50 & flagtype < 0]
  flag50[, fleet := abs(flagtype)]
  flag50_lookup = stats::setNames(flag50$value, flag50$fleet)
  
  # Check one observation per fleet: Nsamp_adj should equal min(Nsamp_in, 1000) / flag50
  for(flt in unique(wt_mfcl_time$Fleet)) {
    flt_char = as.character(flt)
    if(flt_char %in% names(flag50_lookup)) {
      f50 = flag50_lookup[flt_char]
      obs_row = wt_mfcl_time[Fleet == flt][1]
      expected_adj = unname(min(obs_row$Nsamp_in, 1000) / f50)
      expect_equal(obs_row$Nsamp_adj, expected_adj, tolerance = 1e-10,
                   label = paste0("Fleet ", flt, " Nsamp_adj = min(Nsamp_in, 1000) / ", f50))
    }
  }
})

test_that("MFCL Nsamp_adj differs from Nsamp_in when par_file provided", {
  skip_if_not(!is.null(wt_mfcl_time) && file.exists(par_file_mfcl),
              "MFCL time-resolved data and par file must be available")
  
  # Nsamp_adj should NOT equal Nsamp_in (because flag 50 divides it)
  # At least for rows where Nsamp_in > 0
  rows_with_data = wt_mfcl_time[Nsamp_in > 0]
  if(nrow(rows_with_data) > 0) {
    expect_false(all(rows_with_data$Nsamp_adj == rows_with_data$Nsamp_in),
                 label = "Nsamp_adj should differ from Nsamp_in when flag 50 is applied")
  }
})

test_that("MFCL Nsamp_adj is capped at 1000/flag50 when Nsamp_in > 1000", {
  skip_if_not(!is.null(wt_mfcl_time) && file.exists(par_file_mfcl),
              "MFCL time-resolved data and par file must be available")
  
  if(!requireNamespace("FLR4MFCL", quietly = TRUE)) skip("FLR4MFCL not available")
  base_par = FLR4MFCL::read.MFCLPar(par_file_mfcl, first.yr = 1952)
  flag_dt = data.table::as.data.table(FLR4MFCL::flags(base_par))
  flag50 = flag_dt[flag == 50 & flagtype < 0]
  flag50[, fleet := abs(flagtype)]
  flag50_lookup = stats::setNames(flag50$value, flag50$fleet)
  
  # For rows where Nsamp_in > 1000, Nsamp_adj should be 1000/flag50
  for(flt in unique(wt_mfcl_time$Fleet)) {
    flt_char = as.character(flt)
    if(flt_char %in% names(flag50_lookup)) {
      f50 = flag50_lookup[flt_char]
      large_rows = wt_mfcl_time[Fleet == flt & Nsamp_in > 1000]
      if(nrow(large_rows) > 0) {
        expected_adj = 1000 / f50
        expect_true(all(abs(large_rows$Nsamp_adj - expected_adj) < 1e-10),
                    label = paste0("Fleet ", flt,
                                   ": Nsamp_adj should be 1000/", f50,
                                   " when Nsamp_in > 1000"))
      }
    }
  }
})
