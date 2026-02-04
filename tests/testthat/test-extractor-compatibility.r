# Tests for compatibility between SS3 and MFCL biomass extractors

library(testthat)
library(data.table)

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

context("Extractor output compatibility")

test_that("SS3 and MFCL extractors return data.tables with compatible structure", {
	ss3_result = extract_ss3_biomass(file.path(dir_ss3, "01-bet-base"), quarterly = FALSE)
	mfcl_result = extract_mfcl_biomass(file.path(dir_mfcl, "v11", "plot-10.par.rep"), quarterly = FALSE)
	
	# Check that both are data.tables
	expect_s3_class(ss3_result, "data.table")
	expect_s3_class(mfcl_result, "data.table")
})

test_that("SS3 and MFCL extractors return compatible column sets", {
	ss3_result = extract_ss3_biomass(file.path(dir_ss3, "01-bet-base"), quarterly = FALSE)
	mfcl_result = extract_mfcl_biomass(file.path(dir_mfcl, "v11", "plot-10.par.rep"), quarterly = FALSE)
	
	# Check that both have the same columns
	expected_cols = c("model", "year", "ts", "season", "ssb", "ssb_se", "depletion", "depletion_se")
	expect_equal(colnames(ss3_result), expected_cols)
	expect_equal(colnames(mfcl_result), expected_cols)
})

test_that("Common columns have consistent data types between extractors", {
	ss3_result = extract_ss3_biomass(file.path(dir_ss3, "01-bet-base"), quarterly = FALSE)
	mfcl_result = extract_mfcl_biomass(file.path(dir_mfcl, "v11", "plot-10.par.rep"), quarterly = FALSE)
	
	# Check data types for all columns
	for(col in c("model", "year", "ts", "season", "ssb")) {
		expect_identical(class(ss3_result[[col]]), class(mfcl_result[[col]]),
		               info = sprintf("Column '%s' has different types", col))
	}
})

test_that("Extractors produce combinable output (rbindlist compatible)", {
	ss3_result = extract_ss3_biomass(file.path(dir_ss3, "01-bet-base"), quarterly = FALSE)
	mfcl_result = extract_mfcl_biomass(file.path(dir_mfcl, "v11", "plot-10.par.rep"), quarterly = FALSE)
	
	# Both should be directly combinable
	combined = rbindlist(list(ss3_result, mfcl_result))
	
	expect_s3_class(combined, "data.table")
	expect_equal(nrow(combined), nrow(ss3_result) + nrow(mfcl_result))
})

test_that("Extractors produce same value ranges for biomass and depletion", {
	ss3_result = extract_ss3_biomass(file.path(dir_ss3, "01-bet-base"), quarterly = FALSE)
	mfcl_result = extract_mfcl_biomass(file.path(dir_mfcl, "v11", "plot-10.par.rep"), quarterly = FALSE)
	
	# Both should have positive SSB
	expect_true(all(ss3_result$ssb > 0))
	expect_true(all(mfcl_result$ssb > 0))
	
	# Both should have depletion in reasonable range (0 to ~1.1 accounting for recovery), excluding NAs
	expect_true(all(ss3_result$depletion[!is.na(ss3_result$depletion)] >= 0 & ss3_result$depletion[!is.na(ss3_result$depletion)] <= 1.1))
	expect_true(all(mfcl_result$depletion[!is.na(mfcl_result$depletion)] >= 0 & mfcl_result$depletion[!is.na(mfcl_result$depletion)] <= 1.1))
})
