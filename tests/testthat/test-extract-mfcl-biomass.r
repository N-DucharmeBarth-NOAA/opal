# Tests for extract_mfcl_biomass function

library(testthat)
library(data.table)
library(magrittr)
library(FLR4MFCL)

# Define project paths
proj_dir = this.path::this.proj()
dir_model = file.path(proj_dir, "model-files")
dir_mfcl = file.path(dir_model, "mfcl")
dir_mfcl_helper_fns = file.path(proj_dir, "code", "mfcl", "helper-fns")

# Source helper functions
sapply(file.path(dir_mfcl_helper_fns, list.files(dir_mfcl_helper_fns)), source)

context("extract_mfcl_biomass")

# Extract data once to reuse across tests (improves test performance)
mfcl_rep = file.path(dir_mfcl, "v11", "plot-10.par.rep")
result_quarterly = extract_mfcl_biomass(mfcl_rep, model_name = "MFCL-test", quarterly = TRUE)

test_that("extract_mfcl_biomass returns data.table with correct structure", {
	expect_s3_class(result_quarterly, "data.table")
	expect_named(result_quarterly, c("model", "year","ts", "season", "ssb", "ssb_se", "depletion", "depletion_se"))
})

test_that("extract_mfcl_biomass returns correct data types", {
	expect_is(result_quarterly$model, "character")
	expect_is(result_quarterly$year, "numeric")
	expect_is(result_quarterly$ts, "numeric")
	expect_is(result_quarterly$season, "numeric")	
	expect_is(result_quarterly$ssb, "numeric")
	expect_is(result_quarterly$ssb_se, "numeric")
	expect_is(result_quarterly$depletion, "numeric")
	expect_is(result_quarterly$depletion_se, "numeric")

})

test_that("extract_mfcl_biomass returns positive biomass values", {
	expect_true(all(result_quarterly$ssb > 0))
})

test_that("extract_mfcl_biomass returns depletion between 0 and 1", {
	expect_true(all(result_quarterly$depletion >= 0 & result_quarterly$depletion <= 1.1), 
	           info = sprintf("Found depletion values: min=%.3f, max=%.3f", 
	                         min(result_quarterly$depletion), max(result_quarterly$depletion)))
})

test_that("extract_mfcl_biomass returns no NA values", {
	expect_false(anyNA(result_quarterly[,.(model, year, ts, season, ssb, depletion)]))
})

test_that("extract_mfcl_biomass model name is preserved", {
	expect_equal(unique(result_quarterly$model), "MFCL-test")
})

test_that("extract_mfcl_biomass years are sorted", {
	expect_equal(result_quarterly$year, sort(result_quarterly$year))
})

test_that("extract_mfcl_biomass stops if report file not found", {
	expect_error(extract_mfcl_biomass("nonexistent/path/to/file.rep", quarterly = TRUE))
})

test_that("extract_mfcl_biomass depletion decreases monotonically (fishing removes biomass)", {
	# Check if mean depletion in recent years is lower than early years
	n = nrow(result_quarterly)
	early_mean = mean(result_quarterly$depletion[1:min(10, n)])
	recent_mean = mean(result_quarterly$depletion[max(n-9, 1):n])
	
	expect_true(recent_mean <= early_mean, 
	           info = sprintf("Early mean depletion: %.3f, Recent mean: %.3f", early_mean, recent_mean))
})
