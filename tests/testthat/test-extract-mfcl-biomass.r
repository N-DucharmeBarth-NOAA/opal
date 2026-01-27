# Tests for extract_mfcl_biomass function

library(testthat)
library(data.table)
library(magrittr)
library(FLR4MFCL)

# Define project paths
proj_dir = this.path::this.proj()
dir_model = file.path(proj_dir, "model-files")
dir_mfcl = file.path(dir_model, "mfcl")
dir_helper_fns = file.path(proj_dir, "code", "ss3", "helper-fns")

# Source helper functions
sapply(file.path(dir_helper_fns, list.files(dir_helper_fns)), source)

context("extract_mfcl_biomass")

test_that("extract_mfcl_biomass returns data.table with correct structure", {
	mfcl_rep = file.path(dir_mfcl, "v11", "plot-10.par.rep")
	result = extract_mfcl_biomass(mfcl_rep, model_name = "MFCL-test")
	
	expect_s3_class(result, "data.table")
	expect_named(result, c("model", "year", "ssb", "ssb_se", "depletion", "depletion_se"))
})

test_that("extract_mfcl_biomass returns correct data types", {
	mfcl_rep = file.path(dir_mfcl, "v11", "plot-10.par.rep")
	result = extract_mfcl_biomass(mfcl_rep, model_name = "MFCL-test")
	
	expect_is(result$model, "character")
	expect_is(result$year, "numeric")
	expect_is(result$ssb, "numeric")
	expect_is(result$ssb_se, "numeric")
	expect_is(result$depletion, "numeric")
	expect_is(result$depletion_se, "numeric")
})

test_that("extract_mfcl_biomass returns positive biomass values", {
	mfcl_rep = file.path(dir_mfcl, "v11", "plot-10.par.rep")
	result = extract_mfcl_biomass(mfcl_rep)
	
	expect_true(all(result$ssb > 0))
})

test_that("extract_mfcl_biomass returns depletion between 0 and 1", {
	mfcl_rep = file.path(dir_mfcl, "v11", "plot-10.par.rep")
	result = extract_mfcl_biomass(mfcl_rep)
	
	expect_true(all(result$depletion >= 0 & result$depletion <= 1))
})

test_that("extract_mfcl_biomass returns no NA values", {
	mfcl_rep = file.path(dir_mfcl, "v11", "plot-10.par.rep")
	result = extract_mfcl_biomass(mfcl_rep)
	
	expect_false(anyNA(result))
})

test_that("extract_mfcl_biomass model name is preserved", {
	mfcl_rep = file.path(dir_mfcl, "v11", "plot-10.par.rep")
	result = extract_mfcl_biomass(mfcl_rep, model_name = "MFCL-custom")
	
	expect_equal(unique(result$model), "MFCL-custom")
})

test_that("extract_mfcl_biomass years are sorted", {
	mfcl_rep = file.path(dir_mfcl, "v11", "plot-10.par.rep")
	result = extract_mfcl_biomass(mfcl_rep)
	
	expect_equal(result$year, sort(result$year))
})

test_that("extract_mfcl_biomass stops if report file not found", {
	expect_error(extract_mfcl_biomass("nonexistent/path/to/file.rep"))
})

test_that("extract_mfcl_biomass depletion decreases monotonically (fishing removes biomass)", {
	mfcl_rep = file.path(dir_mfcl, "v11", "plot-10.par.rep")
	result = extract_mfcl_biomass(mfcl_rep)
	
	# Check if mean depletion in recent years is lower than early years
	n = nrow(result)
	early_mean = mean(result$depletion[1:min(10, n)])
	recent_mean = mean(result$depletion[max(n-9, 1):n])
	
	expect_true(recent_mean <= early_mean, 
	           info = sprintf("Early mean depletion: %.3f, Recent mean: %.3f", early_mean, recent_mean))
})
