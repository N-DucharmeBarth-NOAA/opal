# Tests for extract_ss3_biomass function

library(testthat)
library(data.table)
library(magrittr)
library(r4ss)

# Define project paths
proj_dir = this.path::this.proj()
dir_model = file.path(proj_dir, "model-files")
dir_ss3 = file.path(dir_model, "ss3")
dir_helper_fns = file.path(proj_dir, "code", "ss3", "helper-fns")

# Source helper functions
sapply(file.path(dir_helper_fns, list.files(dir_helper_fns)), source)

context("extract_ss3_biomass")

# Compute result once and reuse across tests to avoid repeatedly parsing Report.sso
result = extract_ss3_biomass(file.path(dir_ss3, "01-bet-base"))

test_that("extract_ss3_biomass returns data.table with correct structure", {
	
	expect_s3_class(result, "data.table")
	expect_named(result, c("model", "year", "ts", "season", "ssb", "ssb_se", "depletion", "depletion_se"))
})

test_that("extract_ss3_biomass returns correct data types", {
	
	expect_is(result$model, "character")
	expect_is(result$year, "numeric")
	expect_is(result$ssb, "numeric")
	expect_is(result$ssb_se, "numeric")
	expect_is(result$depletion, "numeric")
	expect_is(result$depletion_se, "numeric")
})

test_that("extract_ss3_biomass returns positive biomass values", {
	
	expect_true(all(result$ssb > 0), 
	            info = paste("Found SSB values <= 0:", paste(result[ssb <= 0]$ssb, collapse = ", ")))
})

test_that("extract_ss3_biomass returns depletion between 0 and 1", {
	
	expect_true(all(result$depletion[!is.na(result$depletion)] >= 0 & result$depletion[!is.na(result$depletion)] <= 1), 
	            info = paste("Found depletion outside [0,1]:", 
	                         paste(result[!is.na(depletion) & (depletion < 0 | depletion > 1)]$depletion, collapse = ", ")))
})

test_that("extract_ss3_biomass returns no NA values", {
	result_subset = result[,.(model, year, ts, season, ssb, ssb_se)]
	
	expect_false(anyNA(result_subset), 
	            info = paste("Found NA values in columns:", 
	                         paste(colnames(result_subset)[colSums(is.na(result_subset)) > 0], collapse = ", ")))
})

test_that("extract_ss3_biomass model name matches directory", {
	
	expect_equal(unique(result$model), "01-bet-base")
})

test_that("extract_ss3_biomass years are sorted", {
	
	expect_equal(result$ts, sort(result$ts))
})

test_that("extract_ss3_biomass standard errors are positive", {
	result_subset = result[,.(model, year, ts, season, ssb_se, depletion_se)]
	
	expect_true(all(result_subset$ssb_se > 0))
	expect_true(all(result_subset$depletion_se[-1] > 0))
})

test_that("extract_ss3_biomass stops if Report.sso not found", {
	expect_error(extract_ss3_biomass(file.path(dir_ss3, "nonexistent-model")))
})
