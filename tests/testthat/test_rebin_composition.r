# Tests for rebin_composition function

library(testthat)

# Source the function from project root so tests run from `tests/testthat` correctly
source(file.path(this.path::this.proj(), "code", "ss3", "helper-fns", "rebin_composition.r"))

context("rebin_composition")

test_that("Total count is preserved", {
  src_edges = c(10, 20, 30, 40)
  src_counts = c(5, 10, 15)
  dest_edges = c(10, 25, 35, 40)
  
  result = rebin_composition(src_edges, src_counts, dest_edges)
  
  expect_equal(sum(src_counts), sum(result))
})

test_that("Identical bins return identical counts", {
  src_edges = c(0, 10, 20, 30)
  src_counts = c(3, 7, 12)
  
  result = rebin_composition(src_edges, src_counts, src_edges)
  
  expect_equal(src_counts, result)
})

test_that("Known rebinning example produces correct result", {
  # Two 10cm bins (10-20, 20-30) with 10 fish each
  src_edges = c(10, 20, 30)
  src_counts = c(10, 10)
  
  # One 20cm bin (10-30)
  dest_edges = c(10, 30)
  
  result = rebin_composition(src_edges, src_counts, dest_edges)
  
  expect_equal(result, c(20))
})

test_that("Proportions sum to 1 after rebinning", {
  src_edges = c(50, 60, 70, 80, 90)
  src_props = c(0.2, 0.3, 0.4, 0.1)
  dest_edges = c(50, 65, 85, 90)
  
  result = rebin_composition(src_edges, src_props, dest_edges)
  
  expect_equal(sum(result), 1.0, tolerance = 1e-10)
})

test_that("Non-overlapping bins produce zeros", {
  src_edges = c(10, 20, 30)
  src_counts = c(5, 10)
  dest_edges = c(40, 50, 60)  # No overlap
  
  result = rebin_composition(src_edges, src_counts, dest_edges)
  
  expect_equal(result, c(0, 0))
})

test_that("Error on mismatched input lengths", {
  expect_error(
    rebin_composition(c(1, 2, 3), c(10, 20, 30), c(1, 2)),
    "src_edges must have length"
  )
})

test_that("Error on non-increasing bin edges", {
  expect_error(
    rebin_composition(c(10, 20, 15), c(5, 10), c(10, 20)),
    "src_edges must be strictly increasing"
  )
})

test_that("Error on negative counts", {
  expect_error(
    rebin_composition(c(1, 2, 3), c(5, -3), c(1, 2, 3)),
    "src_counts must be non-negative"
  )
})

test_that("Half-overlap distributes correctly", {
  # Source: 0-10 with 10 fish
  src_edges = c(0, 10)
  src_counts = c(10)
  
  # Dest: 0-5, 5-10 (should split evenly)
  dest_edges = c(0, 5, 10)
  
  result = rebin_composition(src_edges, src_counts, dest_edges)
  
  expect_equal(result, c(5, 5))
})

test_that("Handles very small numbers without precision loss", {
  src_edges = c(100, 110, 120)
  src_counts = c(1e-10, 2e-10)
  dest_edges = c(100, 115, 120)
  
  result = rebin_composition(src_edges, src_counts, dest_edges)
  
  expect_equal(sum(src_counts), sum(result), tolerance = 1e-15)
})

test_that("Single source bin to multiple destination bins", {
  # One source bin split into multiple destination bins
  src_edges = c(0, 30)
  src_counts = c(30)
  dest_edges = c(0, 10, 20, 30)
  
  result = rebin_composition(src_edges, src_counts, dest_edges)
  
  expect_equal(result, c(10, 10, 10))
  expect_equal(sum(src_counts), sum(result))
})

test_that("Multiple source bins to single destination bin", {
  # Multiple source bins combined into one destination bin
  src_edges = c(0, 10, 20, 30)
  src_counts = c(5, 10, 15)
  dest_edges = c(0, 30)
  
  result = rebin_composition(src_edges, src_counts, dest_edges)
  
  expect_equal(result, c(30))
  expect_equal(sum(src_counts), sum(result))
})

test_that("Zero counts handled correctly", {
  src_edges = c(10, 20, 30, 40)
  src_counts = c(0, 10, 0)
  dest_edges = c(10, 25, 40)
  
  result = rebin_composition(src_edges, src_counts, dest_edges)
  
  expect_equal(sum(src_counts), sum(result))
  expect_true(all(result >= 0))
})

test_that("Partial overlap distributes proportionally", {
  # Source bin 10-20 with count 20
  src_edges = c(10, 20)
  src_counts = c(20)
  
  # Dest bins: 5-15 and 15-25 (each overlaps 5cm of the 10cm source bin)
  dest_edges = c(5, 15, 25)
  
  result = rebin_composition(src_edges, src_counts, dest_edges)
  
  # First dest bin overlaps [10,15] = 5cm out of 10cm = 50% = 10 fish
  # Second dest bin overlaps [15,20] = 5cm out of 10cm = 50% = 10 fish
  expect_equal(result, c(10, 10))
  expect_equal(sum(src_counts), sum(result))
})

test_that("Complex overlapping scenario", {
  # Source: 3 bins with different widths
  src_edges = c(10, 15, 25, 30)
  src_counts = c(10, 20, 5)  # Total = 35
  
  # Dest: 2 bins with different boundaries
  dest_edges = c(10, 20, 30)
  
  result = rebin_composition(src_edges, src_counts, dest_edges)
  
  # First dest bin [10,20]:
  #   - Source bin 1 [10,15]: full overlap 5cm, contributes 10
  #   - Source bin 2 [15,25]: overlap [15,20] = 5cm out of 10cm = 50%, contributes 10
  #   Total: 20
  # Second dest bin [20,30]:
  #   - Source bin 2 [15,25]: overlap [20,25] = 5cm out of 10cm = 50%, contributes 10
  #   - Source bin 3 [25,30]: full overlap 5cm, contributes 5
  #   Total: 15
  
  expect_equal(result, c(20, 15))
  expect_equal(sum(src_counts), sum(result))
})

test_that("Very large numbers maintained", {
  src_edges = c(0, 100, 200)
  src_counts = c(1e9, 2e9)
  dest_edges = c(0, 150, 200)
  
  result = rebin_composition(src_edges, src_counts, dest_edges)
  
  expect_equal(sum(src_counts), sum(result), tolerance = 1e-6)
})

test_that("Error on insufficient dest_edges", {
  expect_error(
    rebin_composition(c(1, 2, 3), c(5, 10), c(1)),
    "dest_edges must have at least 2 elements"
  )
})

test_that("Error on non-increasing dest_edges", {
  expect_error(
    rebin_composition(c(1, 2, 3), c(5, 10), c(1, 3, 2)),
    "dest_edges must be strictly increasing"
  )
})

test_that("Dest bins extend beyond source range", {
  # Source bins don't cover entire dest range
  src_edges = c(20, 30, 40)
  src_counts = c(10, 20)
  dest_edges = c(10, 25, 35, 50)
  
  result = rebin_composition(src_edges, src_counts, dest_edges)
  
  # Dest bin 1 [10,25]: overlaps [20,25] from source bin 1 = 5cm out of 10cm = 50% = 5
  # Dest bin 2 [25,35]: overlaps [25,30] from source bin 1 = 5cm out of 10cm = 50% = 5
  #                     overlaps [30,35] from source bin 2 = 5cm out of 10cm = 50% = 10
  #                     Total: 15
  # Dest bin 3 [35,50]: overlaps [35,40] from source bin 2 = 5cm out of 10cm = 50% = 10
  
  expect_equal(result, c(5, 15, 10))
  expect_equal(sum(src_counts), sum(result))
})
