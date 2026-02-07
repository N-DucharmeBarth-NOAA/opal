# Tests for plot_cpue_modular function

library(testthat)
library(data.table)
library(ggplot2)

# Define project paths
proj_dir = this.path::this.proj()
dir_model = file.path(proj_dir, "model-files")
dir_ss3 = file.path(dir_model, "ss3")
dir_mfcl = file.path(dir_model, "mfcl")
dir_helper_fns = file.path(proj_dir, "code", "ss3", "helper-fns")

# Source helper functions
source(file.path(dir_helper_fns, "plot-cpue-modular.r"))

context("plot_cpue_modular")

# Test with file path input
test_that("plot_cpue_modular accepts file path input", {
  cpue_file = file.path(dir_ss3, "01-bet-base", "cpue.csv")
  skip_if_not(file.exists(cpue_file), "CPUE file not found")
  
  p = plot_cpue_modular(cpue_file, verbose = FALSE)
  
  expect_s3_class(p, "ggplot")
})

# Test with data.table input
test_that("plot_cpue_modular accepts data.table input", {
  cpue_file = file.path(dir_ss3, "01-bet-base", "cpue.csv")
  skip_if_not(file.exists(cpue_file), "CPUE file not found")
  
  cpue_dt = fread(cpue_file)
  p = plot_cpue_modular(cpue_dt, verbose = FALSE)
  
  expect_s3_class(p, "ggplot")
})

# Test with data.frame input
test_that("plot_cpue_modular accepts data.frame input", {
  cpue_file = file.path(dir_ss3, "01-bet-base", "cpue.csv")
  skip_if_not(file.exists(cpue_file), "CPUE file not found")
  
  cpue_df = as.data.frame(fread(cpue_file))
  p = plot_cpue_modular(cpue_df, verbose = FALSE)
  
  expect_s3_class(p, "ggplot")
})

# Test column name normalization (lowercase)
test_that("plot_cpue_modular handles lowercase column names", {
  cpue_dt = data.table(
    fleet = c(1, 1, 1),
    time = c(1, 2, 3),
    obs = c(1.5, 1.8, 2.0),
    exp = c(1.4, 1.7, 1.9),
    se = c(0.2, 0.2, 0.2)
  )
  
  p = plot_cpue_modular(cpue_dt, verbose = FALSE)
  expect_s3_class(p, "ggplot")
})

# Test column name normalization (uppercase)
test_that("plot_cpue_modular handles uppercase column names", {
  cpue_dt = data.table(
    Fleet = c(1, 1, 1),
    Time = c(1, 2, 3),
    Obs = c(1.5, 1.8, 2.0),
    Exp = c(1.4, 1.7, 1.9),
    SE = c(0.2, 0.2, 0.2)
  )
  
  p = plot_cpue_modular(cpue_dt, verbose = FALSE)
  expect_s3_class(p, "ggplot")
})

# Test missing Fleet_name column creates default names
test_that("plot_cpue_modular creates default fleet names when missing", {
  cpue_dt = data.table(
    fleet = c(1, 1, 2, 2),
    time = c(1, 2, 1, 2),
    obs = c(1.5, 1.8, 2.0, 2.2),
    exp = c(1.4, 1.7, 1.9, 2.1),
    se = c(0.2, 0.2, 0.2, 0.2)
  )
  
  p = plot_cpue_modular(cpue_dt, verbose = FALSE)
  expect_s3_class(p, "ggplot")
  
  # Check that default fleet names were created
  plot_data = ggplot2::ggplot_build(p)$data[[1]]
  expect_true("PANEL" %in% names(plot_data))
})

# Test fleet_names parameter
test_that("plot_cpue_modular accepts fleet_names parameter", {
  cpue_dt = data.table(
    fleet = c(1, 1, 2, 2),
    time = c(1, 2, 1, 2),
    obs = c(1.5, 1.8, 2.0, 2.2),
    exp = c(1.4, 1.7, 1.9, 2.1),
    se = c(0.2, 0.2, 0.2, 0.2)
  )
  
  p = plot_cpue_modular(cpue_dt, fleet_names = c("Fleet A", "Fleet B"), verbose = FALSE)
  expect_s3_class(p, "ggplot")
})

# Test fleet_names length validation
test_that("plot_cpue_modular validates fleet_names length", {
  cpue_dt = data.table(
    fleet = c(1, 1, 2, 2),
    time = c(1, 2, 1, 2),
    obs = c(1.5, 1.8, 2.0, 2.2),
    exp = c(1.4, 1.7, 1.9, 2.1),
    se = c(0.2, 0.2, 0.2, 0.2)
  )
  
  expect_error(
    plot_cpue_modular(cpue_dt, fleet_names = c("Fleet A"), verbose = FALSE),
    "Length of fleet_names must match number of unique fleet IDs"
  )
})

# Test model_id parameter
test_that("plot_cpue_modular uses model_id when model_label absent", {
  cpue_dt = data.table(
    fleet = c(1, 1, 1),
    time = c(1, 2, 3),
    obs = c(1.5, 1.8, 2.0),
    exp = c(1.4, 1.7, 1.9),
    se = c(0.2, 0.2, 0.2)
  )
  
  p = plot_cpue_modular(cpue_dt, model_id = "test_model", verbose = FALSE)
  expect_s3_class(p, "ggplot")
})

# Test log scale transformation
test_that("plot_cpue_modular applies log scale transformation", {
  cpue_dt = data.table(
    fleet = c(1, 1, 1),
    time = c(1, 2, 3),
    obs = c(1.5, 1.8, 2.0),
    exp = c(1.4, 1.7, 1.9),
    se = c(0.2, 0.2, 0.2)
  )
  
  p = plot_cpue_modular(cpue_dt, use_log_scale = TRUE, verbose = FALSE)
  expect_s3_class(p, "ggplot")
  
  # Check y-axis label
  expect_match(p$labels$y, "log-scale")
})

# Test show_se parameter
test_that("plot_cpue_modular respects show_se parameter", {
  cpue_dt = data.table(
    fleet = c(1, 1, 1),
    time = c(1, 2, 3),
    obs = c(1.5, 1.8, 2.0),
    exp = c(1.4, 1.7, 1.9),
    se = c(0.2, 0.2, 0.2)
  )
  
  p_with_se = plot_cpue_modular(cpue_dt, show_se = TRUE, verbose = FALSE)
  p_without_se = plot_cpue_modular(cpue_dt, show_se = FALSE, verbose = FALSE)
  
  expect_s3_class(p_with_se, "ggplot")
  expect_s3_class(p_without_se, "ggplot")
  
  # With SE should have more layers (error bars)
  expect_true(length(p_with_se$layers) > length(p_without_se$layers))
})

# Test multiple models
test_that("plot_cpue_modular handles multiple models", {
  cpue_dt1 = data.table(
    fleet = c(1, 1, 1),
    time = c(1, 2, 3),
    obs = c(1.5, 1.8, 2.0),
    exp = c(1.4, 1.7, 1.9),
    se = c(0.2, 0.2, 0.2),
    model_label = "Model A"
  )
  
  cpue_dt2 = data.table(
    fleet = c(1, 1, 1),
    time = c(1, 2, 3),
    obs = c(1.5, 1.8, 2.0),
    exp = c(1.5, 1.8, 2.0),
    se = c(0.2, 0.2, 0.2),
    model_label = "Model B"
  )
  
  cpue_combined = rbindlist(list(cpue_dt1, cpue_dt2), use.names = TRUE)
  p = plot_cpue_modular(cpue_combined, verbose = FALSE)
  
  expect_s3_class(p, "ggplot")
})

# Test observed data aggregation (no duplicates)
test_that("plot_cpue_modular aggregates observed data correctly", {
  # Create data with same observations for two models
  cpue_dt1 = data.table(
    fleet = c(1, 1),
    time = c(1, 2),
    obs = c(1.5, 1.8),
    exp = c(1.4, 1.7),
    se = c(0.2, 0.2),
    model_label = "Model A"
  )
  
  cpue_dt2 = data.table(
    fleet = c(1, 1),
    time = c(1, 2),
    obs = c(1.5, 1.8),
    exp = c(1.5, 1.9),
    se = c(0.2, 0.2),
    model_label = "Model B"
  )
  
  cpue_combined = rbindlist(list(cpue_dt1, cpue_dt2), use.names = TRUE)
  p = plot_cpue_modular(cpue_combined, verbose = FALSE)
  
  expect_s3_class(p, "ggplot")
  
  # Check that observed points layer has correct number of points
  # (should be unique fleet/time combinations, not duplicated for each model)
  built_plot = ggplot2::ggplot_build(p)
  
  # Find the points layer (should be layer with shape 21)
  points_layer = NULL
  for (i in seq_along(built_plot$data)) {
    if ("shape" %in% names(built_plot$data[[i]]) && 
        all(built_plot$data[[i]]$shape == 21)) {
      points_layer = built_plot$data[[i]]
      break
    }
  }
  
  expect_false(is.null(points_layer))
  expect_equal(nrow(points_layer), 2)  # Should have 2 unique time points, not 4
})

# Test required columns validation
test_that("plot_cpue_modular validates required columns", {
  cpue_dt = data.table(
    fleet = c(1, 1, 1),
    time = c(1, 2, 3),
    obs = c(1.5, 1.8, 2.0)
    # Missing exp and se columns
  )
  
  expect_error(
    plot_cpue_modular(cpue_dt, verbose = FALSE),
    "missing required columns"
  )
})

# Test invalid input type
test_that("plot_cpue_modular rejects invalid input type", {
  expect_error(
    plot_cpue_modular(list(a = 1, b = 2), verbose = FALSE),
    "cpue_input must be a file path or a data.table/data.frame"
  )
})

# Test nonexistent file path
test_that("plot_cpue_modular rejects nonexistent file path", {
  expect_error(
    plot_cpue_modular("/nonexistent/path/to/cpue.csv", verbose = FALSE),
    "CPUE file not found"
  )
})

# Test n_col parameter
test_that("plot_cpue_modular respects n_col parameter", {
  cpue_dt = data.table(
    fleet = c(1, 1, 2, 2),
    time = c(1, 2, 1, 2),
    obs = c(1.5, 1.8, 2.0, 2.2),
    exp = c(1.4, 1.7, 1.9, 2.1),
    se = c(0.2, 0.2, 0.2, 0.2)
  )
  
  p = plot_cpue_modular(cpue_dt, n_col = 1, verbose = FALSE)
  expect_s3_class(p, "ggplot")
})

# Integration test with real data from both SS3 and MFCL
test_that("plot_cpue_modular works with real SS3 and MFCL data", {
  ss3_file = file.path(dir_ss3, "01-bet-base", "cpue.csv")
  mfcl_file = file.path(dir_mfcl, "v11", "cpue.csv")
  
  skip_if_not(file.exists(ss3_file), "SS3 CPUE file not found")
  skip_if_not(file.exists(mfcl_file), "MFCL CPUE file not found")
  
  # Read both files
  ss3_cpue = fread(ss3_file)
  ss3_cpue[, model_label := "SS3"]
  
  mfcl_cpue = fread(mfcl_file)
  mfcl_cpue[, model_label := "MFCL"]
  
  # Combine
  cpue_combined = rbindlist(list(ss3_cpue, mfcl_cpue), use.names = TRUE, fill = TRUE)
  
  # Plot
  p = plot_cpue_modular(cpue_combined, verbose = FALSE)
  
  expect_s3_class(p, "ggplot")
})
