composition_reporting_fixture <- function(type = "lf", family = 1L, legacy = FALSE) {
  data <- list(years = 2001:2002, first_yr = 2001L, removal_switch_f = c(0L, 0L),
    len_bin_start = 10, len_bin_width = 2, n_len = 4L,
    wt_bin_start = 10, wt_bin_width = 2, n_wt = 4L)
  predictions <- if (legacy) rep(list(c(.4, .6)), 4L) else
    list(c(.4, .6), c(.2, .3, .5), c(.2, .3, .5), c(.4, .6))
  fields <- list(switch = family, fishery_f = 1:2, n_f = c(2L, 2L),
    year_fi = list(1:2, 1:2), n_fi = list(c(100, 100), c(100, 100)),
    n_int_fi = list(c(100L, 100L), c(100L, 100L)),
    row_fi = list(c(2L, 4L), c(1L, 3L)), minbin = c(2L, 2L), maxbin = c(3L, 3L),
    obs_flat = unlist(predictions) * 100, obs_ints = unlist(predictions) * 100,
    obs_prop = unlist(predictions))
  if (!legacy) {
    fields$scheme_edges <- list(c(10, 12, 18), c(10, 11, 14, 18))
    fields$scheme <- c(1L, 2L, 2L, 1L)
    fields$nbins_i <- c(2L, 3L, 3L, 2L)
  }
  names(fields) <- paste0(type, "_", names(fields))
  data[names(fields)] <- fields
  data[[paste0("n_", type)]] <- 4L
  parameters <- setNames(list(log(c(20, 20))), paste0("log_", type, "_tau"))
  report <- if (legacy) setNames(list(list(matrix(c(.4, .6), 2, 2, byrow = TRUE),
    matrix(c(.4, .6), 2, 2, byrow = TRUE))), paste0(type, "_pred")) else
    setNames(list(predictions), paste0(type, "_pred_row"))
  object <- list(env = list(last.par.best = numeric(), parList = function(par) parameters),
    report = function(par) report)
  list(data = data, parameters = parameters, report = report, object = object,
    x = list(data = data, random = character()))
}

test_that("mixed-scheme rows retain input ordering and grouped observation offsets", {
  for (type in c("lf", "wf")) {
    fixture <- composition_reporting_fixture(type)
    rows <- opal:::.opal_composition_rows(fixture$data, fixture$parameters, fixture$report, type)
    expect_identical(vapply(rows, `[[`, integer(1), "row"), 1:4)
    expect_equal(lengths(lapply(rows, `[[`, "bin")), c(3L, 2L, 2L, 3L))
    expect_equal(rows[[1]]$observation_index, 6:8)
    expect_equal(rows[[2]]$observation_index, 1:2)
    expect_equal(rows[[1]]$bin_lower, c(10, 11, 14))
    expect_equal(rows[[1]]$bin_upper, c(11, 14, 18))
    expect_equal(rows[[1]]$predicted, c(.2, .3, .5))
  }
})

test_that("ragged OSA preserves row order and excludes each constrained last bin", {
  for (type in c("lf", "wf")) for (family in 1:3) {
    fixture <- composition_reporting_fixture(type, family)
    local_mocked_bindings(opal_rtmb = function(...) fixture$object, .package = "opal")
    result <- opal:::.opal_with_seed(42, opal:::.opal_residual_result(fixture$x, .95))
    df <- result$residuals
    expect_equal(nrow(df), 10L)
    expect_equal(df$row, rep(1:4, c(3, 2, 2, 3)))
    expect_equal(df$bin, c(1:3, 1:2, 1:2, 1:3))
    expect_equal(which(df$reason == "constrained final bin"), c(3, 5, 7, 10))
    expect_equal(sum(is.finite(df$residual)), 6L)
    expect_equal(df$bin_upper - df$bin_lower, c(1, 3, 4, 2, 6, 2, 6, 1, 3, 4))
  }
})

test_that("composition scales build at physical midpoints with visible extents", {
  for (type in c("lf", "wf")) {
    fixture <- composition_reporting_fixture(type)
    local_mocked_bindings(opal_rtmb = function(...) fixture$object, .package = "opal")
    for (scale in c("auto", "proportion", "density")) {
      plot <- plot_composition(fixture$x, type, scale = scale)
      built <- ggplot2::ggplot_build(plot)
      expect_length(built$data, 4L)
      expect_equal(plot$data$bin_mid, (plot$data$bin_lower + plot$data$bin_upper) / 2)
      expect_equal(sort(built$data[[3]]$xend - built$data[[3]]$x),
        sort(plot$data$bin_upper - plot$data$bin_lower))
      expect_equal(sort(built$data[[4]]$xend - built$data[[4]]$x),
        sort(plot$data$bin_upper - plot$data$bin_lower))
      expect_equal(built$data[[4]]$y, built$data[[1]]$y)
      density <- scale != "proportion"
      expect_equal(plot$labels$y, if (density) "Density" else "Proportion")
      totals <- tapply(plot$data$observed * if (density)
        plot$data$bin_upper - plot$data$bin_lower else 1, plot$data$row, sum)
      expect_equal(as.numeric(totals), rep(1, 4))
      if (density) expect_match(plot$labels$caption, "bin width")
    }
    expect_equal(unique(plot_composition(fixture$x, type, fishery = 2)$data$fishery), 2L)
  }
})

test_that("legacy reports retain bin indices, probabilities and OSA results", {
  for (type in c("lf", "wf")) {
    fixture <- composition_reporting_fixture(type, legacy = TRUE)
    local_mocked_bindings(opal_rtmb = function(...) fixture$object, .package = "opal")
    rows <- opal:::.opal_composition_rows(fixture$data, fixture$parameters, fixture$report, type)
    expect_equal(rows[[1]]$bin, 2:3)
    expect_equal(rows[[1]]$bin_lower, c(12, 14))
    expect_equal(rows[[1]]$bin_upper, c(14, 16))
    plot <- plot_composition(fixture$x, type)
    expect_equal(plot$labels$y, "Proportion")
    expect_equal(plot$data$bin_mid, rep(c(13, 15), 4))
    expect_equal(plot$data$observed, rep(c(.4, .6), 4))
    expect_no_error(ggplot2::ggplot_build(plot))
    expected <- opal:::.opal_with_seed(42, replicate(4,
      opal:::.opal_composition_pit(c(40, 60), c(.4, .6), 1)[1]))
    result <- opal:::.opal_with_seed(42, opal:::.opal_residual_result(fixture$x, .95))
    expect_equal(result$residuals$residual[result$residuals$bin == 2], expected)
    expect_equal(sum(is.na(result$residuals$residual)), 4L)
  }
})

test_that("auto density detects different equal-width schemes in a fishery-year", {
  fixture <- composition_reporting_fixture()
  fixture$data$lf_scheme_edges <- list(c(10, 12, 14), c(11, 13, 15, 17))
  fixture$data$lf_year_fi <- list(c(1L, 1L), c(1L, 1L))
  fixture$x$data <- fixture$data
  local_mocked_bindings(opal_rtmb = function(...) fixture$object, .package = "opal")
  expect_equal(plot_composition(fixture$x)$labels$y, "Density")
  expect_no_error(ggplot2::ggplot_build(plot_composition(fixture$x)))
  expect_equal(plot_composition(fixture$x, scale = "proportion")$labels$y, "Proportion")
})

test_that("uniform rows with different widths in separate panels use proportions", {
  fixture <- composition_reporting_fixture()
  fixture$data$lf_scheme_edges <- list(c(10, 12, 14), c(10, 13, 16, 19))
  fixture$x$data <- fixture$data
  local_mocked_bindings(opal_rtmb = function(...) fixture$object, .package = "opal")
  plot <- plot_composition(fixture$x)
  expect_equal(plot$labels$y, "Proportion")
  expect_no_error(ggplot2::ggplot_build(plot))
})

test_that("native prepared schemes preserve historical aggregated bin indices", {
  fixture <- composition_reporting_fixture(legacy = TRUE)
  fixture$data$lf_scheme_edges <- list(c(12, 14, 16))
  fixture$data$lf_scheme <- rep(1L, 4)
  fixture$data$lf_nbins_i <- rep(2L, 4)
  fixture$report$lf_pred_row <- rep(list(c(.4, .6)), 4)
  rows <- opal:::.opal_composition_rows(fixture$data, fixture$parameters, fixture$report, "lf")
  expect_equal(rows[[1]]$bin, 2:3)
  expect_equal(rows[[1]]$bin_lower, c(12, 14))
  expect_equal(rows[[1]]$bin_upper, c(14, 16))
})

test_that("the coarse-grid example fits and plots all scales", {
  inputs <- opal_example_inputs()
  compositions <- data.frame(fishery = 1L, ts = c(1L, 1L, 2L, 2L, 2L),
    bin_lower = c(20, 60, 20, 40, 60), bin_upper = c(60, 80, 40, 60, 80),
    count = c(85, 35, 20, 60, 40))
  compositions$year <- inputs$data$first_yr + compositions$ts - 1L
  compositions$month <- 1L
  data <- prep_lf_data(inputs$data, compositions)
  fit <- opal_fit(opal_obj(data, inputs$parameters, inputs$map), check = FALSE)
  for (scale in c("auto", "proportion", "density")) {
    plot <- plot_composition(fit, scale = scale)
    expect_equal(plot$data$bin_lower, compositions$bin_lower)
    expect_equal(plot$data$bin_upper, compositions$bin_upper)
    expect_no_error(ggplot2::ggplot_build(plot))
    expect_equal(plot$labels$y, if (scale == "proportion") "Proportion" else "Density")
  }
  residuals <- opal_derived(opal_osa(fit), "osa")$residuals
  lf <- residuals[residuals$type == "lf", ]
  expect_equal(nrow(lf), 5L)
  expect_equal(lf$bin_lower, compositions$bin_lower)
  expect_equal(sum(is.finite(lf$residual)), 3L)
})
