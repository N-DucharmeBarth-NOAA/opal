test_that("selectivity age initialisation pads outside old source range", {
  source_ages <- 2:9

  expect_equal(.match_selectivity_age_indices(2:9, source_ages), 2:9)
  expect_equal(.match_selectivity_age_indices(2:10, source_ages), c(2:9, 9))
  expect_equal(.match_selectivity_age_indices(0:3, source_ages), c(2, 2, 2, 3))
})

test_that("selectivity year initialisation pads outside old source range", {
  expect_equal(.match_selectivity_year_indices(1931:1933, 1931, 92), 1:3)
  expect_equal(.match_selectivity_year_indices(c(1930, 1931, 2022, 2023), 1931, 92),
               c(1L, 1L, 92L, 92L))
  expect_equal(.match_selectivity_year_indices(integer(), 1931, 92), integer())
})

test_that("get_map only maps supported current parameters", {
  parameters <- list(
    log_B0 = 20,
    log_h = log(0.95),
    log_sigma_r = log(0.6),
    log_cpue_q = c(0, 0),
    log_lf_tau = rep(0, 3),
    par_sel = matrix(0, 3, 6),
    rdev_y = rep(0, 10)
  )

  map <- get_map(parameters)

  expect_named(map, c("log_h", "log_sigma_r", "log_lf_tau", "par_sel"))
  expect_true(all(vapply(map, function(x) all(is.na(x)), logical(1))))
  expect_length(map$log_lf_tau, 3)
  expect_length(map$par_sel, 18)
  expect_false(any(grepl("^par_log_", names(map))))
})

test_that("get_bounds applies bounds to current parameter names", {
  par <- c(
    log_B0 = 20,
    log_h = log(0.95),
    log_sigma_r = log(0.6),
    log_cpue_q = 0,
    rdev_y = 0,
    par_sel = 0,
    log_L1 = log(30)
  )
  bounds <- get_bounds(list(par = par), parameters = as.list(par))

  expect_equal(bounds$lower[bounds$parameter == "log_B0"], 0)
  expect_equal(bounds$upper[bounds$parameter == "log_B0"], 22)
  expect_equal(bounds$lower[bounds$parameter == "rdev_y"], -5)
  expect_equal(bounds$upper[bounds$parameter == "rdev_y"], 5)
  expect_equal(bounds$lower[bounds$parameter == "par_sel"], -7)
  expect_equal(bounds$upper[bounds$parameter == "par_sel"], 7)
  expect_equal(bounds$lower[bounds$parameter == "log_L1"], -Inf)
  expect_equal(bounds$upper[bounds$parameter == "log_L1"], Inf)
})
