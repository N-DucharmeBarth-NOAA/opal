test_that("get_data loads bundled model data by name", {
  expect_equal(get_data("wcpo_bet")$n_fishery, 15L)
  expect_equal(get_data("bet")$n_year, 268L)
  expect_equal(get_data("opakapaka")$n_fishery, 3L)
  expect_equal(get_data("baseline")$n_len, 95L)
})

test_that("get_data can return matching parameters", {
  x <- get_data("opaka", include_parameters = TRUE)

  expect_named(x, c("data", "parameters"))
  expect_equal(x$data$n_fishery, 3L)
  expect_true(all(c("log_B0", "par_sel", "rdev_y") %in% names(x$parameters)))
})

test_that("get_parameters loads by model and bundled data", {
  data(wcpo_bet_data, package = "opal", envir = environment())
  data(wcpo_bet_parameters, package = "opal", envir = environment())
  data(opal_baseline_data, package = "opal", envir = environment())
  data(opal_baseline_parameters, package = "opal", envir = environment())

  expect_equal(get_parameters(model = "wcpo_bet"), wcpo_bet_parameters)
  expect_equal(get_parameters("bet"), wcpo_bet_parameters)
  expect_equal(get_parameters(data = wcpo_bet_data), wcpo_bet_parameters)
  expect_equal(get_parameters(data = opal_baseline_data), opal_baseline_parameters)
})

test_that("get_data rejects legacy list input clearly", {
  expect_error(get_data(list(last_yr = 2022)), "loads bundled model data by name")
})

test_that("bundled aliases, defaults, and invalid requests are resolved clearly", {
  expect_identical(get_data(), get_data("opal_baseline"))
  expect_identical(get_data("WCPO-BET"), get_data("bet"))
  expect_identical(get_data("opal"), get_data("baseline"))
  expect_error(get_data("unknown"), "model.*must be one of")
  expect_error(get_data(include_parameters = 1), "single logical")
  expect_error(get_data(include_parameters = c(TRUE, FALSE)), "single logical")
  expect_identical(get_parameters(), get_parameters(model = "baseline"))
  expect_identical(get_parameters(data = get_data("opaka")),
                   get_parameters(model = "opaka"))
  expect_error(get_parameters(data = 2), "model data list")
  expect_error(get_parameters(data = list()), "Could not infer")
  for (model in c("opaka", "bet")) {
    changed <- get_data(model)
    changed$label <- "annotated model"
    expect_identical(get_parameters(data = changed), get_parameters(model = model))
  }
})
