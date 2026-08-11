test_that("evaluate_priors returns zero for empty priors", {
  expect_equal(evaluate_priors(parameters = list(), priors = list()), 0)
  expect_equal(evaluate_priors(parameters = list(), priors = NULL), 0)
})
