test_that("evaluate_priors returns zero for empty priors", {
  expect_equal(evaluate_priors(parameters = list(), priors = list()), 0)
  expect_equal(evaluate_priors(parameters = list(), priors = NULL), 0)
})

test_that("prior lambdas default to one and weight selected elements", {
  parameters <- list(par_sel = matrix(seq_len(6) / 10, 2, 3))
  prior <- list(par_sel = list(type = "normal", par1 = matrix(0, 2, 3),
                               par2 = 1, index = 1))
  base <- evaluate_priors(parameters, prior)
  expect_equal(evaluate_priors(parameters, within(prior, par_sel$lambda <- 1)), base)
  expect_equal(evaluate_priors(parameters, within(prior, par_sel$lambda <- 2)),
               2 * base)
  toggled <- opal_set_prior_lambda(prior, "par_sel", 0, element = 2)
  expected <- -sum(dnorm(parameters$par_sel[-2], 0, 1, log = TRUE))
  expect_equal(evaluate_priors(parameters, toggled), expected)
  expect_error(evaluate_priors(parameters,
    within(prior, par_sel$lambda <- c(1, NA, 1, 1, 1, 1))), "finite")
  expect_error(evaluate_priors(parameters,
    within(prior, par_sel$lambda <- rep(-1, 6))), "non-negative")
  expect_error(evaluate_priors(parameters,
    within(prior, par_sel$lambda <- 1:5)), "parameter length")
})

test_that("mapped-only prior resolution disables only fixed elements", {
  parameters <- list(log_B0 = c(1, 2, 3))
  priors <- list(log_B0 = list(type = "normal", par1 = 0, par2 = 1,
                               index = 1, lambda = c(1, 1, 1)))
  map <- list(log_B0 = factor(c(1, NA, 2)))
  resolved <- opal:::.opal_resolve_prior_lambdas(priors, map, parameters)
  expect_equal(resolved$log_B0$lambda, c(1, 0, 1))
  expect_equal(evaluate_priors(parameters, resolved),
               -sum(dnorm(parameters$log_B0[c(1, 3)], 0, 1, log = TRUE)))
})
