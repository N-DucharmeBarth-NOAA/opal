# Tests for dynamics module functions
library(RTMB)

# Test get_initial_numbers ----

test_that("get_initial_numbers with spawning_potential_a recovers B0 at equilibrium", {
  B0 <- 1e6
  h <- 0.8
  M_a <- rep(0.1, 10)
  maturity_a <- c(0, 0, 0, 0.5, 0.8, 1, 1, 1, 1, 1)
  weight_a <- seq(1, 10, length.out = 10)
  spawning_potential_a <- maturity_a * weight_a

  init <- get_initial_numbers(B0, h, M_a, spawning_potential_a)
  recovered_B0 <- sum(init$Ninit * spawning_potential_a)
  expect_equal(recovered_B0, B0, tolerance = 1e-6)
})

test_that("get_initial_numbers returns correct list elements", {
  B0 <- 1e6
  h <- 0.75
  M_a <- rep(0.2, 5)
  spawning_potential_a <- c(0, 0.2, 0.6, 1.0, 1.0)

  init <- get_initial_numbers(B0, h, M_a, spawning_potential_a)
  expect_named(init, c("Ninit", "R0", "alpha", "beta"))
  expect_equal(length(init$Ninit), length(M_a))
  expect_true(init$R0 > 0)
  expect_true(init$alpha > 0)
  expect_true(init$beta > 0)
})

test_that("get_initial_numbers BH parameters are internally consistent", {
  B0 <- 500000
  h <- 0.9
  M_a <- rep(0.15, 8)
  spawning_potential_a <- c(0, 0, 0.1, 0.4, 0.7, 1, 1, 1)

  init <- get_initial_numbers(B0, h, M_a, spawning_potential_a)
  # Check BH parameter derivation: alpha = 4*h*R0 / (5*h - 1)
  expect_equal(init$alpha, (4 * h * init$R0) / (5 * h - 1), tolerance = 1e-10)
  # Check beta: beta = B0 * (1 - h) / (5*h - 1)
  expect_equal(init$beta, (B0 * (1 - h)) / (5 * h - 1), tolerance = 1e-10)
})
