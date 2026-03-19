test_that("Initial Age Deviations correctly modify equilibrium without breaking R0", {
  M_a <- rep(0.2, 5)
  spawn_pot <- c(0, 0, 10, 20, 30)

  init_base <- get_initial_numbers(
    B0 = 1000,
    h = 0.8,
    M_a = M_a,
    spawning_potential_a = spawn_pot
  )

  init_test <- get_initial_numbers(
    B0 = 1000,
    h = 0.8,
    M_a = M_a,
    spawning_potential_a = spawn_pot,
    init_rdev_a = c(0, 0, 0.5, 0, 0),
    sigma_r = 0.5,
    init_bias_adj_a = rep(0.0, 5)
  )
  init_test_bias <- get_initial_numbers(
    B0 = 1000,
    h = 0.8,
    M_a = M_a,
    spawning_potential_a = spawn_pot,
    init_rdev_a = c(0, 0, 0.5, 0, 0),
    sigma_r = 0.5,
    init_bias_adj_a = c(0, 0, 1.0, 0, 0)
  )

  expect_equal(init_base$R0, init_test$R0)
  expect_equal(init_test$Ninit[3], init_base$Ninit[3] * exp(0.5))
  expect_equal(init_test$Ninit[c(1, 2, 4, 5)], init_base$Ninit[c(1, 2, 4, 5)])
  expect_equal(init_test_bias$Ninit[3], init_base$Ninit[3] * exp(0.5 - 1.0 * 0.5 * 0.5^2))
})
