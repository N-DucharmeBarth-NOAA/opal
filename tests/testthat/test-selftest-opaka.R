# Calibration: 2026-09-16, R 4.5.2, 100 truth-start simulations, seeds
# 20260917--20261016. The check_estimability() gradient gate (max <= 0.01)
# gave 99% convergence. Limits are max(abs(median) + 3 * bootstrap MC-SE,
# floor), rounded upward: B0 0.02, terminal biomass 0.03, terminal depletion
# 0.03, and maximum biomass error 0.11. The quickstart B0 prior is N(9, 1).

test_that("fitted opaka simulation plumbing provides active observations", {
  skip_if_not_selftest()

  inputs <- opaka_inputs()
  om <- opaka_obj(inputs)
  fit_quickstart(om, inputs$parameters)
  par <- om$env$last.par.best

  set.seed(1)
  first_simulation <- om$simulate(par = par)
  set.seed(1)
  second_simulation <- om$simulate(par = par)

  expect_true(all(c("cpue_log_obs", "lf_obs_flat") %in% names(first_simulation)))
  expect_identical(first_simulation, second_simulation)
  expect_false(isTRUE(all.equal(
    exp(first_simulation$cpue_log_obs), inputs$data$cpue_data$value
  )))
})

test_that("opaka quickstart conditional self-test records recovery calibration", {
  skip_if_not_selftest()

  n_sim <- as.integer(Sys.getenv("OPAL_SELFTEST_NSIM", "30"))
  start <- Sys.getenv("OPAL_SELFTEST_START", "truth")
  expect_true(is.finite(n_sim) && n_sim > 0)
  expect_true(start %in% c("truth", "default"))

  results <- run_selftest(n_sim, start)
  output_directory <- Sys.getenv("OPAL_SELFTEST_OUT", "")
  if (nzchar(output_directory)) {
    saveRDS(
      results,
      file.path(output_directory, sprintf("selftest-opaka-%s.rds", start))
    )
  }

  expect_equal(nrow(results), n_sim)
  expect_identical(
    names(results),
    c("sim", "converged", "max_gr", "re_B0", "re_sb_term", "re_dep_term",
      "max_abs_re_sb")
  )

  if (start == "truth") {
    converged <- results[results$converged, , drop = FALSE]
    expect_gte(nrow(converged) / n_sim, 0.90)
    expect_gt(nrow(converged), 0)
    expect_lte(abs(median(converged$re_B0)), 0.02)
    expect_lte(abs(median(converged$re_sb_term)), 0.03)
    expect_lte(abs(median(converged$re_dep_term)), 0.03)
    expect_lte(median(converged$max_abs_re_sb), 0.11)
  }
})