test_that("opaka quickstart objective simulates active observations", {
  inputs <- opaka_inputs()
  obj <- opaka_obj(inputs)

  expect_true(is.finite(obj$fn()))

  set.seed(42)
  sim <- obj$simulate()

  expect_true(all(c("cpue_log_obs", "lf_obs_flat") %in% names(sim)))
  expect_length(sim$cpue_log_obs, nrow(inputs$data$cpue_data))
  expect_length(sim$lf_obs_flat, length(inputs$data$lf_obs_flat))
  expect_true(all(is.finite(sim$cpue_log_obs)))
  expect_true(all(is.finite(sim$lf_obs_flat)))
  expect_true(all(sim$lf_obs_flat >= 0))

  offset <- 0L
  for (fishery_index in seq_along(inputs$data$lf_fishery_f)) {
    fishery <- inputs$data$lf_fishery_f[fishery_index]
    n_bin <- inputs$data$lf_maxbin[fishery] -
      inputs$data$lf_minbin[fishery] + 1L
    for (observation_index in seq_len(inputs$data$lf_n_f[fishery_index])) {
      indices <- offset + seq_len(n_bin)
      expect_lt(abs(
        sum(sim$lf_obs_flat[indices]) - sum(inputs$data$lf_obs_flat[indices])
      ), 1)
      offset <- offset + n_bin
    }
  }

  set.seed(42)
  expect_identical(sim, obj$simulate())
  expect_false(isTRUE(all.equal(
    exp(sim$cpue_log_obs), inputs$data$cpue_data$value
  )))
})