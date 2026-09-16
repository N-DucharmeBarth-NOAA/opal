test_that("opakapaka selectivity types use opal codes", {
  inputs <- opaka_inputs()
  expect_identical(as.integer(inputs$data$sel_type_f), c(1L, 1L, 2L))
})

test_that("sex_ratio directly scales SPR0 and R0 appropriately", {
  inputs <- opaka_inputs()
  data_base <- inputs$data
  data_base$lf_switch <- 0L
  rep_base <- opaka_obj(inputs, data = data_base)$report()
  spr0_base <- rep_base$B0 / rep_base$R0

  data_mod <- data_base
  data_mod$sex_ratio <- rep(0.5, data_base$n_age)
  rep_mod <- opaka_obj(inputs, data = data_mod)$report()
  spr0_mod <- rep_mod$B0 / rep_mod$R0

  expect_equal(rep_mod$B0, rep_base$B0)
  expect_equal(spr0_mod, spr0_base * 0.5)
  expect_equal(rep_mod$R0, rep_base$R0 * 2)
})

test_that("sex_ratio length-based inputs resolve to age-based equivalent via PLA", {
  inputs <- opaka_inputs()
  data_base <- inputs$data
  data_base$lf_switch <- 0L
  rep_base <- opaka_obj(inputs, data = data_base)$report()
  spr0_base <- rep_base$B0 / rep_base$R0

  data_mod <- data_base
  data_mod$sex_ratio <- rep(0.5, data_base$n_len)
  rep_mod <- opaka_obj(inputs, data = data_mod)$report()
  spr0_mod <- rep_mod$B0 / rep_mod$R0

  expect_equal(spr0_mod, spr0_base * 0.5)
  expect_equal(rep_mod$R0, rep_base$R0 * 2)
})
