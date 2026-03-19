test_that("sex_ratio directly scales SPR0 and R0 appropriately", {
  library(RTMB)
  
  data <- get(data("opal_baseline_data", envir = environment()))
  params <- get(data("opal_baseline_parameters", envir = environment()))
  map <- get(data("opal_baseline_map", envir = environment()))
  
  obj_base <- MakeADFun(func = cmb(opal_model, data), parameters = params, map = map, silent = TRUE)
  rep_base <- obj_base$report()
  
  SPR0_base <- rep_base$B0 / rep_base$R0
  
  data_mod <- data
  data_mod$sex_ratio <- rep(0.5, data$n_age)
  
  obj_mod <- MakeADFun(func = cmb(opal_model, data_mod), parameters = params, map = map, silent = TRUE)
  rep_mod <- obj_mod$report()
  
  SPR0_mod <- rep_mod$B0 / rep_mod$R0
  
  # When B0 is fixed, halving sex_ratio should exactly halve SPR0 and double R0
  expect_equal(rep_mod$B0, rep_base$B0)
  expect_equal(SPR0_mod, SPR0_base * 0.5)
  expect_equal(rep_mod$R0, rep_base$R0 * 2)
})

test_that("sex_ratio length-based inputs resolve to age-based equivalent via PLA", {
  library(RTMB)
  
  data <- get(data("opal_baseline_data", envir = environment()))
  params <- get(data("opal_baseline_parameters", envir = environment()))
  map <- get(data("opal_baseline_map", envir = environment()))
  
  # Uniform ratio, but supplied as length-based vector (length = n_len)
  data_mod <- data
  data_mod$sex_ratio <- rep(0.5, data$n_len)
  
  obj_mod <- MakeADFun(func = cmb(opal_model, data_mod), parameters = params, map = map, silent = TRUE)
  rep_mod <- obj_mod$report()
  
  SPR0_mod <- rep_mod$B0 / rep_mod$R0
  
  # Baseline for comparison
  obj_base <- MakeADFun(func = cmb(opal_model, data), parameters = params, map = map, silent = TRUE)
  rep_base <- obj_base$report()
  SPR0_base <- rep_base$B0 / rep_base$R0
  
  expect_equal(SPR0_mod, SPR0_base * 0.5)
  expect_equal(rep_mod$R0, rep_base$R0 * 2)
})
