data <- list(
  last_yr = 2022, age_increase_M = 25, length_m50 = 150, length_m95 = 180, 
  catch_UR_on = 0, catch_surf_case = 1, catch_LL1_case = 1, 
  scenarios_surf = data_csv1$scenarios_surface, 
  scenarios_LL1 = data_csv1$scenarios_LL1,
  removal_switch_f = c(0, 0, 0, 1, 0, 0), # 0=harvest rate, 1=direct removals
  sel_min_age_f = c(2, 2, 2, 8, 6, 0, 2),
  sel_max_age_f = c(17, 9, 17, 22, 25, 7, 17),
  sel_end_f = c(1, 0, 1, 1, 1, 0, 1),
  sel_LL1_yrs = c(1952, 1957, 1961, 1965, 1969, 1973, 1977, 1981, 1985, 1989, 
                  1993, 1997, 2001, 2006, 2007, 2008, 2011, 2014, 2017, 2020),
  sel_LL2_yrs = c(1969, 2001, 2005, 2008, 2011, 2014, 2017, 2020),
  sel_LL3_yrs = c(1954, 1961, 1965, 1969, 1970, 1971, 2005, 2006, 2007),
  sel_LL4_yrs = c(1953),
  sel_Ind_yrs = c(1976, 1995, 1997, 1999, 2002, 2004, 2006, 2008, 2010, 2012, 
                  2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022),
  sel_Aus_yrs = c(1952, 1969, 1973, 1977, 1981, 1985, 1989, 1993, 1997, 1998, 
                  1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 
                  2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 
                  2019, 2020, 2021, 2022),
  sel_CPUE_yrs = c(1969, 1973, 1977, 1981, 1985, 1989, 1993, 1997, 2001, 2006, 
                   2007, 2008, 2011, 2014, 2017, 2020),
  af_switch = 9,
  lf_switch = 9, lf_minbin = c(1, 1, 1, 11),
  cpue_switch = 1, cpue_a1 = 5, cpue_a2 = 17,
  troll_switch = 1, 
  aerial_switch = 4, aerial_tau = 0.3,
  tag_switch = 1, tag_var_factor = 1.82,
  hsp_switch = 1, hsp_false_negative = 0.7467647, 
  pop_switch = 1, 
  gt_switch = 1
)

data <- get_data(data_in = data)
parameters <- get_parameters(data = data)
data$priors <- get_priors(parameters = parameters)
map <- get_map(parameters = parameters)
obj <- MakeADFun(func = cmb(sbt_model, data), parameters = parameters, map = map)

phi_ya <- get_phi(log_psi = obj$report()$par_log_psi, 
                  data$length_m50, data$length_m95, 
                  data$length_mu_ysa, data$length_sd_a, 
                  data$dl_yal)
dimnames(phi_ya) <- list(year = data$first_yr:(data$last_yr + 1), age = data$age_a)

init <- get_initial_numbers(B0 = obj$report()$B0, h = obj$report()$par_h, 
                            M_a = obj$report()$M_a, phi_ya = phi_ya)


rec <- get_recruitment(sbio = obj$report()$spawning_biomass_y[2:93], rdev = obj$report()$rdev_y, 
                       B0 = obj$report()$B0, obj$report()$alpha, obj$report()$beta, obj$report()$sigma_r)

dyn0 <- do_dynamics(
  first_yr = data$first_yr, 
  first_yr_catch = data$first_yr_catch, 
  B0 = obj$report()$B0, R0 = obj$report()$R0,
  alpha = obj$report()$alpha, beta = obj$report()$beta, h = obj$report()$par_h, 
  sigma_r = obj$report()$sigma_r,
  rdev_y = obj$report()$rdev_y,
  M_a = obj$report()$M_a, 
  phi_ya = phi_ya, 
  init_number_a = init$Ninit,
  removal_switch_f = data$removal_switch_f, 
  catch_obs_ysf = data$catch_obs_ysf, 
  sel_fya = obj$report()$sel_fya,
  weight_fya = data$weight_fya, 
  af_sliced_ysfa = data$af_sliced_ysfa)

# plot(data$first_yr:(data$last_yr + 1), obj$report()$spawning_biomass_y, xlim = c(1931, max(proj_years)))
# points(data$first_yr:2023, dyn0$spawning_biomass_y, col = 2, pch = 3)

test_that("Check if we can recover the reconstruction using functions", {
  expect_equal(obj$report()$phi_ya, phi_ya)
  # expect_equal(obj$report()$recruitment_y[2:93], rec)
  expect_equal(obj$report()$number_ysa[1, 1,], init$Ninit)
  expect_equal(obj$report()$number_ysa, dyn0$number_ysa)
})

test_that("Check if it is OK to use a 2D AR1 for a single year", {
  sigma2 <- 0.06
  rho_y <- 0.7
  rho_a <- 0.9
  x <- matrix(rnorm(10), nrow = 1, ncol = 10)
  scale <- sqrt(sigma2) / sqrt(1 - rho_y^2) / sqrt(1 - rho_a^2)
  f1 <- function(x) dautoreg(x, phi = rho_y, log = TRUE)
  f2 <- function(x) dautoreg(x, phi = rho_a, log = TRUE)
  expect_equal(dseparable(f1, f2)(x, scale = scale), dautoreg(x, phi = rho_a, log = TRUE, scale = scale))
})

test_that("Check some priors", {
  parameters <- list(x = log(0.167), y = 0.167)
  priors <- list(
    x = list(type = "normal", par1 = log(0.12), par2 = 0.4, index = which("x" == names(parameters))),
    y = list(type = "student", par1 = 0.12, par2 = 0.4, index = which("y" == names(parameters)))
  )
  pr1 <- evaluate_priors(parameters, priors)
  pr2 <- -1 * sum(c(dnorm(parameters$x, log(0.12), 0.4, log = TRUE), 
                    # dt2(parameters$y, 0.12, 0.4, df = 3, log = TRUE)
                    dt((parameters$y - 0.12) / 0.4, df = 3, log = TRUE) - log(0.4)))
  expect_equal(pr1, pr2)
})
