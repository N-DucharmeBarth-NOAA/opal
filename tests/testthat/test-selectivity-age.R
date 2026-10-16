# Tests for age-based selectivity (sel_basis_f = 2)
library(RTMB)

age_test_inputs <- function() {
  ages <- 0:19
  mu_a <- 20 + (90 - 20) * (1 - exp(-0.25 * ages))
  len_lower <- seq(10, 118, by = 2)
  len_upper <- len_lower + 2
  list(
    ages = ages, mu_a = mu_a,
    len_lower = len_lower, len_upper = len_upper,
    len_mid = len_lower + 1,
    pla = get_pla(len_lower, len_upper, mu_a, 0.1 * mu_a)
  )
}

age_test_par_sel <- function() {
  rbind(c(-0.3, -0.5, 0, 0, 0, 0),
        c(-0.5, -1, 0, 0.5, -9, -2),
        c(-0.5, -0.5, 0.2, 0.3, -0.5, 0.1))
}

test_that("age-based selectivity evaluates each curve family on age", {
  x <- age_test_inputs()
  par_sel <- age_test_par_sel()
  data <- list(n_fishery = 3, n_year = 4, n_age = length(x$ages),
               min_age = 0, sel_type_f = c(1L, 2L, 3L),
               sel_basis_f = c(2L, 2L, 2L))
  sel <- get_selectivity(data, par_sel, x$pla, x$len_mid, age_a = x$ages)
  expect_equal(dim(sel), c(3, 4, length(x$ages)))
  for (f in 1:3) {
    expected <- sel_length(x$ages, par_sel[f, ], data$sel_type_f[f])
    for (y in 1:4) expect_equal(sel[f, y, ], expected)
  }
  expect_true(all(sel >= 0 & sel <= 1))

  # age_a defaults to seq(min_age, length.out = n_age)
  expect_equal(get_selectivity(data, par_sel, x$pla, x$len_mid), sel)
  data$min_age <- NULL
  expect_equal(
    get_selectivity(data, par_sel, x$pla, x$len_mid)[1, 1, ],
    sel_logistic(seq_along(x$ages), par_sel[1, ])
  )
})

test_that("explicit length basis matches omitting sel_basis_f", {
  x <- age_test_inputs()
  par_sel <- age_test_par_sel()
  data <- list(n_fishery = 3, n_year = 2, n_age = length(x$ages),
               sel_type_f = c(1L, 2L, 3L))
  sel_default <- get_selectivity(data, par_sel, x$pla, x$len_mid)
  data$sel_basis_f <- rep(1L, 3)
  expect_identical(get_selectivity(data, par_sel, x$pla, x$len_mid), sel_default)
})

test_that("mixed bases give the right values for each fishery", {
  x <- age_test_inputs()
  par_sel <- age_test_par_sel()
  data <- list(n_fishery = 3, n_year = 2, n_age = length(x$ages),
               sel_type_f = c(1L, 2L, 3L), sel_basis_f = c(2L, 1L, 2L))
  sel <- get_selectivity(data, par_sel, x$pla, x$len_mid, age_a = x$ages)
  expect_equal(sel[1, 2, ], sel_logistic(x$ages, par_sel[1, ]))
  expect_equal(sel[2, 2, ],
               c(t(x$pla) %*% sel_double_normal(x$len_mid, par_sel[2, ])))
  expect_equal(sel[3, 2, ], sel_double_richards(x$ages, par_sel[3, ]))
})

test_that("age-based selectivity does not depend on growth", {
  x <- age_test_inputs()
  par_sel <- age_test_par_sel()
  data <- list(n_fishery = 2, n_year = 2, n_age = length(x$ages),
               sel_type_f = c(1L, 1L), sel_basis_f = c(2L, 1L))
  pla2 <- get_pla(x$len_lower, x$len_upper, 0.8 * x$mu_a, 0.2 * x$mu_a)
  sel1 <- get_selectivity(data, par_sel[1:2, ], x$pla, x$len_mid, x$ages)
  sel2 <- get_selectivity(data, par_sel[1:2, ], pla2, x$len_mid, x$ages)
  expect_identical(sel1[1, , ], sel2[1, , ])
  expect_false(isTRUE(all.equal(sel1[2, , ], sel2[2, , ])))

  # Same check through the model, changing growth parameters
  d <- synth_full_data(wf_switch = 0L)
  d$sel_basis_f <- c(2L, 1L)
  p <- synth_parameters(d)
  p$par_sel[1, ] <- c(-0.3, -0.5, 0, 0, 0, 0)
  p$par_sel[2, ] <- c(0, 0, 0, 0, 0, 0)
  rep1 <- synth_obj(d, p)$report()
  p$log_L1 <- log(25)
  p$log_CV1 <- log(0.2)
  rep2 <- synth_obj(d, p)$report()
  expect_identical(rep1$sel_fya[1, , ], rep2$sel_fya[1, , ])
  expect_equal(rep1$sel_fya[1, 1, ], sel_logistic(1:5, p$par_sel[1, ]))
  expect_false(isTRUE(all.equal(rep1$sel_fya[2, , ], rep2$sel_fya[2, , ])))
})

test_that("sel_basis_f is validated", {
  x <- age_test_inputs()
  par_sel <- age_test_par_sel()
  data <- list(n_fishery = 3, n_year = 2, n_age = length(x$ages),
               sel_type_f = c(1L, 2L, 3L))
  data$sel_basis_f <- c(1L, 2L)
  expect_error(get_selectivity(data, par_sel, x$pla, x$len_mid),
               "'sel_basis_f' must have length n_fishery \\(3\\)")
  data$sel_basis_f <- c(1L, 3L, 2L)
  expect_error(get_selectivity(data, par_sel, x$pla, x$len_mid),
               "only 1 \\(length\\) or 2 \\(age\\); found: 3")
  data$sel_basis_f <- c(1L, NA, 2L)
  expect_error(get_selectivity(data, par_sel, x$pla, x$len_mid),
               "only 1 \\(length\\) or 2 \\(age\\)")
  data$sel_basis_f <- c(1L, 2L, 2L)
  expect_error(get_selectivity(data, par_sel, x$pla, x$len_mid, age_a = 1:3),
               "'age_a' must have length n_age")

  d <- synth_full_data(wf_switch = 0L)
  d$sel_basis_f <- c(2L, 0L)
  expect_error(synth_obj(d), "'sel_basis_f' must contain only")
})

test_that("sel_fa_external takes precedence over sel_basis_f", {
  d <- synth_full_data(wf_switch = 0L)
  d$sel_basis_f <- c(2L, 2L)
  d$sel_fa_external <- matrix(c(0.2, 0.5, 0.8, 1, 1), 2, 5, byrow = TRUE)
  sel <- synth_obj(d)$report()$sel_fya
  expect_equal(sel[1, 2, ], d$sel_fa_external[1, ])
})

test_that("age-based selectivity can be mirrored through the map", {
  d <- synth_full_data(wf_switch = 0L)
  d$sel_basis_f <- c(2L, 2L)
  p <- synth_parameters(d)
  p$par_sel[1, ] <- c(-0.3, -0.5, 0, 0, 0, 0)
  p$par_sel[2, ] <- p$par_sel[1, ]
  map <- synth_map(p)
  lv <- matrix(NA_integer_, 2, 6)
  lv[, 1:2] <- matrix(1:2, 2, 2, byrow = TRUE)
  map$par_sel <- factor(as.vector(lv))
  obj <- synth_obj(d, p, map)
  is_sel <- names(obj$par) == "par_sel"
  expect_equal(sum(is_sel), 2L)
  par <- obj$par
  par[is_sel] <- c(0.1, -0.2)
  sel <- obj$report(par)$sel_fya
  expect_identical(sel[1, , ], sel[2, , ])
  expect_equal(sel[1, 1, ], sel_logistic(1:5, c(0.1, -0.2)))
})

test_that("age-based logistic selectivity is recovered from synthetic data", {
  d <- synth_full_data(wf_switch = 0L)
  d$sel_basis_f <- c(2L, 1L)
  p <- synth_parameters(d)
  truth <- c(-0.4, -0.3)
  p$par_sel[1, 1:2] <- truth
  map <- synth_map(p)
  lv <- matrix(NA_integer_, 2, 6)
  lv[1, 1:2] <- 1:2
  map$par_sel <- factor(as.vector(lv))

  # Expected-value length compositions generated from the true selectivity
  n_obs <- 1000
  obj_true <- synth_obj(d, p, map)
  rep_true <- obj_true$report()
  pred <- rep_true$lf_pred[[1]]
  d$lf_obs_flat <- as.vector(t(pred)) * n_obs
  d$lf_obs <- d$lf_obs_flat
  d$lf_n <- rep(n_obs, d$n_lf)
  d$lf_obs_prop <- as.vector(t(pred))

  p_start <- p
  p_start$par_sel[1, 1:2] <- c(0.3, 0.3)
  obj <- synth_obj(d, p_start, map)
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  expect_equal(opt$convergence, 0L)
  sel_hat <- obj$report(opt$par)$sel_fya[1, 1, ]
  expect_equal(sel_hat, rep_true$sel_fya[1, 1, ], tolerance = 1e-3)
  expect_equal(unname(opt$par[names(opt$par) == "par_sel"]), truth,
               tolerance = 1e-2)
})

test_that("SS3 converters handle age-based patterns 12 and 20", {
  ages <- 0:28
  len <- seq(11, 199, by = 2)
  ss3 <- rbind(c(4.5, 2.5, NA, NA, NA, NA),
               c(3, -2, 1.2, 2.5, -9, -3),
               c(100, 10, NA, NA, NA, NA))
  type <- c(1L, 2L, 1L)
  basis <- c("age", "age", "length")
  par_sel <- convert_ss3_selex_to_rtmb(ss3, type, len, basis = basis,
                                       sel_ages = ages)
  expect_equal(par_sel[1, 1:2],
               c((4.5 - mean(ages)) / sd(ages), log(2.5 / sd(ages))))
  expect_equal(par_sel[3, 1:2],
               c((100 - mean(len)) / sd(len), log(10 / sd(len))))
  expect_equal(par_sel[2, 3], 1.2 - 2 * log(sd(ages)))

  # The age-based logistic reproduces SS3's inflection and width on age
  expect_equal(sel_logistic(ages, par_sel[1, ]),
               1 / (1 + exp(-log(19) * (ages - 4.5) / 2.5)))

  back <- convert_rtmb_selex_to_ss3(par_sel, type, len, basis = basis,
                                    sel_ages = ages)
  expect_equal(back$peak_or_inflection, ss3[, 1])
  expect_equal(back$top_logit_or_width, ss3[, 2])
  expect_equal(back$ascend_se[2], 1.2)
  expect_equal(back$descend_se[2], 2.5)

  # Default basis is unchanged
  expect_identical(convert_ss3_selex_to_rtmb(ss3[3, , drop = FALSE], 1L, len),
                   convert_ss3_selex_to_rtmb(ss3[3, , drop = FALSE], 1L, len,
                                             basis = "length"))
  expect_error(convert_ss3_selex_to_rtmb(ss3, type, len, basis = basis),
               "'sel_ages' must be supplied")
  expect_error(convert_ss3_selex_to_rtmb(ss3, type, len, basis = "weight"),
               "'basis' must be")
  expect_error(convert_ss3_selex_to_rtmb(ss3, type, len, basis = c("age", "length"),
                                         sel_ages = ages),
               "'basis' must have length 1 or n_fishery")
})

test_that("explicit length basis leaves the opaka model unchanged", {
  inputs <- opaka_inputs()
  obj0 <- opaka_obj(inputs)
  data <- inputs$data
  data$sel_basis_f <- rep(1L, data$n_fishery)
  obj1 <- opaka_obj(inputs, data = data)
  expect_identical(obj1$fn(obj0$par), obj0$fn(obj0$par))
  expect_identical(obj1$report(obj0$par), obj0$report(obj0$par))
})
