initial_harvest_args <- function(n_season = 1L, F = c(1, 1),
                                 spawning = rep(1, 10)) {
  list(B0 = 1e6, h = 0.8, M_a = rep(0.2, 10),
       spawning_potential_a = spawning, init_F_f = F,
       sel_fa = matrix(1, length(F), 10), n_season = n_season)
}

test_that("overlapping initial fisheries cannot produce negative abundance", {
  for (seasons in c(1L, 2L, 4L)) {
    args <- initial_harvest_args(seasons, F = c(3, 3))
    init <- do.call(get_initial_numbers, args)
    expect_true(all(is.finite(init$Ninit)))
    expect_true(all(init$Ninit > 0))
    expect_gt(init$lp_penalty, 0)
    args$init_F_f <- c(0, 0)
    unfished <- do.call(get_initial_numbers, args)
    expect_equal(init$Ninit0, unfished$Ninit0)
    expect_equal(init$R0, unfished$R0)
  }
})

test_that("infeasible equilibrium recruitment is positive and penalised", {
  args <- initial_harvest_args(F = 3, spawning = c(0, 0, rep(1, 8)))
  init <- do.call(get_initial_numbers, args)
  expect_true(all(is.finite(init$Ninit)))
  expect_true(all(init$Ninit > 0))
  expect_gt(init$lp_penalty, 0)
})

test_that("feasible multi-fishery initialisation retains analytical equilibrium", {
  for (seasons in c(1L, 4L)) {
    args <- initial_harvest_args(seasons, F = c(0.05, 0.1))
    args$sel_fa <- rbind(plogis(seq_len(10) - 3), plogis(seq_len(10) - 5))
    init <- do.call(get_initial_numbers, args)
    u <- 1 - exp(-args$init_F_f / seasons)
    survival <- exp(-args$M_a) * (1 - colSums(args$sel_fa * u))^seasons
    per_recruit <- c(1, cumprod(survival[-length(survival)]))
    per_recruit[10] <- per_recruit[10] / (1 - survival[10])
    recruitment <- init$alpha - init$beta / sum(per_recruit)
    expect_equal(init$Ninit, recruitment * per_recruit, tolerance = 1e-12)
    expect_equal(init$lp_penalty, 0)
  }
})

test_that("initialisation AD tape crosses the harvest constraint correctly", {
  objective <- function(p) {
    args <- initial_harvest_args()
    args$init_F_f <- exp(p$log_F)
    init <- do.call(get_initial_numbers, args)
    init$lp_penalty + sum(init$Ninit) / 1e6
  }
  obj <- RTMB::MakeADFun(objective, list(log_F = log(c(0.1, 0.1))), silent = TRUE)
  for (F in list(c(0.1, 0.2), c(1, 1), c(3, 3))) {
    par <- log(F)
    fresh <- RTMB::MakeADFun(objective, list(log_F = par), silent = TRUE)
    expect_true(is.finite(obj$fn(par)))
    expect_true(all(is.finite(obj$gr(par))))
    expect_equal(obj$fn(par), fresh$fn(), tolerance = 1e-10)
    fd <- vapply(seq_along(par), function(i) {
      e <- replace(numeric(length(par)), i, 1e-5)
      (obj$fn(par + e) - obj$fn(par - e)) / 2e-5
    }, numeric(1))
    expect_equal(as.vector(obj$gr(par)), fd, tolerance = 1e-6)
  }
})

test_that("the full model includes the initialisation penalty", {
  d <- synth_data()
  d$catch_obs_ysf[] <- 0
  d$sel_fa_external <- matrix(1, d$n_fishery, d$n_age)
  d$maturity <- rep(1, d$n_age)
  d$fecundity <- rep(1, d$n_age)
  p <- synth_parameters(d)
  p$log_h <- log(0.8)
  p$log_init_F_f <- log(c(1, 1))
  obj <- synth_obj(d, p)
  report <- obj$report()
  expect_true(is.finite(obj$fn()))
  expect_true(all(is.finite(obj$gr())))
  expect_gt(report$lp_init_penalty, 0)
  expect_equal(report$lp_penalty, report$lp_init_penalty)
  expect_equal(obj$fn(), report$lp_prior + report$lp_penalty +
                 report$lp_rec + report$lp_init_rec + sum(report$lp_cpue) +
                 sum(report$lp_lf) + sum(report$lp_wf))
})
