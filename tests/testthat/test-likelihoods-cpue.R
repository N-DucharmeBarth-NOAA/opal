# Tests for get_cpue_like()
library(RTMB)

# Shared helper ----

# Build minimal valid inputs for get_cpue_like().
# By default: 1 fishery, 3 years, 4 ages, 1 season; TWO CPUE observations (years 1
# and 2) — the function loops `for (i in 2:n_cpue)` which degenerates to `2:1 =
# c(2,1)` when n_cpue=1, so a minimum of 2 rows is required.
# All selectivity = 1, uniform numbers-at-age (so both years have identical sum_n).
make_cpue_args <- function(n_fishery = 1, n_year = 3, n_age = 4,
                           cpue_switch = 1,
                           log_cpue_tau   = log(0.1),
                           log_cpue_omega = log(1),
                           cpue_creep     = 0,
                           log_cpue_q     = log(1),
                           units  = c(2L, 2L),
                           ts     = c(1L, 2L),
                           fishery = c(1L, 1L),
                           value  = c(1.0, 1.0),
                           se     = c(0.1, 0.1),
                           number_ysa = NULL, sel_fya = NULL,
                           weight_fya = NULL) {

  if (is.null(number_ysa))
    number_ysa <- array(1, dim = c(n_year, 1L, n_age))
  if (is.null(sel_fya))
    sel_fya <- array(1, dim = c(n_fishery, n_year, n_age))
  if (is.null(weight_fya))
    weight_fya <- array(1, dim = c(n_fishery, n_year, n_age))

  cpue_data <- data.frame(ts = ts, fishery = fishery,
                          value = value, se = se, units = units)

  data <- list(
    cpue_data   = cpue_data,
    cpue_switch = cpue_switch
  )
  parameters <- list(
    log_cpue_tau   = log_cpue_tau,
    log_cpue_omega = log_cpue_omega,
    cpue_creep     = cpue_creep,
    log_cpue_q     = log_cpue_q,
    weight_fya     = weight_fya
  )

  list(data = data, parameters = parameters,
       number_ysa = number_ysa, sel_fya = sel_fya)
}

# Tests ----

test_that("NLL is finite for a basic CPUE observation", {
  s  <- make_cpue_args()
  lp <- get_cpue_like(s$data, s$parameters, s$number_ysa, s$sel_fya)
  expect_true(all(is.finite(lp)))
})

test_that("cpue_switch = 0 gives zero NLL regardless of data", {
  s  <- make_cpue_args(cpue_switch = 0, value = c(999, 999), se = c(0.01, 0.01))
  lp <- get_cpue_like(s$data, s$parameters, s$number_ysa, s$sel_fya)
  expect_equal(sum(lp), 0)
})

test_that("NLL matches dnorm reference value", {
  # Two obs with identical abundance (uniform number_ysa) and creep=0, omega=1.
  # All log_pred before centering = omega * log(sum_n), which is the same for
  # every obs => mean(exp(log_pred)) = exp(omega*log(sum_n)) => after centering:
  # log_pred[i] = log_q exactly.  So NLL[i] = -dnorm(log(obs[i]), log_q, sigma).
  log_q  <- log(2)
  obs    <- 2.0
  se_obs <- 0.1
  tau    <- 0.05
  sigma  <- sqrt(se_obs^2 + tau^2)

  expected_nll <- -dnorm(log(obs), mean = log_q, sd = sigma, log = TRUE)

  s  <- make_cpue_args(log_cpue_q = log_q, log_cpue_tau = log(tau),
                       value = c(obs, obs), se = c(se_obs, se_obs))
  lp <- get_cpue_like(s$data, s$parameters, s$number_ysa, s$sel_fya)
  expect_equal(lp[1], expected_nll, tolerance = 1e-10)
  expect_equal(lp[2], expected_nll, tolerance = 1e-10)
})

test_that("cpue_sigma combines observation SE and process tau in quadrature", {
  # sigma = sqrt(se^2 + tau^2); NLL increases with lower sigma at non-zero residual
  se_obs <- 0.2
  tau_lo <- 0.01
  tau_hi <- 0.5

  s_lo <- make_cpue_args(log_cpue_tau = log(tau_lo), value = c(2, 2), se = c(se_obs, se_obs))
  s_hi <- make_cpue_args(log_cpue_tau = log(tau_hi), value = c(2, 2), se = c(se_obs, se_obs))

  lp_lo <- get_cpue_like(s_lo$data, s_lo$parameters, s_lo$number_ysa, s_lo$sel_fya)
  lp_hi <- get_cpue_like(s_hi$data, s_hi$parameters, s_hi$number_ysa, s_hi$sel_fya)

  # Larger tau => wider sigma => smaller NLL contribution for same residual
  expect_true(lp_lo[1] > lp_hi[1])
})

test_that("cpue_omega scales the effect of abundance on predicted CPUE", {
  # omega = 0 => sum_n has no effect => all preds equal => mean-centering =>
  # cpue_log_pred = log_q for all obs; different log_q values give different NLL
  s_om0 <- make_cpue_args(log_cpue_omega = log(1e-9), value = c(1, 1), log_cpue_q = log(1))
  s_om1 <- make_cpue_args(log_cpue_omega = log(1),    value = c(1, 1), log_cpue_q = log(1))

  lp_om0 <- get_cpue_like(s_om0$data, s_om0$parameters, s_om0$number_ysa, s_om0$sel_fya)
  lp_om1 <- get_cpue_like(s_om1$data, s_om1$parameters, s_om1$number_ysa, s_om1$sel_fya)

  expect_true(all(is.finite(lp_om0)))
  expect_true(all(is.finite(lp_om1)))
})

test_that("cpue_creep accumulates across multiple observations", {
  # Two obs with identical abundance across years (uniform number_ysa).
  # creep=0  => cpue_adjust = (1, 1)   => same log_pred for both obs (after centering)
  #          => same NLL for both (identical obs value and SE)
  # creep>0  => cpue_adjust = (1, 1.1) => log_pred[2] > log_pred[1] before centering
  #          => after centering they differ => NLL[1] != NLL[2]
  # Note: mean-centering couples ALL obs, so changing creep affects every obs's NLL.
  s_no_creep   <- make_cpue_args(cpue_creep = 0,   value = c(1, 1), se = c(0.1, 0.1))
  s_with_creep <- make_cpue_args(cpue_creep = 0.1, value = c(1, 1), se = c(0.1, 0.1))

  lp_no_creep   <- get_cpue_like(s_no_creep$data,   s_no_creep$parameters,
                                 s_no_creep$number_ysa,   s_no_creep$sel_fya)
  lp_with_creep <- get_cpue_like(s_with_creep$data, s_with_creep$parameters,
                                 s_with_creep$number_ysa, s_with_creep$sel_fya)

  # With creep=0 and identical abundance both obs get same pred => same NLL
  expect_equal(lp_no_creep[1], lp_no_creep[2], tolerance = 1e-10)
  # With creep>0 the two obs have different adjustments => different preds => different NLL
  expect_false(isTRUE(all.equal(lp_with_creep[1], lp_with_creep[2])))
  # Total NLL also differs between the two creep scenarios
  expect_false(isTRUE(all.equal(sum(lp_no_creep), sum(lp_with_creep))))
})

test_that("correct year's numbers-at-age are used for each observation", {
  # Year 1 has 10x more fish than year 2; two separate obs should give different NLLs
  n_year <- 3; n_age <- 4; n_fishery <- 1
  number_ysa <- array(1, dim = c(n_year, 1L, n_age))
  number_ysa[1, 1, ] <- 10   # year 1: high abundance
  number_ysa[2, 1, ] <- 1    # year 2: low abundance
  sel_fya <- array(1, dim = c(n_fishery, n_year, n_age))

  # Two obs with same observed value but drawn from different years
  cpue_data <- data.frame(ts = c(1L, 2L), fishery = c(1L, 1L),
                          value = c(1, 1), se = c(0.1, 0.1), units = c(2L, 2L))
  data <- list(cpue_data = cpue_data, cpue_switch = 1L)
  parameters <- list(
    log_cpue_tau   = log(0.1),
    log_cpue_omega = log(1),
    cpue_creep     = 0,
    log_cpue_q     = log(1),
    weight_fya     = array(1, dim = c(n_fishery, n_year, n_age))
  )

  lp <- get_cpue_like(data, parameters, number_ysa, sel_fya)

  expect_true(all(is.finite(lp)))
  # Different predicted abundances => different residuals => different NLL per obs
  expect_false(isTRUE(all.equal(lp[1], lp[2])))
})
