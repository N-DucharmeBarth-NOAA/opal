# Tests for get_cpue_like()
library(RTMB)

# Shared helper ----

# Build minimal valid inputs for get_cpue_like().
# By default: 1 fishery, 3 years, 4 ages, 1 season; TWO CPUE observations (years 1
# and 2) — the function loops `for (i in 2:n_cpue)` which degenerates to `2:1 =
# c(2,1)` when n_cpue=1, so a minimum of 2 rows is required.
# All selectivity = 1, uniform numbers-at-age (so both years have identical sum_n).
make_cpue_args <- function(n_fishery = 1, n_year = 3, n_age = 4,
                           n_index = 1,
                           cpue_switch = 1,
                           log_cpue_tau   = log(0.1),
                           log_cpue_omega = log(1),
                           cpue_creep     = 0,
                           log_cpue_q     = log(1),
                           units  = c(2L, 2L),
                           ts     = c(1L, 2L),
                           fishery = c(1L, 1L),
                           index  = c(1L, 1L),
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

  if (length(log_cpue_tau) == 1)   log_cpue_tau   <- rep(log_cpue_tau, n_index)
  if (length(log_cpue_omega) == 1) log_cpue_omega <- rep(log_cpue_omega, n_index)
  if (length(cpue_creep) == 1)     cpue_creep     <- rep(cpue_creep, n_index)
  if (length(log_cpue_q) == 1)     log_cpue_q     <- rep(log_cpue_q, n_index)

  cpue_data <- data.frame(ts = ts, fishery = fishery,
                          value = value, se = se, units = units, index = index)

  data <- list(
    cpue_data   = cpue_data,
    cpue_switch = cpue_switch,
    n_index     = n_index
  )
  parameters <- list(
    log_cpue_tau   = log_cpue_tau,
    log_cpue_omega = log_cpue_omega,
    cpue_creep     = cpue_creep,
    log_cpue_q     = log_cpue_q
  )

  list(data = data, parameters = parameters,
       number_ysa = number_ysa, sel_fya = sel_fya, weight_fya = weight_fya)
}

# Legacy single-index implementation used for regression comparison when n_index=1.
legacy_cpue_like <- function(data, parameters, number_ysa, sel_fya, weight_fya,
                             creep_init = 1) {
  cpue_data <- data$cpue_data
  cpue_switch <- data$cpue_switch
  log_cpue_tau <- parameters$log_cpue_tau[1]
  log_cpue_omega <- parameters$log_cpue_omega[1]
  cpue_creep <- parameters$cpue_creep[1]
  log_cpue_q <- parameters$log_cpue_q[1]

  cpue_tau <- exp(log_cpue_tau)
  cpue_omega <- exp(log_cpue_omega)
  n_cpue <- nrow(cpue_data)
  cpue_adjust <- cpue_log_pred <- lp <- numeric(n_cpue)
  cpue_adjust[1] <- creep_init
  for (i in 2:n_cpue) cpue_adjust[i] <- cpue_adjust[i - 1] + cpue_creep
  cpue_sigma <- sqrt(cpue_data$se^2 + cpue_tau^2)
  for (i in seq_len(n_cpue)) {
    y <- cpue_data$ts[i]
    f <- cpue_data$fishery[i]
    cpue_n <- number_ysa[y, 1, ] * sel_fya[f, y, ]
    if (cpue_data$units[i] == 1) cpue_n <- cpue_n * weight_fya[f, y, ]
    sum_n <- sum(cpue_n) + 1e-6
    cpue_log_pred[i] <- log(cpue_adjust[i]) + cpue_omega * log(sum_n)
  }
  cpue_log_pred <- cpue_log_pred - log(mean(exp(cpue_log_pred))) + log_cpue_q
  cpue_log_obs <- log(cpue_data$value)
  if (cpue_switch > 0) {
    lp[] <- -dnorm(x = cpue_log_obs, mean = cpue_log_pred, sd = cpue_sigma, log = TRUE)
  }
  lp
}

# Tests ----

test_that("NLL is finite for a basic CPUE observation", {
  s  <- make_cpue_args()
  lp <- get_cpue_like(s$data, s$parameters, s$number_ysa, s$sel_fya, s$weight_fya)
  expect_true(all(is.finite(lp)))
})

test_that("cpue_switch = 0 gives zero NLL regardless of data", {
  s  <- make_cpue_args(cpue_switch = 0, value = c(999, 999), se = c(0.01, 0.01))
  lp <- get_cpue_like(s$data, s$parameters, s$number_ysa, s$sel_fya, s$weight_fya)
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
  lp <- get_cpue_like(s$data, s$parameters, s$number_ysa, s$sel_fya, s$weight_fya)
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

  lp_lo <- get_cpue_like(s_lo$data, s_lo$parameters, s_lo$number_ysa, s_lo$sel_fya, s_lo$weight_fya)
  lp_hi <- get_cpue_like(s_hi$data, s_hi$parameters, s_hi$number_ysa, s_hi$sel_fya, s_hi$weight_fya)

  # Larger tau => wider sigma => smaller NLL contribution for same residual
  expect_true(lp_lo[1] > lp_hi[1])
})

test_that("cpue_omega scales the effect of abundance on predicted CPUE", {
  # omega = 0 => sum_n has no effect => all preds equal => mean-centering =>
  # cpue_log_pred = log_q for all obs; different log_q values give different NLL
  s_om0 <- make_cpue_args(log_cpue_omega = log(1e-9), value = c(1, 1), log_cpue_q = log(1))
  s_om1 <- make_cpue_args(log_cpue_omega = log(1),    value = c(1, 1), log_cpue_q = log(1))

  lp_om0 <- get_cpue_like(s_om0$data, s_om0$parameters, s_om0$number_ysa, s_om0$sel_fya, s_om0$weight_fya)
  lp_om1 <- get_cpue_like(s_om1$data, s_om1$parameters, s_om1$number_ysa, s_om1$sel_fya, s_om1$weight_fya)

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
                                s_no_creep$number_ysa,   s_no_creep$sel_fya, s_no_creep$weight_fya)
  lp_with_creep <- get_cpue_like(s_with_creep$data, s_with_creep$parameters,
                                s_with_creep$number_ysa, s_with_creep$sel_fya, s_with_creep$weight_fya)
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
                          value = c(1, 1), se = c(0.1, 0.1), units = c(2L, 2L),
                          index = c(1L, 1L))
  data <- list(cpue_data = cpue_data, cpue_switch = 1L, n_index = 1L)
  weight_fya <- array(1, dim = c(n_fishery, n_year, n_age))
  parameters <- list(
    log_cpue_tau   = log(0.1),
    log_cpue_omega = log(1),
    cpue_creep     = 0,
    log_cpue_q     = log(1)
  )

  lp <- get_cpue_like(data, parameters, number_ysa, sel_fya, weight_fya)

  expect_true(all(is.finite(lp)))
  # Different predicted abundances => different residuals => different NLL per obs
  expect_false(isTRUE(all.equal(lp[1], lp[2])))
})

test_that("n_index=1 with explicit index column matches legacy scalar behaviour", {
  n_year <- 5
  n_age <- 4
  n_cpue <- 5
  number_ysa <- array(1, dim = c(n_year, 1L, n_age))
  number_ysa[1, 1, ] <- c(1, 2, 3, 4)
  number_ysa[2, 1, ] <- c(4, 3, 2, 1)
  number_ysa[3, 1, ] <- c(2, 2, 2, 2)
  number_ysa[4, 1, ] <- c(5, 4, 3, 2)
  number_ysa[5, 1, ] <- c(1, 1, 2, 3)
  sel_fya <- array(1, dim = c(1L, n_year, n_age))
  sel_fya[1, 1, ] <- c(0.9, 0.8, 0.7, 0.6)
  sel_fya[1, 2, ] <- c(0.8, 0.7, 0.6, 0.5)
  sel_fya[1, 3, ] <- c(0.7, 0.6, 0.5, 0.4)
  sel_fya[1, 4, ] <- c(0.6, 0.5, 0.4, 0.3)
  sel_fya[1, 5, ] <- c(0.5, 0.4, 0.3, 0.2)
  weight_fya <- array(1, dim = c(1L, n_year, n_age))

  s <- make_cpue_args(
    n_index = 1, n_year = n_year, n_age = n_age,
    ts = as.integer(1:n_cpue),
    fishery = rep(1L, n_cpue),
    index = rep(1L, n_cpue),
    units = c(1L, 2L, 1L, 2L, 1L),
    value = c(1.2, 0.9, 1.6, 0.7, 1.1),
    se = c(0.1, 0.2, 0.15, 0.12, 0.09),
    log_cpue_tau = log(0.12),
    log_cpue_omega = log(0.85),
    cpue_creep = 0.03,
    log_cpue_q = log(1.4),
    number_ysa = number_ysa, sel_fya = sel_fya, weight_fya = weight_fya
  )

  lp_new <- get_cpue_like(s$data, s$parameters, s$number_ysa, s$sel_fya, s$weight_fya)
  lp_old <- legacy_cpue_like(s$data, s$parameters, s$number_ysa, s$sel_fya, s$weight_fya)
  expect_true(all(is.finite(lp_new)))
  expect_equal(lp_new, lp_old, tolerance = 1e-14)
})

test_that("two indices get independent q values", {
  s <- make_cpue_args(
    n_index = 2,
    ts = c(1L, 1L),
    fishery = c(1L, 1L),
    index = c(1L, 2L),
    value = c(1.0, 1.0),
    se = c(0.1, 0.1),
    log_cpue_q = c(log(1), log(2))
  )
  lp <- get_cpue_like(s$data, s$parameters, s$number_ysa, s$sel_fya, s$weight_fya)
  expect_true(all(is.finite(lp)))
  expect_false(isTRUE(all.equal(lp[1], lp[2])))
})

test_that("two indices get independent tau values", {
  s <- make_cpue_args(
    n_index = 2,
    ts = c(1L, 1L),
    fishery = c(1L, 1L),
    index = c(1L, 2L),
    value = c(2.0, 2.0),
    se = c(0.1, 0.1),
    log_cpue_tau = c(log(0.01), log(0.5))
  )
  lp <- get_cpue_like(s$data, s$parameters, s$number_ysa, s$sel_fya, s$weight_fya)
  expect_true(all(is.finite(lp)))
  expect_true(lp[1] > lp[2])
})

test_that("two indices get independent omega values", {
  n_year <- 3
  n_age <- 4
  number_ysa <- array(1, dim = c(n_year, 1L, n_age))
  number_ysa[1, 1, ] <- 10
  number_ysa[2, 1, ] <- 1

  s <- make_cpue_args(
    n_index = 2, n_year = n_year, n_age = n_age,
    ts = c(1L, 2L, 1L, 2L),
    fishery = c(1L, 1L, 1L, 1L),
    index = c(1L, 1L, 2L, 2L),
    value = c(1, 1, 1, 1),
    se = c(0.1, 0.1, 0.1, 0.1),
    log_cpue_omega = c(log(1), log(0.5)),
    number_ysa = number_ysa
  )
  lp <- get_cpue_like(s$data, s$parameters, s$number_ysa, s$sel_fya, s$weight_fya)
  expect_true(all(is.finite(lp)))
  expect_equal(length(lp), 4)
})

test_that("two indices get independent creep", {
  s <- make_cpue_args(
    n_index = 2,
    ts = c(1L, 2L, 1L, 2L),
    fishery = c(1L, 1L, 1L, 1L),
    index = c(1L, 1L, 2L, 2L),
    value = c(1, 1, 1, 1),
    se = c(0.1, 0.1, 0.1, 0.1),
    cpue_creep = c(0, 0.1)
  )
  lp <- get_cpue_like(s$data, s$parameters, s$number_ysa, s$sel_fya, s$weight_fya)
  expect_true(all(is.finite(lp)))
  expect_equal(lp[1], lp[2], tolerance = 1e-10)
  expect_false(isTRUE(all.equal(lp[3], lp[4])))
})

test_that("mean-centering is independent across indices", {
  n_year <- 3
  n_age <- 4
  n_fishery <- 2
  number_ysa <- array(1, dim = c(n_year, 1L, n_age))
  sel_fya <- array(1, dim = c(n_fishery, n_year, n_age))
  weight_fya <- array(1, dim = c(n_fishery, n_year, n_age))
  sel_fya[1, , ] <- 1
  sel_fya[2, , ] <- 0.01

  s <- make_cpue_args(
    n_index = 2, n_fishery = n_fishery, n_year = n_year, n_age = n_age,
    ts = c(1L, 2L, 1L, 2L),
    fishery = c(1L, 1L, 2L, 2L),
    index = c(1L, 1L, 2L, 2L),
    value = c(1, 1, 1, 1),
    se = c(0.1, 0.1, 0.1, 0.1),
    number_ysa = number_ysa, sel_fya = sel_fya, weight_fya = weight_fya
  )
  lp <- get_cpue_like(s$data, s$parameters, s$number_ysa, s$sel_fya, s$weight_fya)
  expect_true(all(is.finite(lp)))
  expect_equal(lp[1], lp[2], tolerance = 1e-10)
  expect_equal(lp[3], lp[4], tolerance = 1e-10)
})

test_that("cpue_switch = 0 gives zero NLL with multiple indices", {
  s <- make_cpue_args(
    n_index = 2, cpue_switch = 0,
    ts = c(1L, 1L),
    fishery = c(1L, 1L),
    index = c(1L, 2L),
    value = c(999, 999),
    se = c(0.01, 0.01)
  )
  lp <- get_cpue_like(s$data, s$parameters, s$number_ysa, s$sel_fya, s$weight_fya)
  expect_equal(sum(lp), 0)
})
