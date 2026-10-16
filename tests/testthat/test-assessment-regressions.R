test_that("CPUE decodes every flattened seasonal time step", {
  number <- array(c(10, 100, 1000, 20, 200, 2000), c(3, 2, 1))
  d <- data.frame(ts = 1:4, fishery = 1L, index = c(1L, 1L, 2L, 2L),
                  units = 2L, value = 1, se = 0.1)
  pars <- list(log_cpue_q = c(0, 0), log_cpue_tau = c(-2, -2),
               log_cpue_omega = c(0, 0), cpue_creep = c(0, 0))
  make <- function(d) RTMB::MakeADFun(function(p) sum(get_cpue_like(d, p,
    number, array(1, c(1, 2, 1)), array(1, c(1, 2, 1)))), pars, silent = TRUE)
  obj <- make(d)
  expect_equal(obj$report()$cpue_pred, c(10 + 1e-6, 20 + 1e-6,
    100 + 1e-6, 200 + 1e-6) / rep(c(15 + 1e-6, 150 + 1e-6), each = 2))
  expect_true(all(is.finite(obj$gr(obj$par))))
  for (bad in c(0, 1.5, 5, NA)) {
    d$ts[1] <- bad
    expect_error(make(d), "time step")
  }
})

test_that("bounds diagnostics identify boundary rows and preserve indices", {
  out <- check_bounds(list(par = c(a = 0, b = 1, c = 0.5)), rep(0, 3), rep(1, 3))
  expect_equal(out$par, c("a", "b"))
  expect_equal(out$index, 1:2)
  expect_equal(check_bounds(list(par = c(a = -1, a = 2)), c(-Inf, 0), c(0, Inf))$index,
               integer())
})

test_that("invalid priors cannot silently change the scientific target", {
  p <- list(theta = c(0.2, 0.4))
  prior <- list(theta = list(type = "normal", par1 = 0, par2 = 1, index = 1L))
  expect_equal(evaluate_priors(p, prior), -sum(stats::dnorm(p$theta, log = TRUE)))
  bad <- prior; bad$theta$type <- "normla"
  expect_error(evaluate_priors(p, bad), "distribution")
  bad <- prior; bad$theta$index <- 2L
  expect_error(evaluate_priors(p, bad), "index")
  bad <- prior; bad$theta$par2 <- -1
  expect_error(evaluate_priors(p, bad), "positive")
  bad <- prior; bad$theta$par1 <- 1:3
  expect_error(evaluate_priors(p, bad), "length")
  bad <- prior; bad$theta$type <- "beta"; bad$theta$par1 <- 2
  expect_error(evaluate_priors(p, bad), "between zero and one")
})

test_that("prior labels cannot silently point at a different parameter block", {
  expect_error(evaluate_priors(list(a=0,b=0),
    list(a=list(type="normal",par1=0,par2=1,index=2L))),"name and parameter index")
})
