test_that("conditional composition CDFs match independent analytical cases", {
  pit <- opal:::.opal_composition_pit
  set.seed(42); u <- runif(1)
  set.seed(42); z <- pit(c(2, 3, 5), c(.2, .3, .5), 1)
  expect_equal(z[1], qnorm(pbinom(1, 10, .2) + u * dbinom(2, 10, .2)))
  expect_true(is.na(z[3]))
  z <- pit(c(.2, .3, .5), c(.3, .3, .4), 2, 10)
  expect_equal(z[1:2], qnorm(c(pbeta(.2, 3, 7), pbeta(.3/.8, 3, 4))))
  set.seed(42); z <- pit(c(2, 3, 5), c(.2, .3, .5), 3, 10)
  mass <- function(k) choose(10, k) * beta(k + 2, 10 - k + 8) / beta(2, 8)
  expect_equal(z[1], qnorm(mass(0) + mass(1) + u * mass(2)))
  expect_true(all(is.na(pit(c(10, 0, 0), c(.2, .3, .5), 1)[2:3])))
  expect_true(all(is.finite(pit(c(1000, 1, 1), c(.01, .49, .5), 1)[1:2])))
})

test_that("composition residuals calibrate for all three distributions", {
  set.seed(704)
  n <- 1500L
  p <- c(.2, .3, .5)
  for (family in 1:3) {
    z <- replicate(n, {
      obs <- if (family == 1) as.numeric(rmultinom(1, 60, p)) else {
        q <- rgamma(3, p * 30); q <- q / sum(q)
        if (family == 2) q else as.numeric(rmultinom(1, 60, q))
      }
      opal:::.opal_composition_pit(obs, p, family, 30)[1:2]
    })
    expect_lt(abs(mean(z)), .08)
    expect_lt(abs(sd(as.numeric(z)) - 1), .06)
    expect_lt(abs(cor(t(z))[1, 2]), .08)
  }
})

osa_fixture <- function(family) {
  x <- small_opal_object()
  d <- x$data
  for (type in c("lf", "wf")) {
    fields <- list(switch = family, n_f = c(2L, 2L), fishery_f = 1:2,
      year_fi = list(1:2, 1:2), n_fi = list(c(120,120), c(120,120)),
      minbin = c(1L,1L), maxbin = c(15L,15L), obs_flat = rep(1:15,4),
      obs_ints = rep(1:15,4), obs_prop = rep((1:15)/120,4))
    names(fields) <- paste0(type, "_", names(fields))
    d[names(fields)] <- fields
    d[[paste0("n_",type)]] <- 4L
  }
  d$n_wt <- 15L
  d$wf_rebin_matrix <- diag(15L)
  x <- opal_build(opal_update(x, data = d))
  o <- opal_rtmb(x)
  opal_attach_fit(x, list(par = o$par, objective = o$fn(o$par), convergence = 0L), check = FALSE)
}

test_that("OSA covers both composition streams and all families with exclusions", {
  for (family in 1:3) {
    x <- osa_fixture(family)
    set.seed(100); rng <- .Random.seed
    x <- opal_osa(x, seed = 3)
    expect_identical(.Random.seed, rng)
    result <- opal_derived(x, "osa")
    expect_setequal(result$summary$dataset, c("CPUE 1", "LF 1", "LF 2", "WF 1", "WF 2"))
    expect_equal(result$summary$n, c(2L, 28L, 28L, 28L, 28L))
    expect_true(all(result$summary$failures == 0))
    expect_identical(opal_derived(opal_osa(x, seed=3), "osa"), result)
    expect_s3_class(plot_osa_sdnr(x), "ggplot")
    expect_equal(plot_osa_sdnr(x)$layers[[1]]$data$xintercept, 1)
    expect_s3_class(plot_composition(x, "wf"), "ggplot")
    expect_s3_class(plot_osa_residuals(x, type="qq"), "ggplot")
    d <- x$data; d$removal_switch_f[2] <- 1L
    x <- opal_build(opal_update(x, data=d)); o <- opal_rtmb(x)
    x <- opal_attach_fit(x,list(par=o$par,objective=o$fn(o$par),convergence=0L),check=FALSE)
    result <- opal_derived(opal_osa(x), "osa")
    expect_equal(result$summary$n[result$summary$dataset %in% c("LF 2", "WF 2")], c(0L,0L))
  }
})

test_that("SDNR intervals, cache integrity, and CPUE residuals are explicit", {
  x <- opal_osa(opal_fit(small_opal_object()))
  result <- opal_derived(x, "osa")
  r <- opal_report(x)
  expect_equal(result$residuals$residual, (log(x$data$cpue_data$value)-log(r$cpue_pred))/r$cpue_sigma)
  z <- c(-1, 0, 1)
  s <- opal:::.opal_sdnr(z, .95)
  expect_equal(s$sdnr, 1)
  expect_equal(s$lower, sqrt(2/qchisq(.975,2)))
  expect_equal(s$upper, sqrt(2/qchisq(.025,2)))
  x$derived$osa$result$summary$sdnr <- 1
  expect_error(opal_derived(x,"osa"), "modified")
})

test_that("Dirichlet-multinomial OSA agrees with conditional beta-binomial probabilities", {
  observation <- c(2L,3L,5L)
  model <- function(p) {
    observed <- RTMB::OBS(observation)
    opal:::.opal_ddirmult(observed,10,exp(p$log_alpha),log=TRUE) * -1
  }
  o <- RTMB::MakeADFun(model,list(log_alpha=log(c(2,3,5))),silent=TRUE)
  expect_equal(o$fn(o$par),-RTMBdist::ddirmult(observation,10,c(2,3,5),log=TRUE))
  res <- RTMB::oneStepPredict(o,observation.name="observation",method="oneStepGeneric",
    discrete=TRUE,range=c(0,Inf),subset=1:2,seed=44,trace=FALSE)$residual
  p <- pnorm(res[1:2])
  cdf <- function(q,n,a,b) if(q<0) 0 else sum(choose(n,0:q)*beta((0:q)+a,n-(0:q)+b)/beta(a,b))
  expect_true(p[1] >= cdf(1,10,2,8) - 1e-6 && p[1] <= cdf(2,10,2,8) + 1e-6)
  expect_true(p[2] >= cdf(2,8,3,5) - 1e-6 && p[2] <= cdf(3,8,3,5) + 1e-6)
})

test_that("marginal OSA returns residuals for every composition family", {
  for (family in 1:3) {
    x <- osa_fixture(family)
    map <- x$map; map$rdev_y <- NULL
    x <- opal_fit(opal_update(x,map=map,random="rdev_y"),check=FALSE)
    x <- opal_osa(x)
    result <- opal_derived(x,"osa")
    expect_true(all(result$summary$failures==0))
    expect_equal(result$summary$n,c(2L,28L,28L,28L,28L))
    expect_true(all(result$summary$source=="RTMB marginal OSA"))
  }
})

test_that("OSA honours the model's omitted likelihood defaults", {
  inputs <- opal_example_inputs()
  inputs$data$cpue_switch <- NULL
  x <- opal_fit(opal_obj(inputs$data,inputs$parameters,inputs$map),check=FALSE)
  expect_false(any(opal_derived(opal_osa(x),"osa")$summary$family=="lognormal"))
  inputs$data$n_lf <- NULL
  inputs$data$n_wf <- NULL
  x <- opal_build(opal_obj(inputs$data,inputs$parameters,inputs$map))
  o <- opal_rtmb(x)
  x <- opal_attach_fit(x,list(par=o$par,objective=o$fn(o$par),convergence=0L),check=FALSE)
  expect_error(opal_osa(x),"No active observation")
})
