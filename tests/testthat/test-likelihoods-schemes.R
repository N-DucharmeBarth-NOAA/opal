scheme_like_fixture <- function(prefix, type, scheme = c(1L, 2L, 3L, 2L)) {
  fx <- comp_fixture()
  coarse <- matrix(0, 3, fx$n_len)
  coarse[cbind(rep(1:3, each = 4), seq_len(fx$n_len))] <- 1
  uneven <- matrix(0, 4, fx$n_len)
  uneven[cbind(c(1, 1, 2, 2, 2, 3, 3, 4, 4, 4, 4, 4),
                seq_len(fx$n_len))] <- 1
  uneven[, 3] <- c(0.25, 0.75, 0, 0)
  matrices <- list(diag(fx$n_len), coarse, uneven)
  nbins <- vapply(matrices[scheme], nrow, integer(1))
  counts <- lapply(nbins, function(n) as.integer(seq_len(n) + 2L))
  sizes <- vapply(counts, sum, integer(1))
  args <- list(catch_pred_fya = fx$catch, pla = fx$pla,
               removal_switch_f = c(0L, 0L))
  put <- function(key, value) args[[paste0(prefix, "_", key)]] <<- value
  put("obs_flat", unlist(counts))
  put("obs_ints", unlist(counts))
  put("obs_prop", unlist(lapply(counts, function(x) x / sum(x))))
  put("n_f", c(3L, 1L))
  put("fishery_f", c(1L, 2L))
  put("year_fi", list(c(1L, 2L, 4L), 3L))
  put("n_fi", list(sizes[1:3], sizes[4]))
  put("n_int_fi", list(sizes[1:3], sizes[4]))
  put("minbin", c(1L, 1L))
  put("maxbin", rep(fx$n_len, 2))
  put("switch", type)
  put("rebin_s", matrices)
  put("scheme", scheme)
  put("nbins_i", nbins)
  args[[paste0("n_", prefix)]] <- 4L
  args[[paste0("log_", prefix, "_tau")]] <- log(c(15, 25))
  if (prefix == "lf") {
    args$n_len <- fx$n_len
  } else {
    args$n_wt <- fx$n_len
    args$wf_rebin_matrix <- diag(fx$n_len)
  }
  list(args = args, counts = counts, sizes = sizes,
       fishery = c(1L, 1L, 1L, 2L), year = c(1L, 2L, 4L, 3L))
}

scheme_like_function <- function(prefix) {
  if (prefix == "lf") get_length_like else get_weight_like
}

test_that("exact, coarse, and uneven schemes preserve all three likelihood formulas", {
  for (prefix in c("lf", "wf")) for (type in 1:3) {
    fixture <- scheme_like_fixture(prefix, type)
    args <- fixture$args
    matrices <- args[[paste0(prefix, "_rebin_s")]]
    scheme <- args[[paste0(prefix, "_scheme")]]
    expected_pred <- lapply(seq_len(4), function(i) {
      length_mass <- args$pla %*%
        args$catch_pred_fya[fixture$fishery[i], fixture$year[i], ]
      mass <- c(matrices[[scheme[i]]] %*% length_mass) + 1e-8
      mass / sum(mass)
    })
    expected_lp <- vapply(seq_len(4), function(i) {
      x <- fixture$counts[[i]]
      pred <- expected_pred[[i]]
      tau <- exp(args[[paste0("log_", prefix, "_tau")]][fixture$fishery[i]])
      switch(type,
             ref_multinom_nll(x, pred),
             ref_dirichlet_nll(x / sum(x), pred * fixture$sizes[i] * tau),
             ref_dirmult_nll(x, pred * tau, size = fixture$sizes[i]))
    }, numeric(1))
    fun <- scheme_like_function(prefix)
    expect_equal(do.call(fun, args), expected_lp, tolerance = 1e-10)
    obj <- RTMB::MakeADFun(function(p) {
      a <- args
      a[[paste0("log_", prefix, "_tau")]] <- p$log_tau
      sum(do.call(fun, a))
    }, list(log_tau = args[[paste0("log_", prefix, "_tau")]]), silent = TRUE)
    rpt <- obj$report()
    expect_equal(rpt[[paste0(prefix, "_pred_row")]], expected_pred)
    expect_type(rpt[[paste0(prefix, "_pred")]][[1]], "list")
    expect_equal(dim(rpt[[paste0(prefix, "_pred")]][[2]]), c(1L, 3L))
    expect_equal(vapply(rpt[[paste0(prefix, "_pred_row")]], sum, numeric(1)),
                 rep(1, 4))
    expect_equal(vapply(fixture$counts, sum, integer(1)), fixture$sizes)
  }
})

test_that("scheme gradients include catch, concentration, and estimated growth", {
  for (prefix in c("lf", "wf")) for (type in 1:3) {
    args <- scheme_like_fixture(prefix, type)$args
    fun <- scheme_like_function(prefix)
    obj <- RTMB::MakeADFun(function(p) {
      "[<-" <- RTMB::ADoverload("[<-")
      a <- args
      mu <- get_growth(6L, 1, 6, exp(p$log_L1), exp(p$log_L2), p$log_k)
      sd <- get_sd_at_age(mu, exp(p$log_L1), exp(p$log_L2), log(0.2), log(0.15))
      a$pla <- get_pla(seq(0, 11), seq(1, 12), mu, sd)
      a$catch_pred_fya[, , 1] <- a$catch_pred_fya[, , 1] * exp(p$log_catch)
      a[[paste0("log_", prefix, "_tau")]] <- p$log_tau
      sum(do.call(fun, a))
    }, list(log_L1 = log(2), log_L2 = log(10), log_k = log(0.25),
            log_catch = 0, log_tau = log(c(15, 25))), silent = TRUE)
    par <- obj$par
    ad <- as.numeric(obj$gr(par))
    eps <- 1e-5
    fd <- vapply(seq_along(par), function(i) {
      up <- down <- par
      up[i] <- up[i] + eps
      down[i] <- down[i] - eps
      (obj$fn(up) - obj$fn(down)) / (2 * eps)
    }, numeric(1))
    expect_true(is.finite(obj$fn(par)))
    expect_true(all(is.finite(ad)))
    expect_equal(ad, fd, tolerance = 1e-5)
    expect_gt(max(abs(ad[1:3])), 1e-3)
  }
})

test_that("scheme DM rows validate integer totals without mixing row offsets", {
  for (prefix in c("lf", "wf")) {
    args <- scheme_like_fixture(prefix, 3L)$args
    key <- paste0(prefix, "_obs_ints")
    args[[key]][13] <- args[[key]][13] + 1L
    expect_error(do.call(scheme_like_function(prefix), args),
                 "counts must sum to the integer sample size")
  }
})

test_that("single exact schemes retain legacy predictions and fishery matrices", {
  for (prefix in c("lf", "wf")) for (type in 1:3) {
    args <- scheme_like_fixture(prefix, type, rep(1L, 4))$args
    fun <- scheme_like_function(prefix)
    scheme_lp <- do.call(fun, args)
    legacy <- args
    legacy[paste0(prefix, c("_rebin_s", "_scheme", "_nbins_i"))] <- NULL
    expect_equal(scheme_lp, do.call(fun, legacy), tolerance = 1e-12)
    obj <- RTMB::MakeADFun(function(p) {
      a <- args
      a[[paste0("log_", prefix, "_tau")]] <- p$tau
      sum(do.call(fun, a))
    }, list(tau = args[[paste0("log_", prefix, "_tau")]]), silent = TRUE)
    rpt <- obj$report()
    expect_equal(dim(rpt[[paste0(prefix, "_pred")]][[1]]), c(3L, 12L))
    expect_equal(rpt[[paste0(prefix, "_pred")]][[1]][1, ],
                 rpt[[paste0(prefix, "_pred_row")]][[1]])
  }
})

test_that("opal_model forwards scheme inputs and supports data without them", {
  for (type in 1:3) {
    legacy <- synth_full_data(wf_switch = type, lf_switch = type)
    legacy$wf_n_f <- 2L
    legacy$wf_obs_ints <- rep(1L, 2L * legacy$n_wt)
    d <- legacy
    d$lf_rebin_s <- list(diag(d$n_len))
    d$lf_scheme <- rep(1L, d$n_lf)
    d$lf_nbins_i <- rep(d$n_len, d$n_lf)
    d$wf_rebin_s <- list(d$wf_rebin_matrix)
    d$wf_scheme <- rep(1L, d$n_wf)
    d$wf_nbins_i <- rep(d$n_wt, d$n_wf)
    obj <- synth_obj(d)
    old <- synth_obj(legacy)
    expect_equal(obj$fn(), old$fn(), tolerance = 1e-10)
    expect_equal(obj$gr(), old$gr(), tolerance = 1e-10)
    for (prefix in c("lf", "wf")) {
      rpt <- obj$report()
      expect_equal(rpt[[paste0(prefix, "_pred")]],
                   old$report()[[paste0(prefix, "_pred")]], tolerance = 1e-10)
      expect_length(rpt[[paste0(prefix, "_pred_row")]], d[[paste0("n_", prefix)]])
    }
  }
})

test_that("opal_model routes mixed schemes to prediction-only aggregation", {
  for (type in 1:3) {
    d <- synth_full_data(wf_switch = type, lf_switch = type)
    d$wf_n_f <- 2L
    coarse <- matrix(0, 3L, d$n_len)
    coarse[cbind(rep(1:3, each = 5L), seq_len(d$n_len))] <- 1
    for (prefix in c("lf", "wf")) {
      d[[paste0(prefix, "_rebin_s")]] <- list(coarse, diag(d$n_len))
      d[[paste0(prefix, "_scheme")]] <- c(1L, 2L)
      d[[paste0(prefix, "_nbins_i")]] <- c(3L, d$n_len)
      d[[paste0(prefix, "_obs_flat")]] <- c(rep(5, 3), rep(1, d$n_len))
      d[[paste0(prefix, "_obs_ints")]] <- as.integer(d[[paste0(prefix, "_obs_flat")]])
      d[[paste0(prefix, "_obs_prop")]] <- d[[paste0(prefix, "_obs_flat")]] / 15
    }
    parameters <- synth_parameters(d)
    obj <- synth_obj(d, parameters)
    rpt <- obj$report()
    expect_true(is.finite(obj$fn()))
    for (prefix in c("lf", "wf")) {
      pred <- rpt[[paste0(prefix, "_pred_row")]]
      expect_equal(lengths(pred), c(3L, d$n_len))
      tau <- exp(parameters[[paste0("log_", prefix, "_tau")]][1])
      counts <- list(rep(5, 3), rep(1, d$n_len))
      expected <- vapply(1:2, function(i) {
        switch(type,
               ref_multinom_nll(counts[[i]], pred[[i]]),
               ref_dirichlet_nll(counts[[i]] / 15, pred[[i]] * 15 * tau),
               ref_dirmult_nll(counts[[i]], pred[[i]] * tau, size = 15))
      }, numeric(1))
      expect_equal(rpt[[paste0("lp_", prefix)]], expected, tolerance = 1e-10)
    }
  }
})

test_that("fleet lambdas zero every row across mixed LF and WF schemes", {
  d <- synth_full_data(wf_switch = 1L, lf_switch = 1L)
  d$lambdas <- list(lf = c(0, 1), wf = c(0, 1))
  d$n_lf <- d$n_wf <- 3L
  d$lf_n_f <- d$wf_n_f <- c(2L, 1L)
  d$lf_fishery_f <- d$wf_fishery_f <- c(1L, 2L)
  d$lf_year_fi <- d$wf_year_fi <- list(c(1L, 2L), 1L)
  d$lf_n_fi <- d$wf_n_fi <- list(c(15, 15), 15)
  d$lf_n_int_fi <- d$wf_n_int_fi <- list(c(15L, 15L), 15L)
  coarse <- matrix(0, 3L, d$n_len)
  coarse[cbind(rep(1:3, each = 5L), seq_len(d$n_len))] <- 1
  for (prefix in c("lf", "wf")) {
    d[[paste0(prefix, "_rebin_s")]] <- list(coarse, diag(d$n_len))
    d[[paste0(prefix, "_scheme")]] <- c(1L, 2L, 1L)
    d[[paste0(prefix, "_nbins_i")]] <- c(3L, d$n_len, 3L)
    counts <- c(rep(5, 3), rep(1, d$n_len), rep(5, 3))
    d[[paste0(prefix, "_obs_flat")]] <- counts
    d[[paste0(prefix, "_obs_ints")]] <- as.integer(counts)
    d[[paste0(prefix, "_obs_prop")]] <- c(rep(1 / 3, 3),
      rep(1 / d$n_len, d$n_len), rep(1 / 3, 3))
  }
  obj <- synth_obj(d)
  report <- obj$report()
  for (prefix in c("lf", "wf")) {
    expect_equal(report[[paste0("lp_", prefix, "_w")]],
      report[[paste0("lp_", prefix)]] * c(0, 0, 1))
  }
})

test_that("prepared mixed observations remain aligned across fisheries", {
  make_rows <- function(f, y, edges, counts) {
    data.frame(fishery = f, year = y, month = 1L, ts = y,
               bin_lower = head(edges, -1L), bin_upper = tail(edges, -1L),
               count = counts)
  }
  lf <- rbind(make_rows(2L, 1L, c(20, 30, 50), c(12, 18)),
              make_rows(1L, 1L, c(20, 24, 36, 50), c(6, 8, 16)),
              make_rows(2L, 2L, c(20, 24, 36, 50), c(5, 10, 15)))
  wf <- rbind(make_rows(2L, 1L, c(0.08, 0.3, 1.25), c(12, 18)),
              make_rows(1L, 1L, c(0.08, 0.2, 0.5, 1.25), c(6, 8, 16)),
              make_rows(2L, 2L, c(0.08, 0.2, 0.5, 1.25), c(5, 10, 15)))
  for (type in 1:3) {
    d <- synth_data()
    d$wt_bin_start <- 0
    d$wt_bin_width <- 0.1
    d$n_wt <- 15L
    d <- prep_lf_data(d, lf, lf_switch = type)
    d <- prep_wf_data(d, wf, wf_switch = type)
    obj <- synth_obj(d)
    rpt <- obj$report()
    expect_true(is.finite(obj$fn()))
    expect_true(all(is.finite(obj$gr())))
    for (prefix in c("lf", "wf")) {
      expect_equal(unlist(d[[paste0(prefix, "_row_fi")]], use.names = FALSE),
                   c(2L, 1L, 3L))
      pred <- rpt[[paste0(prefix, "_pred_row")]]
      expect_equal(lengths(pred), c(3L, 2L, 3L))
      expect_equal(lengths(pred), d[[paste0(prefix, "_nbins_i")]])
      if (type == 3L) {
        ends <- cumsum(lengths(pred))
        starts <- c(1L, head(ends, -1L) + 1L)
        counts <- d[[paste0(prefix, "_obs_ints")]]
        totals <- vapply(1:3, function(i) sum(counts[starts[i]:ends[i]]), numeric(1))
        expect_equal(totals,
                     unlist(d[[paste0(prefix, "_n_int_fi")]], use.names = FALSE))
      }
    }
  }
})

test_that("prepared native min/max schemes retain the legacy two-stage mapping", {
  make_wide <- function(lower) {
    x <- data.frame(fishery = c(2L, 1L, 2L), year = c(1L, 1L, 2L),
                    month = 1L, ts = c(1L, 1L, 2L))
    for (i in seq_along(lower)) x[[as.character(lower[i])]] <- c(i, 16L - i, i + 2L)
    x
  }
  for (type in 1:3) {
    d <- synth_data()
    d$wt_bin_start <- 0
    d$wt_bin_width <- 0.1
    d$n_wt <- 15L
    d <- prep_lf_data(d, make_wide(d$len_lower), lf_switch = type,
                      lf_minbin = c(3L, 2L), lf_maxbin = c(12L, 10L))
    expect_warning(
      d <- prep_wf_data(d, make_wide(seq(0, by = 0.1, length.out = d$n_wt)),
                        wf_switch = type, wf_minbin = c(2L, 3L),
                        wf_maxbin = c(11L, 10L)),
      "narrower"
    )
    legacy <- d
    for (prefix in c("lf", "wf")) {
      legacy[paste0(prefix, c("_rebin_s", "_scheme", "_nbins_i"))] <- NULL
    }
    obj <- synth_obj(d)
    old <- synth_obj(legacy)
    expect_equal(obj$fn(), old$fn(), tolerance = 1e-10)
    expect_equal(obj$gr(), old$gr(), tolerance = 1e-10)
    for (prefix in c("lf", "wf")) {
      expect_equal(obj$report()[[paste0(prefix, "_pred")]],
                   old$report()[[paste0(prefix, "_pred")]], tolerance = 1e-10)
      expect_equal(obj$report()[[paste0(prefix, "_pred_row")]],
                   old$report()[[paste0(prefix, "_pred_row")]], tolerance = 1e-10)
    }
  }
})
