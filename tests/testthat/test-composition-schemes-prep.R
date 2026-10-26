scheme_prep_data <- function() {
  list(len_bin_start = 0, len_bin_width = 2, n_len = 4L, n_fishery = 2L,
       wt_bin_start = 0, wt_bin_width = 2, n_wt = 4L, lw_a = 1, lw_b = 1)
}

scheme_wide <- function(fishery = 1L, lower = c(0, 4), counts = c(3, 7),
                        ts = 1L) {
  d <- data.frame(fishery = fishery, year = 2000L, month = 1L, ts = ts)
  for (i in seq_along(lower)) d[[as.character(lower[i])]] <- counts[i]
  d
}

test_that("mixed schemes preserve counts and align with fishery-grouped flat order", {
  x <- list(scheme_wide(2L), scheme_wide(1L, c(0, 2, 6), c(1, 2, 3)),
            scheme_wide(2L, ts = 2L))
  d <- prep_lf_data(scheme_prep_data(), x, lf_upper_edge = c(8, 8, 8),
                    lf_addtocomp = 0)
  expect_equal(d$lf_row_fi, list(`1` = 2L, `2` = c(1L, 3L)))
  expect_equal(d$lf_scheme, c(1L, 2L, 2L))
  expect_equal(d$lf_nbins_i, c(3L, 2L, 2L))
  expect_equal(d$lf_obs_flat, c(1, 2, 3, 3, 7, 3, 7))
  expect_equal(d$lf_scheme_edges, list(c(0, 2, 6, 8), c(0, 4, 8)))
  expect_equal(d$lf_rebin_s[[1]], rbind(c(1, 0, 0, 0),
                                      c(0, 1, 1, 0), c(0, 0, 0, 1)))
  expect_equal(d$lf_obs_ints, as.integer(d$lf_obs_flat))
})

test_that("long input uses source grouping and retains independently observed bins", {
  x <- data.frame(fishery = 1L, year = 2000L, month = 1L, ts = 1L,
                  source = rep(c("a", "b"), each = 2),
                  bin_lower = rep(c(0, 4), 2), bin_upper = rep(c(4, 8), 2),
                  count = c(3, 7, 2, 8))
  d <- prep_lf_data(scheme_prep_data(), x, lf_addtocomp = 0)
  expect_equal(d$n_lf, 2L)
  expect_equal(d$lf_scheme, c(1L, 1L))
  expect_equal(d$lf_obs_flat, c(3, 7, 2, 8))
  x$bin_upper[1] <- 3
  expect_error(prep_lf_data(scheme_prep_data(), x), "contiguous")
})

test_that("caps adjustments and the two legacy smoothing steps are preserved", {
  x <- scheme_wide(counts = c(0, 100))
  d <- prep_lf_data(scheme_prep_data(), x, lf_upper_edge = 8,
                    lf_cap = 20, lf_var_adjust = c(2, 1), lf_addtocomp = 0.01)
  p <- (c(0, 1) + 0.01) / 1.02
  expect_equal(d$lf_n, 10)
  expect_equal(d$lf_obs_flat, p * 10)
  expect_equal(d$lf_obs_prop, (p + 0.01) / 1.02)
  expect_equal(sum(d$lf_obs_ints), 10L)
})

test_that("native min/max mappings reproduce legacy weight row pooling exactly", {
  d <- scheme_prep_data()
  d$lw_a <- 0.7
  d$lw_b <- 1.2
  x <- scheme_wide(lower = c(0, 2, 4, 6), counts = c(1, 2, 3, 4))
  z <- prep_wf_data(d, x, wf_minbin = c(2L, 1L), wf_maxbin = c(3L, 4L),
                    wf_addtocomp = 0)
  m <- z$wf_rebin_matrix
  expect_equal(z$wf_rebin_s[[1]], rbind(colSums(m[1:2, ]),
                                       colSums(m[3:4, ])), tolerance = 0)
  expect_equal(z$wf_scheme_edges, list(c(2, 4, 6)))
  expect_equal(z$wf_obs_flat, c(3, 7))
  expect_equal(z$wf_nbins_i, 2L)
})

test_that("native schemes retain legacy smoothing and fishery ordering", {
  x <- rbind(scheme_wide(2L, c(0, 2, 4, 6), c(0, 2, 5, 3)),
             scheme_wide(1L, c(0, 2, 4, 6), c(1, 2, 3, 4)))
  d <- scheme_prep_data()
  for (kind in c("lf", "wf")) {
    prep <- get(paste0("prep_", kind, "_data"))
    current <- prep(d, x)
    p <- rbind(c(1, 2, 3, 4), c(0, 2, 5, 3)) / 10
    p <- (p + 1e-8) / (1 + 4e-8)
    expect_equal(current[[paste0(kind, "_obs_flat")]], as.numeric(t(p * 10)))
    expect_equal(current[[paste0(kind, "_obs_prop")]],
                 as.numeric(t((p + 1e-8) / (1 + 4e-8))))
    expect_equal(current[[paste0(kind, "_row_fi")]], list(`1` = 2L, `2` = 1L))
  }
})

test_that("accumulation changes prediction matrix edges but not registry edges", {
  x <- scheme_wide(lower = c(2, 4), counts = c(2, 8))
  a <- prep_lf_data(scheme_prep_data(), x, lf_upper_edge = 6, lf_addtocomp = 0)
  t <- prep_lf_data(scheme_prep_data(), x, lf_upper_edge = 6,
                    lf_tails = "truncate", lf_addtocomp = 0)
  expect_equal(a$lf_scheme_edges, list(c(2, 4, 6)))
  expect_equal(colSums(a$lf_rebin_s[[1]]), rep(1, 4))
  expect_equal(colSums(t$lf_rebin_s[[1]]), c(0, 1, 1, 0))
  expect_identical(a$lf_obs_flat, t$lf_obs_flat)
})

test_that("misalignment diagnoses the affected observations and supports opt-in", {
  x <- scheme_wide(lower = c(0, 3), counts = c(1, 9))
  expect_error(prep_lf_data(scheme_prep_data(), x, lf_upper_edge = 8),
               "scheme 1, fisheries 1, rows 1.*edges 3.*split share.*grid width.*par_sel")
  expect_warning(d <- prep_lf_data(scheme_prep_data(), x, lf_upper_edge = 8,
                                   lf_allow_misaligned = TRUE), "misaligned edges")
  expect_equal(d$lf_rebin_s[[1]][, 2L], c(0.5, 0.5))
  expect_equal(d$lf_scheme_edges, list(c(0, 3, 8)))
})

test_that("near edges snap silently and schemes deduplicate within tolerance", {
  x <- list(scheme_wide(), scheme_wide(lower = c(0, 4 + 1e-9)))
  expect_silent(d <- prep_lf_data(scheme_prep_data(), x, lf_upper_edge = 8))
  expect_equal(d$lf_scheme, c(1L, 1L))
  expect_equal(length(d$lf_rebin_s), 1L)
  expect_equal(d$lf_rebin_s[[1]], rbind(c(1, 1, 0, 0), c(0, 0, 1, 1)))
})

test_that("narrow bins and entirely outside weight schemes warn", {
  x <- scheme_wide(lower = c(0, 1), counts = c(2, 8))
  expect_warning(expect_warning(prep_lf_data(scheme_prep_data(), x,
    lf_upper_edge = 8, lf_allow_misaligned = TRUE), "narrower"), "misaligned")
  w <- scheme_wide(lower = c(20, 22))
  expect_warning(prep_wf_data(scheme_prep_data(), w, wf_upper_edge = 24),
                 "entirely outside")
})

test_that("invalid counts and edges are rejected", {
  for (counts in list(c(-1, 2), c(NA_real_, 1), c(Inf, 1))) {
    expect_error(prep_lf_data(scheme_prep_data(), scheme_wide(counts = counts),
                             lf_upper_edge = 8), "finite and nonnegative")
  }
  expect_error(prep_lf_data(scheme_prep_data(), scheme_wide(), lf_upper_edge = 3),
               "increasing")
  expect_error(prep_lf_data(scheme_prep_data(), scheme_wide()), "upper_edge")
})

test_that("singleton zero-row filtered and empty inputs have valid layouts", {
  x <- scheme_wide(lower = 0, counts = 10)
  for (kind in c("lf", "wf")) {
    prep <- get(paste0("prep_", kind, "_data"))
    args <- setNames(list(8, 0), paste0(kind, c("_upper_edge", "_addtocomp")))
    d <- do.call(prep, c(list(scheme_prep_data(), x), args))
    expect_equal(d[[paste0(kind, "_obs_flat")]], 10)
    expect_equal(dim(d[[paste0(kind, "_rebin_s")]][[1]]), c(1L, 4L))
    x0 <- x
    x0[["0"]] <- 0
    for (empty in list(x[FALSE, ], x0, list())) {
      z <- do.call(prep, c(list(scheme_prep_data(), empty), args))
      expect_identical(z[[paste0("n_", kind)]], 0L)
      expect_identical(z[[paste0(kind, "_obs_flat")]], numeric())
      expect_identical(z[[paste0(kind, "_obs_ints")]], integer())
      expect_identical(z[[paste0(kind, "_scheme")]], integer())
    }
  }
})
