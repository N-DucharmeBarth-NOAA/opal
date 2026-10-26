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
  expect_identical(d$lf_native_i, c(FALSE, FALSE, FALSE))
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
  expect_warning(z <- prep_wf_data(d, x, wf_minbin = c(2L, 1L),
    wf_maxbin = c(3L, 4L), wf_addtocomp = 0), "narrower")
  m <- z$wf_rebin_matrix
  expect_equal(z$wf_rebin_s[[1]], rbind(colSums(m[1:2, ]),
                                       colSums(m[3:4, ])), tolerance = 0)
  expect_equal(z$wf_scheme_edges, list(c(2, 4, 6)))
  expect_equal(z$wf_obs_flat, c(3, 7))
  expect_equal(z$wf_nbins_i, 2L)
  expect_identical(z$wf_native_i, TRUE)
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
  expect_equal(prep_lf_data(scheme_prep_data(), scheme_wide())$lf_scheme_edges,
               list(c(0, 4, 8)))
  expect_error(prep_lf_data(scheme_prep_data(), list(scheme_wide(), scheme_wide()),
                            lf_upper_edge = c(8, 8, 8)), "length one")
})

test_that("long observation grouping ignores changing calendar metadata", {
  x <- data.frame(fishery = 1L, year = c(2000L, 2001L), month = c(1L, 2L),
                  ts = 1L, bin_lower = c(0, 4), bin_upper = c(4, 8), count = c(2, 8))
  d <- prep_lf_data(scheme_prep_data(), x, lf_addtocomp = 0)
  expect_equal(d$n_lf, 1L)
  expect_equal(d$lf_obs_flat, c(2, 8))
})

test_that("nonmodel LF edges outside the source range require explicit opt-in", {
  x <- scheme_wide(lower = c(-2, 4))
  expect_error(prep_lf_data(scheme_prep_data(), x),
               "misaligned edges -2.*Check bin labeling/units.*extend its range")
  expect_warning(d <- prep_lf_data(scheme_prep_data(), x,
    lf_allow_misaligned = TRUE), "misaligned")
  expect_equal(d$lf_scheme_edges, list(c(-2, 4, 8)))
  expect_equal(colSums(d$lf_rebin_s[[1]]), rep(1, 4))
})

test_that("outside-only LF misalignment has a finite zero actual split diagnostic", {
  x <- scheme_wide(lower = c(-1.2, 4))
  expect_error(prep_lf_data(scheme_prep_data(), x, lf_upper_edge = 9.2),
               "misaligned edges -1\\.2,9\\.2; maximum split share 0;")
  expect_warning(d <- prep_lf_data(scheme_prep_data(), x, lf_upper_edge = 9.2,
    lf_allow_misaligned = TRUE), "maximum split share 0;")
  expect_equal(d$lf_scheme_edges, list(c(-1.2, 4, 9.2)))
  expect_true(all(is.finite(d$lf_rebin_s[[1]])))
})

test_that("custom weight bins compare widths with every overlapping source bin", {
  d <- scheme_prep_data()
  d$lw_a <- 2
  x <- scheme_wide(lower = c(0, 2), counts = c(2, 8))
  expect_warning(w <- prep_wf_data(d, x), "narrower than overlapping source")
  expect_equal(w$wf_scheme_edges, list(c(0, 2, 4)))
})

test_that("shifted full-size weight grids default to their observed upper edge", {
  x <- scheme_wide(lower = c(1, 3, 5, 7), counts = c(1, 2, 3, 4))
  d <- prep_wf_data(scheme_prep_data(), x)
  expect_equal(d$wf_scheme_edges, list(c(1, 3, 5, 7, 9)))
  expect_equal(colSums(d$wf_rebin_s[[1]]), rep(1, 4))
})

test_that("equal displayed weight edges do not merge different native tail mappings", {
  d <- scheme_prep_data()
  d$lw_a <- 2
  x <- list(scheme_wide(lower = c(0, 2, 4, 6), counts = c(1, 2, 3, 4)),
            scheme_wide(lower = c(2, 4), counts = c(3, 7)))
  expect_warning(expect_warning(z <- prep_wf_data(d, x, wf_upper_edge = c(8, 6),
    wf_minbin = c(2L, 1L), wf_maxbin = c(3L, 4L)), "narrower"), "narrower")
  expect_equal(z$wf_scheme_edges, list(c(2, 4, 6), c(2, 4, 6)))
  expect_identical(z$wf_native_i, c(TRUE, FALSE))
  expect_equal(z$wf_scheme, c(1L, 2L))
  expect_equal(colSums(z$wf_rebin_s[[1]]), c(1, 1, 0, 0))
  expect_equal(colSums(z$wf_rebin_s[[2]]), rep(1, 4))
  expect_warning(expect_warning(reverse <- prep_wf_data(d, rev(x),
    wf_upper_edge = c(6, 8), wf_minbin = c(2L, 1L), wf_maxbin = c(3L, 4L)),
    "narrower"), "narrower")
  expect_identical(reverse$wf_native_i, c(FALSE, TRUE))
  expect_equal(reverse$wf_rebin_s, rev(z$wf_rebin_s))
})

test_that("long list and native inputs agree on a native observation grid", {
  wide <- scheme_wide(lower = c(0, 2, 4, 6), counts = c(1, 2, 3, 4))
  long <- data.frame(fishery = 1L, year = 2000L, month = 1L, ts = 1L,
                     bin_lower = c(0, 2, 4, 6), bin_upper = c(2, 4, 6, 8),
                     count = c(1, 2, 3, 4))
  for (kind in c("lf", "wf")) {
    prep <- get(paste0("prep_", kind, "_data"))
    args <- setNames(list(c(2L, 1L), c(3L, 4L)),
                     paste0(kind, c("_minbin", "_maxbin")))
    reference <- do.call(prep, c(list(scheme_prep_data(), wide), args))
    for (input in list(list(wide), long)) {
      d <- do.call(prep, c(list(scheme_prep_data(), input), args))
      for (field in c("obs_flat", "obs_ints", "obs_prop", "scheme_edges", "rebin_s",
                       "scheme", "nbins_i", "native_i"))
        expect_identical(d[[paste0(kind, "_", field)]],
                         reference[[paste0(kind, "_", field)]])
    }
  }
})

test_that("filtering interleaved rows retains grouped sizes and DM allocations", {
  x <- list(scheme_wide(2L, ts = 10L),
            scheme_wide(1L, c(0, 2, 6), c(1, 2, 3), ts = 11L),
            scheme_wide(2L, ts = 12L),
            scheme_wide(1L, counts = c(0, 0), ts = 13L))
  d <- prep_lf_data(scheme_prep_data(), x, lf_upper_edge = 8,
                    lf_var_adjust = c(2, 4), lf_addtocomp = 0, lf_switch = 3L)
  expect_equal(d$lf_n_int, c(3L, 3L, 3L))
  expect_equal(d$lf_obs_ints, c(1L, 1L, 1L, 1L, 2L, 1L, 2L))
  expect_equal(d$lf_nbins_i, c(3L, 2L, 2L))
  f <- prep_lf_data(scheme_prep_data(), x, lf_keep_fisheries = 2L,
                    lf_addtocomp = 0)
  expect_equal(f$lf_year_fi, list(`2` = c(10L, 12L)))
  expect_equal(f$lf_row_fi, list(`2` = c(1L, 2L)))
  expect_equal(f$lf_obs_flat, c(3, 7, 3, 7))
})

test_that("custom weight matrices overlap converted length edges analytically", {
  d <- scheme_prep_data()
  d$lw_b <- 2
  x <- scheme_wide(lower = c(0, 8), counts = c(4, 6))
  expect_warning(z <- prep_wf_data(d, x, wf_upper_edge = 64), "narrower")
  expect_equal(z$wf_rebin_s[[1]],
               rbind(c(1, 1 / 3, 0, 0), c(0, 2 / 3, 1, 1)))
})

test_that("finer grid suggestions use the GCD of all decimal observed offsets", {
  d <- scheme_prep_data()
  d$len_bin_width <- 1
  x <- scheme_wide(lower = c(0, 1.2), counts = c(4, 6))
  expect_error(prep_lf_data(d, x),
               "suggested largest grid width dividing offsets: 0\\.4")
  expect_equal(opal:::.opal_composition_grid_width(c(0, 1.2, 4), 0:4, 1e-8),
               0.4)
  expect_equal(opal:::.opal_composition_grid_width(c(0, 2, 4), 0:4, 1e-8), 1)
  expect_equal(opal:::.opal_composition_grid_width(c(0, 1.2 + 1e-10, 4), 0:4, 1e-8),
               0.4, tolerance = 1e-8)
})

test_that("integer-valued double native bin indices are coerced after validation", {
  x <- scheme_wide(lower = c(0, 2, 4, 6), counts = c(1, 2, 3, 4))
  for (kind in c("lf", "wf")) {
    args <- setNames(list(c(2, 1), c(3, 4)), paste0(kind, c("_minbin", "_maxbin")))
    d <- do.call(get(paste0("prep_", kind, "_data")),
                 c(list(scheme_prep_data(), x), args))
    expect_identical(d[[paste0(kind, "_minbin")]], c(2L, 1L))
    expect_identical(d[[paste0(kind, "_maxbin")]], c(3L, 4L))
    expect_identical(d[[paste0(kind, "_nbins_i")]], 2L)
  }
})

test_that("native weight truncate selects intermediate rows without pooling tails", {
  d <- scheme_prep_data()
  d$lw_a <- 2
  x <- scheme_wide(lower = c(0, 2, 4, 6), counts = c(1, 2, 3, 4))
  expect_warning(a <- prep_wf_data(d, x, wf_minbin = c(2, 1),
    wf_maxbin = c(3, 4)), "narrower")
  expect_warning(t <- prep_wf_data(d, x, wf_minbin = c(2, 1),
    wf_maxbin = c(3, 4), wf_tails = "truncate"), "narrower")
  expect_identical(t$wf_rebin_s[[1]], t$wf_rebin_matrix[2:3, , drop = FALSE])
  expect_equal(a$wf_rebin_s[[1]], rbind(c(1, 0, 0, 0), c(0, 1, 0, 0)))
  expect_equal(t$wf_rebin_s[[1]], rbind(c(0.5, 0, 0, 0), c(0, 0.5, 0, 0)))
  expect_identical(t$wf_obs_flat, a$wf_obs_flat)
})

test_that("native weight bins warn when finer than overlapping source bins", {
  d <- scheme_prep_data()
  d$lw_a <- 2
  x <- scheme_wide(lower = c(0, 2, 4, 6), counts = c(1, 2, 3, 4))
  expect_warning(z <- prep_wf_data(d, x), "narrower than overlapping source")
  expect_identical(z$wf_native_i, TRUE)
  expect_equal(z$wf_rebin_s[[1]], z$wf_rebin_matrix)
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
      expect_identical(z[[paste0(kind, "_native_i")]], logical())
    }
  }
})
