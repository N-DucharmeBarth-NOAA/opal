# Tests for rebin_counts() and rebin_matrix()
library(RTMB)

# 1. Identity rebin ----
test_that("rebin_counts: identity rebin returns input unchanged", {
  edges  <- 0:5
  counts <- c(10, 20, 15, 5, 30)
  result <- rebin_counts(edges, counts, edges)
  expect_equal(result, counts)
})

test_that("rebin_matrix: identity rebin gives identity matrix", {
  edges <- 0:5
  W <- rebin_matrix(edges, edges)
  expect_equal(W, diag(5))
})

# 2. 2:1 coarsening ----
test_that("rebin_counts: 2:1 coarsening sums pairs of bins", {
  src_edges  <- 0:4
  src_counts <- c(10, 20, 15, 5)
  dest_edges <- c(0, 2, 4)
  result <- rebin_counts(src_edges, src_counts, dest_edges)
  expect_equal(result, c(30, 20))
})

# 3. 1:2 refinement ----
test_that("rebin_counts: 1:2 refinement splits each bin in half", {
  src_edges  <- c(0, 2, 4)
  src_counts <- c(10, 20)
  dest_edges <- 0:4
  result <- rebin_counts(src_edges, src_counts, dest_edges)
  expect_equal(result, c(5, 5, 10, 10))
})

# 4. Irregular to regular ----
test_that("rebin_counts: irregular source to regular destination", {
  src_edges  <- c(0, 1, 4, 10)
  src_counts <- c(6, 9, 12)
  dest_edges <- c(0, 2, 4, 6, 8, 10)
  result <- rebin_counts(src_edges, src_counts, dest_edges)
  # bin [0,1]: 6/1 * 1=6 from src[1], 0 from others → dest[1] gets 6 + 9*(1/3)=9
  # [0,2]: 6*1 + 9*(1/3)=9  → 6 + 3 = 9
  # [2,4]: 9*(2/3) = 6
  # [4,6]: 12*(2/6)=4
  # [6,8]: 12*(2/6)=4
  # [8,10]: 12*(2/6)=4
  expect_equal(result, c(9, 6, 4, 4, 4))
  expect_equal(sum(result), sum(src_counts))
})

# 5. Partial overlap ----
test_that("rebin_counts: partial overlap — sum(dest) < sum(src)", {
  src_edges  <- c(0, 5, 10)
  src_counts <- c(20, 30)
  dest_edges <- c(2, 8)  # narrower than source
  result <- rebin_counts(src_edges, src_counts, dest_edges)
  expect_true(sum(result) < sum(src_counts))
  # [2,8]: from [0,5] → overlap [2,5]=3/5 * 20=12; from [5,10] → overlap [5,8]=3/5 * 30=18
  expect_equal(result, c(30))
})

# 6. No overlap ----
test_that("rebin_counts: disjoint source and dest gives all zeros", {
  src_edges  <- c(0, 5)
  src_counts <- c(100)
  dest_edges <- c(10, 20)
  result <- rebin_counts(src_edges, src_counts, dest_edges)
  expect_equal(result, c(0))
})

# 7. Single source bin ----
test_that("rebin_counts: single source bin split across many dest bins", {
  src_edges  <- c(0, 10)
  src_counts <- c(100)
  dest_edges <- 0:10
  result <- rebin_counts(src_edges, src_counts, dest_edges)
  expect_equal(result, rep(10, 10))
  expect_equal(sum(result), 100)
})

# 8. Single destination bin ----
test_that("rebin_counts: many source bins collapsed into one dest bin", {
  src_edges  <- 0:5
  src_counts <- c(10, 20, 15, 5, 30)
  dest_edges <- c(0, 5)
  result <- rebin_counts(src_edges, src_counts, dest_edges)
  expect_equal(result, c(sum(src_counts)))
})

# 9. Matrix equivalence ----
test_that("rebin_matrix %*% counts equals rebin_counts for all cases", {
  cases <- list(
    list(src = 0:5,         cnt = c(10, 20, 15, 5, 30), dst = 0:5),
    list(src = 0:4,         cnt = c(10, 20, 15, 5),     dst = c(0, 2, 4)),
    list(src = c(0, 2, 4),  cnt = c(10, 20),             dst = 0:4),
    list(src = c(0, 1, 4, 10), cnt = c(6, 9, 12),       dst = c(0, 2, 4, 6, 8, 10)),
    list(src = c(0, 10),    cnt = c(100),                dst = 0:10)
  )
  for (case in cases) {
    direct <- rebin_counts(case$src, case$cnt, case$dst)
    via_W  <- as.vector(rebin_matrix(case$src, case$dst) %*% case$cnt)
    expect_equal(via_W, direct)
  }
})

# 10. Mass conservation ----
test_that("rebin_counts: mass conserved when dest covers full src range", {
  src_edges  <- c(0, 1, 3, 6, 10)
  src_counts <- c(5, 15, 10, 20)
  dest_edges <- c(0, 2, 5, 10)
  result <- rebin_counts(src_edges, src_counts, dest_edges)
  expect_equal(sum(result), sum(src_counts))
})

# 11. Realistic L-W case ----
test_that("rebin_counts: realistic L-W rebinning conserves total counts", {
  # Length bin edges: 2 cm bins from 10 to 200 cm (95 bins, 96 edges)
  len_edges <- seq(10, 200, by = 2)
  lw_a <- 2.5e-5
  lw_b <- 2.9
  # Weight edges from L-W transform
  wt_edges <- lw_a * len_edges^lw_b
  # Destination: regular 1 kg bins, 0 to ceil(max wt)
  max_wt_kg <- ceiling(max(wt_edges))
  dest_edges <- seq(0, max_wt_kg, by = 1)
  # Fake length-bin counts (random positive values)
  set.seed(42)
  src_counts <- runif(length(len_edges) - 1, 0, 100)
  result <- rebin_counts(wt_edges, src_counts, dest_edges)
  expect_equal(length(result), length(dest_edges) - 1L)
  expect_equal(sum(result), sum(src_counts), tolerance = 1e-8)
})

# 12. AD tape test ----
test_that("rebin_counts works inside RTMB::MakeADFun with finite fn and gr", {
  src_edges  <- c(0, 1, 3, 6, 10)
  dest_edges <- c(0, 2, 5, 10)
  N <- length(src_edges) - 1L

  f <- function(par) {
    sum(rebin_counts(src_edges, par$log_counts, dest_edges))
  }

  par_init <- list(log_counts = log(c(5, 15, 10, 20)))
  obj <- RTMB::MakeADFun(f, par_init, silent = TRUE)

  fn_val <- obj$fn()
  gr_val <- obj$gr()

  expect_true(is.finite(fn_val))
  expect_true(all(is.finite(gr_val)))
})

# 13. Edge validation ----
test_that("rebin_counts errors on non-monotonic src_edges", {
  expect_error(
    rebin_counts(c(0, 5, 3, 10), c(1, 2, 3), c(0, 5, 10)),
    "src_edges must be monotonically increasing"
  )
})

test_that("rebin_counts errors on non-monotonic dest_edges", {
  expect_error(
    rebin_counts(c(0, 5, 10), c(1, 2), c(0, 10, 5)),
    "dest_edges must be monotonically increasing"
  )
})

test_that("rebin_counts errors when src_counts length doesn't match src_edges", {
  expect_error(
    rebin_counts(c(0, 5, 10), c(1, 2, 3), c(0, 5, 10)),
    "length\\(src_counts\\) must equal length\\(src_edges\\) - 1"
  )
})

test_that("rebin_matrix errors on non-monotonic src_edges", {
  expect_error(
    rebin_matrix(c(0, 5, 3, 10), c(0, 5, 10)),
    "src_edges must be monotonically increasing"
  )
})

test_that("rebin_matrix errors on non-monotonic dest_edges", {
  expect_error(
    rebin_matrix(c(0, 5, 10), c(0, 10, 5)),
    "dest_edges must be monotonically increasing"
  )
})
