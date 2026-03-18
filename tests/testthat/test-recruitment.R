test_that("Bias adjustment ramp generates correct scalars", {
  years <- 1900:2000
  bias_years <- c(1920, 1950, 1980, 1990)
  max_bias <- 0.8

  adj_vector <- get_bias_adj_vector(
    years = years,
    do_rec_bias_ramp = 1,
    bias_years = bias_years,
    max_bias_adj = max_bias
  )

  expect_equal(adj_vector[years == 1910], 0.0, tolerance = 1e-12)
  expect_equal(adj_vector[years == 1935], 0.4, tolerance = 1e-12)
  expect_equal(adj_vector[years == 1965], 0.8, tolerance = 1e-12)
  expect_equal(adj_vector[years == 1985], 0.4, tolerance = 1e-12)
  expect_equal(adj_vector[years == 1995], 0.0, tolerance = 1e-12)
})

test_that("Bias adjustment ramp defaults to full correction when off or missing", {
  years <- 1990:1995
  expect_equal(
    get_bias_adj_vector(years, do_rec_bias_ramp = 0, bias_years = c(1, 2, 3, 4), max_bias_adj = 0.8),
    rep(1.0, length(years))
  )
  expect_equal(
    get_bias_adj_vector(years, do_rec_bias_ramp = NULL, bias_years = c(1, 2, 3, 4), max_bias_adj = 0.8),
    rep(1.0, length(years))
  )
})

test_that("get_recruitment applies dynamic bias adjustment correctly", {
  sbio <- 1000
  rdev <- 0.0
  B0 <- 1000
  alpha <- 1.5
  beta <- 500
  sigma_r <- 0.5

  det_rec <- (alpha * sbio) / (beta + sbio)
  rec_no_adj <- get_recruitment(sbio, rdev, B0, alpha, beta, sigma_r, bias_adj = 0.0)
  expect_equal(rec_no_adj, det_rec, tolerance = 1e-12)

  rec_full_adj <- get_recruitment(sbio, rdev, B0, alpha, beta, sigma_r, bias_adj = 1.0)
  expect_equal(rec_full_adj, det_rec * exp(-0.5 * sigma_r^2), tolerance = 1e-12)

  rec_part_adj <- get_recruitment(sbio, rdev, B0, alpha, beta, sigma_r, bias_adj = 0.5)
  expect_equal(rec_part_adj, det_rec * exp(-0.25 * sigma_r^2), tolerance = 1e-12)
})
