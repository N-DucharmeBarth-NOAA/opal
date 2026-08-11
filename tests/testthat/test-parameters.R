test_that("selectivity age initialisation pads outside old source range", {
  source_ages <- 2:9

  expect_equal(.match_selectivity_age_indices(2:9, source_ages), 2:9)
  expect_equal(.match_selectivity_age_indices(2:10, source_ages), c(2:9, 9))
  expect_equal(.match_selectivity_age_indices(0:3, source_ages), c(2, 2, 2, 3))
})

test_that("selectivity year initialisation pads outside old source range", {
  expect_equal(.match_selectivity_year_indices(1931:1933, 1931, 92), 1:3)
  expect_equal(.match_selectivity_year_indices(c(1930, 1931, 2022, 2023), 1931, 92),
               c(1L, 1L, 92L, 92L))
  expect_equal(.match_selectivity_year_indices(integer(), 1931, 92), integer())
})
