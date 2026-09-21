test_that("posfun penalty is nonnegative and increases under violation", {
  safe <- posfun(1, eps = 0.001)
  near <- posfun(0.001, eps = 0.001)
  violated <- posfun(-0.1, eps = 0.001)

  expect_gte(safe$penalty, 0)
  expect_gte(near$penalty, 0)
  expect_gte(violated$penalty, 0)
  expect_lt(safe$penalty, 1e-12)
  expect_gt(violated$penalty, near$penalty)
})
