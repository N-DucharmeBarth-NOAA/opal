test_that("parameter tables distinguish fixed, shared, and boundary parameters", {
  parameters <- list(alpha = c(0.01, 0.99, 0.5), beta = c(3, 3, 7),
                     rdev_y = 0.2, init_rdev_a = 0.1)
  map <- list(beta = factor(c(1, 1, NA)))
  obj <- RTMB::MakeADFun(function(p) {
    sum(p$alpha^2) + sum(p$beta^2) + sum(p$rdev_y^2) + sum(p$init_rdev_a^2)
  }, parameters,
                        map = map, silent = TRUE)
  obj$fn(obj$par)
  lower <- c(0, 0, 0, 0, -Inf, -Inf)
  upper <- c(1, 1, 1, 10, Inf, Inf)
  tab <- get_par_table(obj, parameters, map, lower, upper,
                       include = "all", digits = NULL, show_map = TRUE)
  expect_equal(tab$est, unname(unlist(parameters)))
  expect_equal(tab$bd_chk[1:3], c("LO", "HI", "OK"))
  expect_equal(tab$gr[1:3], 2 * parameters$alpha)
  expect_equal(tab$gr[4], 12) # Both beta elements share one estimated value.
  expect_true(is.na(tab$gr[5]))
  expect_true(tab$fixed[6])
  expect_true(is.na(tab$gr_chk[6]))
  expect_equal(tab$map[4:5], c("1", "1"))
  core <- get_par_table(obj, parameters, map)
  expect_false(any(grepl("rdev", core$par)))
  expect_false(any(c("map", "fixed", "group") %in% names(core)))
  estimated <- get_par_table(obj, parameters, map, include = "all_est")
  expect_false("beta3" %in% estimated$par)
  expect_true(all(c("rdev_y", "init_rdev_a") %in% estimated$par))
})

test_that("one-sided bounds and gradient thresholds are reported correctly", {
  p <- list(theta = c(0.01, 0.99, 5))
  obj <- RTMB::MakeADFun(function(p) sum(p$theta^2), p, silent = TRUE)
  obj$fn(obj$par)
  tab <- get_par_table(obj, p, list(), lower = c(0, -Inf, -Inf),
                       upper = c(Inf, 1, Inf), grad_tol = 0.1)
  expect_equal(tab$bd_chk, c("LO", "HI", "OK"))
  expect_equal(tab$gr_chk, c("OK", "BAD", "BAD"))
  unbounded <- get_par_table(obj, p, list())
  expect_equal(unbounded$bd_chk, rep("OK", 3))
  expect_equal(unbounded$lwr, rep(-Inf, 3))
  expect_equal(unbounded$upr, rep(Inf, 3))
})

test_that("correlation diagnostics recover known positive and negative pairs", {
  covariance <- matrix(c(1, 0.95, 0.95, 1), 2)
  precision <- solve(covariance)
  p <- list(theta = c(0, 0), rdev_y = 0)
  obj <- RTMB::MakeADFun(function(p) {
    sum(p$theta * (precision %*% p$theta)) / 2 + p$rdev_y^2 / 2
  }, p, silent = TRUE)
  obj$fn(obj$par)
  pairs <- get_cor_pairs(obj, digits = NULL)
  expect_equal(pairs$par1, "theta1")
  expect_equal(pairs$par2, "theta2")
  expect_equal(pairs$correlation, 0.95, tolerance = 1e-8)
  expect_equal(get_cor_pairs(obj, include = "all"), pairs)
  expect_match(get_cor_pairs(obj, threshold = 0.99), "No pairs")
  numerical <- obj
  numerical["he"] <- list(NULL)
  expect_equal(get_cor_pairs(numerical)$correlation, 0.95)
  negative <- diag(3)
  negative[1:2, 1:2] <- solve(matrix(c(1, -0.96, -0.96, 1), 2))
  expect_equal(get_cor_pairs(obj, h = negative)$correlation, -0.96)
  expect_error(get_cor_pairs(obj, h = matrix(0, 3, 3)), "Hessian is singular")
})

test_that("natural mortality plot uses the model's age-specific report", {
  x <- opal_build(small_opal_object())
  object <- opal_rtmb(x)
  plot <- plot_natural_mortality(x$data, object)
  layer <- ggplot2::ggplot_build(plot)$data[[1]]
  expect_equal(layer$x, x$data$min_age:x$data$max_age)
  expect_equal(layer$y, object$report()$M_a)
  expect_true(all(layer$linetype == "dashed"))
  expect_equal(plot$labels$y, "Natural mortality")
})
