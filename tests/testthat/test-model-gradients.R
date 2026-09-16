gradient_free <- list(
  log_L1 = TRUE, log_L2 = TRUE, log_k = TRUE, log_CV1 = TRUE, log_CV2 = TRUE,
  log_h = TRUE, log_sigma_r = TRUE,
  log_cpue_tau = TRUE, log_cpue_omega = TRUE, cpue_creep = TRUE,
  log_init_F_f = 1L
)

safe_perturb <- function(par, scale = 0.05, phase = 0) {
  p <- par + scale * sin(seq_along(par) + phase)
  h <- names(p) == "log_h"
  p[h] <- pmin(p[h], log(0.95))
  creep <- names(p) == "cpue_creep"
  p[creep] <- 0
  p
}

fd5 <- function(fn, par, rel = 1e-4) {
  vapply(seq_along(par), function(i) {
    h <- rel * max(1, abs(par[i]))
    e <- replace(numeric(length(par)), i, h)
    (-fn(par + 2 * e) + 8 * fn(par + e) - 8 * fn(par - e) + fn(par - 2 * e)) /
      (12 * h)
  }, numeric(1))
}

expect_gradient_match <- function(ga, gf, nms, abs_tol = 1e-6, rel_tol = 1e-5) {
  excess <- abs(ga - gf) - (abs_tol + rel_tol * abs(gf))
  worst <- head(order(excess, decreasing = TRUE), 5)
  expect_true(all(excess <= 0), info = paste(
    sprintf("%s: ad=%.8g fd=%.8g", nms[worst], ga[worst], gf[worst]),
    collapse = "; "
  ))
}

expect_nonzero_groups <- function(g, nms, groups, tol = 1e-8) {
  for (grp in intersect(groups, nms)) {
    expect_true(any(abs(g[nms == grp]) > tol),
                label = sprintf("gradient for %s is identically zero", grp))
  }
}

test_that("recorded tape matches a fresh tape at perturbed parameters", {
  inputs <- opaka_inputs()
  map <- free_in_map(inputs$map, inputs$parameters, gradient_free)
  obj <- opaka_obj(inputs, map = map)
  for (k in 1:3) {
    p <- safe_perturb(obj$par, phase = k)
    fresh <- opaka_obj(inputs, parameters = obj$env$parList(p), map = map)
    lab <- sprintf("perturbation %d", k)
    expect_equal(obj$fn(p), fresh$fn(), tolerance = 1e-10, label = lab)
    expect_equal(obj$report(p)$spawning_biomass_y,
                 fresh$report()$spawning_biomass_y,
                 tolerance = 1e-10, label = lab)
  }
})

test_that("AD gradient matches finite differences (lf_switch = 1)", {
  inputs <- opaka_inputs()
  map <- free_in_map(inputs$map, inputs$parameters, gradient_free)
  obj <- opaka_obj(inputs, map = map)
  p <- safe_perturb(obj$par)
  ga <- as.vector(obj$gr(p))
  gf <- fd5(obj$fn, p)
  expect_gradient_match(ga, gf, names(p))
  expect_nonzero_groups(ga, names(p), c(
    "log_L1", "log_L2", "log_k", "log_CV1", "log_CV2", "log_h",
    "log_sigma_r", "log_cpue_tau", "log_cpue_omega", "log_init_F_f",
    "par_sel"
  ))
})

for (sw in 2:3) {
  test_that(sprintf("AD gradient matches finite differences (lf_switch = %d)", sw), {
    inputs <- opaka_inputs()
    data <- inputs$data
    data$lf_switch <- sw
    map <- free_in_map(inputs$map, inputs$parameters, list(
      log_lf_tau = c(1L, 3L), log_L1 = TRUE, log_k = TRUE, log_CV2 = TRUE
    ))
    map$rdev_y <- factor(rep(NA, length(inputs$parameters$rdev_y)))
    map$init_rdev_a <- factor(rep(NA, length(inputs$parameters$init_rdev_a)))
    obj <- opaka_obj(inputs, data = data, map = map)
    p <- safe_perturb(obj$par)
    ga <- as.vector(obj$gr(p))
    expect_gradient_match(ga, fd5(obj$fn, p), names(p))
    expect_nonzero_groups(ga, names(p), c("log_lf_tau", "log_L1", "log_k", "log_CV2"))
  })
}