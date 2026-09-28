small_opal_object <- function(random = FALSE) {
  data <- make_opal_fit_data()
  parameters <- make_opal_fit_parameters(data)
  parameters$log_h <- log(0.75)
  parameters$log_B0 <- 15
  map <- make_opal_fit_map(parameters)
  if (random) map$rdev_y <- NULL
  opal_obj(data, parameters, map, random = if (random) "rdev_y" else character(),
    bounds = list(lower = c(log_B0 = 13), upper = c(log_B0 = 22)))
}

mock_opal_sampler <- function(obj, chains = 2L, num_samples = 400L,
                              num_warmup = 10L, seed = 42L, ...) {
  set.seed(seed)
  n <- num_samples + num_warmup
  vars <- opal:::.opal_expand_parameter_names(names(obj$par))
  samples <- array(rnorm(n * chains * length(vars), sd = 0.01), c(n, chains, length(vars)),
                    dimnames = list(NULL, NULL, vars))
  samples <- sweep(samples, 3, obj$par, "+")
  list(samples = samples, warmup = num_warmup, iter = num_samples,
    max_treedepth = 10L, sampler_params = replicate(chains,
      cbind(divergent__ = rep(0, n), treedepth__ = rep(3, n)), simplify = FALSE))
}
