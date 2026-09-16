ref_multinom_nll <- function(x, p) {
  -(lgamma(sum(x) + 1) - sum(lgamma(x + 1)) + sum(x * log(p)))
}

ref_dirichlet_nll <- function(x, alpha) {
  -(lgamma(sum(alpha)) - sum(lgamma(alpha)) + sum((alpha - 1) * log(x)))
}

ref_dirmult_nll <- function(x, alpha, size = sum(x)) {
  alpha_sum <- sum(alpha)
  -(lgamma(size + 1) - sum(lgamma(x + 1)) + lgamma(alpha_sum) -
      lgamma(size + alpha_sum) + sum(lgamma(x + alpha) - lgamma(alpha)))
}

ref_comp_pred <- function(catch_a, pla, bmin, bmax, addtocomp = 1e-8,
                          rebin = NULL) {
  pred <- as.vector(pla %*% catch_a)
  if (!is.null(rebin)) pred <- as.vector(rebin %*% pred)
  n <- length(pred)
  if (bmin > 1) pred[bmin] <- sum(pred[1:bmin])
  if (bmax < n) pred[bmax] <- sum(pred[bmax:n])
  pred <- pred[bmin:bmax] + addtocomp
  pred / sum(pred)
}

comp_fixture <- function(n_fishery = 2L, n_year = 4L, n_age = 6L,
                        n_len = 12L) {
  ages <- seq_len(n_age)
  lens <- seq_len(n_len)
  pla <- sapply(ages, function(age) {
    p <- dnorm(lens, 1.6 * age + 1, 1.2)
    p / sum(p)
  })
  catch <- array(0, c(n_fishery, n_year, n_age))
  for (f in seq_len(n_fishery)) for (y in seq_len(n_year)) {
    catch[f, y, ] <- exp(-0.3 * ages) * (1 + 0.2 * f) *
      (1 + 0.1 * y) * plogis(ages - 2)
  }
  list(pla = pla, catch = catch, n_fishery = n_fishery, n_len = n_len)
}

base_groups <- function() list(
  list(f = 1L, ys = c(1L, 3L), bmin = 1L, bmax = 12L, n = c(80, 120)),
  list(f = 2L, ys = c(2L, 3L, 4L), bmin = 3L, bmax = 9L, n = c(40, 60, 50))
)

attach_obs <- function(groups, integer = FALSE) {
  lapply(seq_along(groups), function(k) {
    group <- groups[[k]]
    n_bins <- group$bmax - group$bmin + 1L
    group$obs <- t(vapply(seq_along(group$ys), function(i) {
      weights <- 1 + 0.5 * sin(seq_len(n_bins) + 3 * i + k)
      values <- group$n[i] * weights / sum(weights)
      if (integer) {
        values <- floor(values)
        values[which.max(values)] <- values[which.max(values)] +
          group$n[i] - sum(values)
      }
      values
    }, numeric(n_bins)))
    group
  })
}

comp_args <- function(prefix, groups, fx, type, log_tau, rebin = NULL,
                      removal_switch_f = rep(0L, fx$n_fishery),
                      addtocomp = 1e-8) {
  n_bins <- if (is.null(rebin)) nrow(fx$pla) else nrow(rebin)
  name <- function(suffix) paste0(prefix, "_", suffix)
  minbin <- rep(1L, fx$n_fishery)
  maxbin <- rep(n_bins, fx$n_fishery)
  for (group in groups) {
    minbin[group$f] <- group$bmin
    maxbin[group$f] <- group$bmax
  }
  flat <- unlist(lapply(groups, function(group) as.vector(t(group$obs))))
  prop <- unlist(lapply(groups, function(group) {
    as.vector(t(group$obs / rowSums(group$obs)))
  }))

  args <- list(catch_pred_fya = fx$catch, pla = fx$pla,
               removal_switch_f = removal_switch_f)
  args[[name("obs_flat")]] <- flat
  args[[name("obs_ints")]] <- as.integer(round(flat))
  args[[name("obs_prop")]] <- prop
  args[[name("n_f")]] <- vapply(groups, function(group) length(group$ys), integer(1))
  args[[name("fishery_f")]] <- vapply(groups, `[[`, integer(1), "f")
  args[[name("year_fi")]] <- lapply(groups, `[[`, "ys")
  args[[name("n_fi")]] <- lapply(groups, `[[`, "n")
  args[[name("minbin")]] <- minbin
  args[[name("maxbin")]] <- maxbin
  args[[name("switch")]] <- as.integer(type)
  args[[name("addtocomp")]] <- addtocomp
  args[[paste0("n_", prefix)]] <- sum(args[[name("n_f")]])
  args[[paste0("log_", prefix, "_tau")]] <- log_tau
  if (prefix == "lf") {
    args$n_len <- n_bins
  } else {
    args$n_wt <- n_bins
    args$wf_rebin_matrix <- rebin
  }
  args
}

expected_comp_nll <- function(groups, fx, type, log_tau, rebin = NULL,
                              removal_switch_f = rep(0L, fx$n_fishery),
                              addtocomp = 1e-8) {
  unlist(lapply(groups, function(group) vapply(seq_along(group$ys), function(i) {
    if (removal_switch_f[group$f] == 1L || group$n[i] == 0) return(0)
    pred <- ref_comp_pred(fx$catch[group$f, group$ys[i], ], fx$pla,
                          group$bmin, group$bmax, addtocomp, rebin)
    x <- group$obs[i, ]
    tau <- exp(log_tau[group$f])
    switch(type,
           ref_multinom_nll(x, pred),
           ref_dirichlet_nll(x / sum(x), pred * group$n[i] * tau),
           ref_dirmult_nll(round(x), pred * tau, size = group$n[i]))
  }, numeric(1))), use.names = FALSE)
}