# Fit the SS3-matched Indian Ocean yellowfin models for the io-yft-ss3-matched
# article. Each adds one SS3 setting to fit A-SS3 (age-based selectivity,
# every selectivity parameter SS3 estimated). Run from the package root;
# takes about 30 minutes on 4 cores. Afterwards, rerun
# generate-io-yft-vignette-cache.R.
#
#   inst/extdata/io_yft_matched_s1_fit.rds   S1: no priors (SS3's lambda of 0
#                                            on parameter priors)
#   inst/extdata/io_yft_matched_s2_fit.rds   S2: S1 plus SS3's index size
#                                            compositions, lambda 0.51
#   inst/extdata/io_yft_matched_s1_free.rds  identifiability check: S1 with
#                                            A-SS3's map, unpolished
#   inst/extdata/io_yft_matched_starts.rds   start sensitivity: S1 and S2
#                                            from SS3's solution, from each
#                                            other's estimates and from four
#                                            jittered starts
#
# S1 and S2 start from A-SS3's estimates. Started from SS3's own solution
# they stop at a local optimum about 20 likelihood units worse, with a
# larger stock; the start-sensitivity fits record this.
pkgload::load_all()
source("data-raw/io-yft-ss3-matched-setup.R")

inputs <- io_yft_matched_inputs()
control <- list(eval.max = 10000, iter.max = 10000)
metadata <- list(
  stock = "Indian Ocean yellowfin (simulated, replicate 4)",
  source = "SS3 model by M.N. Maunder, H. Xu and C. Minte-Vera (IATTC), converted to opal",
  converter = "https://github.com/N-DucharmeBarth/wcpfc-bet-single-region (code/rtmb)")
build <- function(data, parameters, map, model = NULL) {
  opal_obj(data, parameters, map, bounds = io_yft_age_bounds(parameters, map),
           control = control, metadata = c(metadata, list(model = model)))
}
# Parameter list at an estimate
par_list <- function(x, par) opal_rtmb(x, fresh = TRUE)$env$parList(par)
n_ts <- inputs$data$n_year
# Two nlminb passes, then Newton steps only while they lower both the
# objective and the largest gradient: without priors the Hessian need not be
# positive definite, and an unguarded Newton step can leave the optimum.
fit_polished_safe <- function(x) {
  x <- opal_fit(x, n_passes = 2L, check = FALSE)
  object <- opal_rtmb(x, fresh = TRUE)
  opt <- x$fit$opt
  g <- as.numeric(object$gr(opt$par))
  f <- object$fn(opt$par)
  for (i in seq_len(5L)) {
    if (max(abs(g)) < 1e-8) break
    step <- tryCatch(as.numeric(solve(object$he(opt$par), g)), error = function(e) NULL)
    if (is.null(step)) break
    par <- stats::setNames(opt$par - step, names(opt$par))
    f_new <- object$fn(par)
    g_new <- as.numeric(object$gr(par))
    if (!is.finite(f_new) || f_new > f || max(abs(g_new)) >= max(abs(g))) break
    opt$par <- par
    f <- f_new
    g <- g_new
  }
  opt$objective <- f
  opal_attach_fit(x, opt, check = TRUE, check_args = list(gradient_tolerance = 1e-4))
}
summarise_fit <- function(x, opt) {
  object <- opal_rtmb(x, fresh = TRUE)
  report <- object$report(opt$par)
  list(objective = opt$objective, B0 = report$B0,
       max_gradient = max(abs(object$gr(opt$par))),
       ssb = report$spawning_biomass_y[seq_len(n_ts)],
       static_dep = report$static_depletion_y[seq_len(n_ts)])
}

# A-SS3's estimates: the starting point of S1 and S2
a_ss3 <- opal_read(file.path("inst", "extdata", "io_yft_age_ss3_fit.rds"))
p_a <- par_list(a_ss3, a_ss3$fit$opt$par)

data_m <- list(S1 = inputs$data_s1, S2 = inputs$data_s2)
model_m <- c(S1 = "S1: as A-SS3, with SS3's lambda of 0 on parameter priors",
             S2 = "S2: as S1, with SS3's index size compositions (lambda 0.51)")

# Phase 1, in parallel: S1 and S2 ------------------------------------------
# Without the prior, several selectivity parameters are not identified: the
# likelihood is flat in them. Each model is fitted freely first (A-SS3's
# map, every parameter SS3 estimated); parameters with a Hessian diagonal
# below 1e-3 at that estimate are then held at it, which leaves the
# objective unchanged, and any parameter still dominating a near-zero
# eigenvalue is added, before the final fit is polished and checked.
sel_labels <- function(map) {
  lv <- matrix(as.integer(map$par_sel), nrow = inputs$data$n_fishery)
  first <- vapply(seq_len(max(lv, na.rm = TRUE)), function(k) which(lv == k)[1], integer(1))
  list(lv = lv, first = first, pos = arrayInd(first, dim(lv)))
}
hold <- function(map, cells) {
  lv <- matrix(as.integer(map$par_sel), nrow = inputs$data$n_fishery)
  lv[lv %in% lv[cells]] <- NA_integer_
  map$par_sel <- factor(as.vector(lv))
  map
}
fit_identified <- function(job) {
  x <- opal_fit(build(data_m[[job]], p_a, inputs$map_ss3), n_passes = 2L, check = FALSE)
  object <- opal_rtmb(x, fresh = TRUE)
  opt <- x$fit$opt
  hess <- object$he(opt$par)
  sl <- sel_labels(inputs$map_ss3)
  sel <- grep("^par_sel", names(opt$par))
  flat <- data.frame(fishery = sl$pos[, 1], col = sl$pos[, 2],
                     start = p_a$par_sel[sl$first], estimate = unname(opt$par[sel]),
                     hessian = diag(hess)[sel])
  free <- list(objective = opt$objective, n_par = length(opt$par),
               max_gradient = max(abs(object$gr(opt$par))),
               n_flat_eigen = sum(abs(eigen(hess, symmetric = TRUE,
                                            only.values = TRUE)$values) < 1e-6))
  p_free <- par_list(x, opt$par)
  map <- hold(inputs$map_ss3, sl$first[abs(flat$hessian) < 1e-3])
  for (round in 1:3) {
    y <- opal_fit(build(data_m[[job]], p_free, map, model_m[[job]]), n_passes = 2L,
                  check = FALSE)
    o <- opal_rtmb(y, fresh = TRUE)
    e <- eigen(o$he(y$fit$opt$par), symmetric = TRUE)
    small <- which(abs(e$values) < 1e-6)
    if (!length(small)) break
    sly <- sel_labels(map)
    sel_y <- grep("^par_sel", names(y$fit$opt$par))
    add <- vapply(small, function(k) sel_y[which.max(abs(e$vectors[sel_y, k]))], integer(1))
    p_free <- par_list(y, y$fit$opt$par)
    map <- hold(map, sly$first[match(add, sel_y)])
  }
  held_cells <- which(!is.na(as.integer(inputs$map_ss3$par_sel)) & is.na(as.integer(map$par_sel)))
  flat$held <- sl$first %in% held_cells |
    vapply(seq_len(nrow(flat)), function(k) any(sl$lv[held_cells] == k), logical(1))
  free$flat <- flat
  y <- fit_polished_safe(build(data_m[[job]], p_free, map, model_m[[job]]))
  cat(sprintf("%s: objective %.4f (free %.4f), %d parameters, %d held, check %s\n", job,
              y$fit$opt$objective, free$objective, length(y$fit$opt$par), sum(flat$held),
              isTRUE(y$validation$fit$passes)))
  list(fit = y, free = free)
}
phase1 <- parallel::mclapply(c(S1 = "S1", S2 = "S2"), fit_identified,
                             mc.cores = 2L, mc.preschedule = FALSE)
for (nm in c("S1", "S2")) {
  if (inherits(phase1[[nm]], "try-error")) stop(nm, ": ", phase1[[nm]])
  # Report the fit check rather than stop on it: without priors the Hessian
  # need not be positive definite (see the article)
  v <- phase1[[nm]]$fit$validation$fit
  cat(sprintf("%s fit check: passes %s, max gradient %.2e, positive Hessian %s\n", nm,
              isTRUE(v$passes), v$metrics$max_gradient, isTRUE(v$metrics$positive_hessian)))
  opal_save(phase1[[nm]]$fit, file.path("inst", "extdata",
                                       sprintf("io_yft_matched_%s_fit.rds", tolower(nm))),
            compress = "xz", overwrite = TRUE)
}
saveRDS(phase1$S1$free, file.path("inst", "extdata", "io_yft_matched_s1_free.rds"),
        compress = "xz", version = 3L)
fitted <- lapply(phase1, `[[`, "fit")

# Phase 2, in parallel: start sensitivity ---------------------------------
# Same data and map as each fit; unpolished. Starts: SS3's solution, the
# other model's estimates, and four jitters of the fit's own estimates
# (normal, SD 0.2, on the estimated selectivity parameters and log B0).
# Held parameters stay at the fit's values in every start.
fitted_par <- lapply(fitted, function(x) par_list(x, x$fit$opt$par))
start_from <- function(p, nm) {
  out <- p
  est <- !is.na(fitted[[nm]]$map$par_sel)
  out$par_sel <- fitted_par[[nm]]$par_sel
  out$par_sel[est] <- p$par_sel[est]
  out
}
jitter <- function(nm, seed, sd = 0.2) {
  p <- fitted_par[[nm]]
  est <- !is.na(fitted[[nm]]$map$par_sel)
  set.seed(seed)
  p$par_sel[est] <- p$par_sel[est] + stats::rnorm(sum(est), 0, sd)
  p$log_B0 <- p$log_B0 + stats::rnorm(1, 0, sd)
  p
}
starts <- list(
  S1 = c(list(ss3 = start_from(inputs$parameters, "S1"), cross = start_from(fitted_par$S2, "S1")),
         stats::setNames(lapply(1:4, function(i) jitter("S1", 100 + i)), paste0("jitter", 1:4))),
  S2 = c(list(ss3 = start_from(inputs$parameters, "S2"), cross = start_from(fitted_par$S1, "S2")),
         stats::setNames(lapply(1:4, function(i) jitter("S2", 200 + i)), paste0("jitter", 1:4))))
jobs <- do.call(rbind, lapply(names(starts), function(m)
  data.frame(model = m, start = names(starts[[m]]))))
runs <- parallel::mclapply(seq_len(nrow(jobs)), function(i) {
  m <- jobs$model[i]
  x <- build(data_m[[m]], starts[[m]][[jobs$start[i]]], fitted[[m]]$map)
  x <- tryCatch(opal_fit(x, n_passes = 2L, check = FALSE), error = function(e) e)
  if (inherits(x, "error")) return(list(error = conditionMessage(x)))
  s <- summarise_fit(x, x$fit$opt)
  cat(sprintf("%s from %s: objective %.3f, B0 %.4g\n", m, jobs$start[i], s$objective, s$B0))
  s
}, mc.cores = 4L, mc.preschedule = FALSE)
reference <- lapply(fitted, function(x) summarise_fit(x, x$fit$opt))
saveRDS(list(jobs = jobs, runs = runs, reference = reference, jitter_sd = 0.2),
        file.path("inst", "extdata", "io_yft_matched_starts.rds"),
        compress = "xz", version = 3L)
