# Inputs for the Indian Ocean yellowfin example, shared by
# generate-io-yft-fits.R and generate-io-yft-vignette-cache.R. The
# io-yft-ss3-comparison article builds AGE-SF with the same steps, inline.
#
# io_yft_inputs() returns the data, parameters and map of each fit in the
# article's bridge table. Each makes one change from the one before:
#
#   LEN        the bundled inputs: SS3's growth, opal's length-based
#              selectivity, io_yft_map (49 selectivity parameters);
#   AGE        age-based selectivity (sel_basis_f = 2) started from SS3's
#              age_selex parameters, estimating every selectivity parameter
#              SS3 estimated (58), with SS3's mirrors;
#   AGE-SF     AGE plus SS3's generalised size compositions for the index
#              (fleet 17, 14 bins from 30 cm with a 160 cm plus bin) on their
#              own bin scheme, with SS3's lambda of 0.51 on the four
#              quarterly index fisheries. The article's model;
#   AGE-SF-NP  AGE-SF with SS3's lambda of 0 on parameter priors, i.e. SS3's
#              configuration. Without the priors some selectivity parameters
#              are not identified; generate-io-yft-fits.R holds them fixed.
#
# The parameters returned are SS3's solution; generate-io-yft-fits.R starts
# AGE-SF and AGE-SF-NP from AGE's estimates instead.

# Unadjusted descending limb of opal's double normal at the last x, as in
# sel_double_normal()
dn_unadjusted_end <- function(par, x) {
  n <- length(x)
  bin_width <- x[2] - x[1]
  peak <- mean(x) + par[1] * sd(x)
  peak2 <- peak + bin_width + (0.99 * x[n] - peak - bin_width) / (1 + exp(-par[2]))
  exp(-(x[n] - peak2)^2 / (exp(par[4]) * sd(x)^2))
}

# The bundled 5 cm length compositions as one wide data frame of counts, in
# their input order. lf_n already carries SS3's variance adjustment
# (n / lf_var_adjust), so it is undone here for prep_lf_data() to reapply.
io_yft_lf_wide <- function(d = io_yft_data) {
  n_raw <- d$lf_n * d$lf_var_adjust[d$lf_fishery]
  counts <- d$lf_obs_in * n_raw
  colnames(counts) <- as.character(d$len_lower)
  data.frame(fishery = d$lf_fishery, year = d$lf_year, month = 1L, ts = d$lf_year,
             counts, check.names = FALSE)
}

# SS3's index size compositions as a wide data frame of counts
# (proportion x SS3 sample size), assigned to the quarterly index fishery of
# each time step. The last bin is SS3's plus bin: it ends at the top of the
# model's length grid.
io_yft_sizefreq_wide <- function(ref = io_yft_reference) {
  sf <- ref$ss3$sizefreq$data
  wide <- stats::reshape(sf[, c("ts", "fishery", "nsamp", "bin", "obs")],
                         idvar = c("ts", "fishery", "nsamp"), timevar = "bin",
                         direction = "wide")
  counts <- as.matrix(wide[, paste0("obs.", ref$ss3$sizefreq$bins)]) * wide$nsamp
  colnames(counts) <- as.character(ref$ss3$sizefreq$bins)
  data.frame(fishery = wide$fishery, year = wide$ts, month = 1L, ts = wide$ts,
             counts, check.names = FALSE, row.names = NULL)
}

io_yft_inputs <- function() {
  d <- io_yft_data
  p <- io_yft_parameters
  ref <- io_yft_reference
  sp <- ref$ss3$sel_par
  pat <- ref$ss3$sel_pattern

  # AGE ----------------------------------------------------------------------
  # SS3 fleet whose parameters each opal fishery uses: fisheries 1-16 are
  # SS3 fleets 1-16 (15 and 16 mirror 14 through pattern 15), and the four
  # quarterly index series 17-20 are SS3's survey, fleet 17.
  ss3_fleet <- c(seq_len(16), rep(17L, 4))
  mirror <- pat$pattern[ss3_fleet] == 15L
  ss3_fleet[mirror] <- pat$special[ss3_fleet[mirror]]

  ss3_pars <- t(vapply(ss3_fleet, function(fl) {
    v <- sp$value[sp$fleet == fl]
    c(v, rep(0, 6 - length(v)))
  }, numeric(6)))

  d_age <- d
  d_age$sel_basis_f <- rep(2L, d$n_fishery)
  p_age <- p
  p_age$par_sel <- convert_ss3_selex_to_rtmb(ss3_pars, d$sel_type_f,
                                             sel_lengths = d$len_mid,
                                             basis = "age", sel_ages = d$age_a)
  # SS3's final selectivity of -999 ("no adjustment") becomes -9 above.
  # Setting it to the unadjusted curve's own value at the oldest age
  # reproduces "no adjustment" exactly.
  for (f in which(d$sel_type_f == 2L & ss3_pars[, 6] <= -999)) {
    p_age$par_sel[f, 6] <- max(stats::qlogis(dn_unadjusted_end(p_age$par_sel[f, ], d$age_a)), -20)
  }

  # Every parameter SS3 estimated, with SS3's mirrors
  est <- matrix(NA_integer_, d$n_fishery, 6)
  k <- 0L
  for (fl in unique(ss3_fleet)) {
    cols <- sp$par[sp$fleet == fl & sp$phase > 0]
    ids <- k + seq_along(cols)
    k <- k + length(cols)
    for (f in which(ss3_fleet == fl)) est[f, cols] <- ids
  }
  map_age <- io_yft_map
  map_age$par_sel <- factor(as.vector(est))

  # AGE-SF -------------------------------------------------------------------
  # The index size compositions as a second bin scheme. SS3 puts the
  # population below 30 cm into the first bin (SIZEFREQ_TRANSLATION in
  # Report.sso), and the last bin is an open plus bin: opal's default
  # "accumulate" tails do both.
  index_f <- sort(unique(d$cpue_data$fishery))
  top <- d$len_bin_start + d$n_len * d$len_bin_width
  d_sf <- prep_lf_data(
    d_age,
    list(io_yft_lf_wide(d_age), io_yft_sizefreq_wide(ref)),
    lf_switch     = d_age$lf_switch,
    lf_minbin     = d_age$lf_minbin,
    lf_maxbin     = d_age$lf_maxbin,
    lf_var_adjust = d_age$lf_var_adjust,
    lf_addtocomp  = d_age$lf_addtocomp,
    lf_upper_edge = c(top, top))
  lf_lambda <- rep(1, d$n_fishery)
  lf_lambda[index_f] <- ref$ss3$sizefreq$lambda
  d_sf$lambdas <- list(lf = lf_lambda)

  # AGE-SF-NP: SS3's lambda of 0 on parameter priors -------------------------
  d_sf_np <- d_sf
  d_sf_np$lambdas <- list(prior = 0, lf = lf_lambda)

  list(
    LEN = list(data = d, parameters = p, map = io_yft_map),
    AGE = list(data = d_age, parameters = p_age, map = map_age),
    `AGE-SF` = list(data = d_sf, parameters = p_age, map = map_age),
    `AGE-SF-NP` = list(data = d_sf_np, parameters = p_age, map = map_age),
    ss3_fleet = ss3_fleet, ss3_pars = ss3_pars, index_fisheries = index_f)
}

# A compact numeric fingerprint of a model's inputs: the length and sum of
# every numeric, integer or logical element. The article computes it for the
# AGE-SF inputs it builds inline and compares it, with all.equal(), with the
# one the vignette cache stores, so that the article fails to render if its
# code no longer builds the stored fit's inputs.
io_yft_fingerprint <- function(data, parameters, map) {
  rapply(list(data = data, parameters = parameters, map = lapply(map, as.integer)),
         function(v) c(length(v), sum(as.numeric(v), na.rm = TRUE)),
         classes = c("numeric", "integer", "logical"), how = "unlist")
}

# opal's default bounds, with par_sel widened from +/-7 to +/-20 (SS3's
# bounds for P2-P4): on the age scale, SS3's plateau and ascending widths
# for several fisheries are near -10 to -15.
io_yft_age_bounds <- function(parameters, map, lim = 20) {
  active <- unlist(lapply(names(parameters), function(n)
    rep(n, if (is.null(map[[n]])) length(parameters[[n]]) else nlevels(droplevels(map[[n]])))))
  b <- get_bounds(list(par = stats::setNames(numeric(length(active)), active)), parameters)
  b$lower[b$parameter == "par_sel"] <- -lim
  b$upper[b$parameter == "par_sel"] <- lim
  b
}

# Two nlminb passes, Newton steps to take the gradient below 1e-8, and
# opal_check() with a gradient tolerance of 1e-4. With `guarded = TRUE` a
# Newton step is kept only if it lowers both the objective and the largest
# gradient: without priors the Hessian need not be positive definite, and an
# unguarded step can leave the optimum.
io_yft_fit_polished <- function(x, guarded = FALSE) {
  x <- opal_fit(x, n_passes = 2L, check = FALSE)
  object <- opal_rtmb(x, fresh = TRUE)
  opt <- x$fit$opt
  g <- as.numeric(object$gr(opt$par))
  f <- object$fn(opt$par)
  for (i in seq_len(5L)) {
    if (max(abs(g)) < 1e-8) break
    step <- tryCatch(as.numeric(solve(object$he(opt$par), g)),
                     error = function(e) if (guarded) NULL else stop(e))
    if (is.null(step)) break
    par <- stats::setNames(opt$par - step, names(opt$par))
    f_new <- object$fn(par)
    g_new <- as.numeric(object$gr(par))
    if (guarded && (!is.finite(f_new) || f_new > f || max(abs(g_new)) >= max(abs(g)))) break
    opt$par <- par
    f <- f_new
    g <- g_new
  }
  opt$objective <- f
  opal_attach_fit(x, opt, check = TRUE,
                  check_args = list(gradient_tolerance = 1e-4))
}

# One row of the article's bridge table, from a fitted model
io_yft_bridge_row <- function(x, model, change) {
  object <- opal_rtmb(x, fresh = TRUE)
  report <- object$report(x$fit$opt$par)
  n_ts <- x$data$n_year
  v <- x$validation$fit
  data.frame(model = model, change = change,
             sel_free = sum(grepl("^par_sel", names(x$fit$opt$par))),
             n_par = length(x$fit$opt$par),
             B0 = report$B0,
             terminal_depletion = report$static_depletion_y[n_ts],
             max_gradient = v$metrics$max_gradient,
             positive_hessian = isTRUE(v$metrics$positive_hessian),
             passes = isTRUE(v$passes))
}
