# Age-based selectivity inputs for the Indian Ocean yellowfin example, shared
# by generate-io-yft-age-fits.R and generate-io-yft-vignette-cache.R. The
# io-yft-age-selectivity article shows the same steps.
#
# io_yft_age_inputs() returns the data with sel_basis_f = 2 for every
# fishery, par_sel converted from SS3's age_selex parameters, and two maps:
#   map_c2   C2's map: the same 49 selectivity parameters as the
#            length-based fit, so the basis is the only difference;
#   map_ss3  every selectivity parameter SS3 estimated (phase > 0): C2's
#            49 plus the nine C2 fixes.

# Unadjusted descending limb of opal's double normal at the last x, as in
# sel_double_normal()
dn_unadjusted_end <- function(par, x) {
  n <- length(x)
  bin_width <- x[2] - x[1]
  peak <- mean(x) + par[1] * sd(x)
  peak2 <- peak + bin_width + (0.99 * x[n] - peak - bin_width) / (1 + exp(-par[2]))
  exp(-(x[n] - peak2)^2 / (exp(par[4]) * sd(x)^2))
}

io_yft_age_inputs <- function() {
  d <- io_yft_data
  p <- io_yft_parameters
  ref <- io_yft_reference
  sp <- ref$ss3$sel_par
  pat <- ref$ss3$sel_pattern

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

  d$sel_basis_f <- rep(2L, d$n_fishery)
  p$par_sel <- convert_ss3_selex_to_rtmb(ss3_pars, d$sel_type_f,
                                         sel_lengths = d$len_mid,
                                         basis = "age", sel_ages = d$age_a)
  # SS3's final selectivity of -999 ("no adjustment") becomes -9 above.
  # Setting it to the unadjusted curve's own value at the oldest age
  # reproduces "no adjustment" exactly.
  for (f in which(d$sel_type_f == 2L & ss3_pars[, 6] <= -999)) {
    p$par_sel[f, 6] <- max(stats::qlogis(dn_unadjusted_end(p$par_sel[f, ], d$age_a)), -20)
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
  map_ss3 <- io_yft_map
  map_ss3$par_sel <- factor(as.vector(est))

  list(data = d, parameters = p, map_c2 = io_yft_map, map_ss3 = map_ss3,
       ss3_fleet = ss3_fleet, ss3_pars = ss3_pars)
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
# opal_check(), as in generate-io-yft-c2-fit.R.
io_yft_fit_polished <- function(x) {
  x <- opal_fit(x, n_passes = 2L, check = FALSE)
  object <- opal_rtmb(x, fresh = TRUE)
  opt <- x$fit$opt
  for (i in seq_len(5L)) {
    g <- as.numeric(object$gr(opt$par))
    if (max(abs(g)) < 1e-8) break
    opt$par <- stats::setNames(opt$par - as.numeric(solve(object$he(opt$par), g)),
                               names(opt$par))
  }
  opt$objective <- object$fn(opt$par)
  opal_attach_fit(x, opt, check = TRUE,
                  check_args = list(gradient_tolerance = 1e-4))
}
