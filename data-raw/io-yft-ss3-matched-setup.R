# SS3-matched inputs for the Indian Ocean yellowfin example, shared by
# generate-io-yft-matched-fits.R and generate-io-yft-vignette-cache.R. The
# io-yft-ss3-matched article shows the same steps.
#
# Starts from io_yft_age_inputs() (age-based selectivity at SS3's age_selex
# parameters, map_ss3 estimating every selectivity parameter SS3 estimated)
# and adds, one step at a time, the two SS3 settings fit A-SS3 lacks:
#
#   S1  SS3's lambda of 0 on parameter priors (lambdas$prior = 0);
#   S2  as S1, plus SS3's generalised size compositions for the index
#       (fleet 17, 14 bins from 30 cm with a 160 cm plus bin), on their own
#       bin scheme alongside the 5 cm length compositions, with SS3's lambda
#       of 0.51 on the four quarterly index fisheries;
#   S3  as S2 with the priors kept (A-SS3 plus the index size compositions),
#       which identifies every selectivity parameter SS3 estimated.
#
# Without the prior, some selectivity parameters are not identified: the
# likelihood is flat in them. generate-io-yft-matched-fits.R finds them and
# holds them at their estimates. SS3 estimated them, but was run without a
# Hessian.
#
# io_yft_matched_inputs() returns A-SS3's inputs and the data for S1-S3.
source("data-raw/io-yft-age-setup.R")

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

io_yft_matched_inputs <- function() {
  inputs <- io_yft_age_inputs()
  d <- inputs$data
  ref <- io_yft_reference
  index_f <- sort(unique(d$cpue_data$fishery))

  # S1: no priors, as SS3's lambda of 0 on parameter priors
  d_s1 <- d
  d_s1$lambdas <- list(prior = 0)

  # S2: add the index size compositions as a second bin scheme. SS3 puts the
  # population below 30 cm into the first bin (SIZEFREQ_TRANSLATION in
  # Report.sso), and the last bin is an open plus bin: opal's default
  # "accumulate" tails do both.
  top <- d$len_bin_start + d$n_len * d$len_bin_width
  d_s2 <- prep_lf_data(
    d,
    list(io_yft_lf_wide(d), io_yft_sizefreq_wide(ref)),
    lf_switch     = d$lf_switch,
    lf_minbin     = d$lf_minbin,
    lf_maxbin     = d$lf_maxbin,
    lf_var_adjust = d$lf_var_adjust,
    lf_addtocomp  = d$lf_addtocomp,
    lf_upper_edge = c(top, top))
  lf_lambda <- rep(1, d$n_fishery)
  lf_lambda[index_f] <- ref$ss3$sizefreq$lambda
  d_s2$lambdas <- list(prior = 0, lf = lf_lambda)
  d_s3 <- d_s2
  d_s3$lambdas <- list(lf = lf_lambda)

  c(inputs, list(data_s1 = d_s1, data_s2 = d_s2, data_s3 = d_s3, index_fisheries = index_f))
}
