#' Prepare length composition data for model input
#'
#' Transforms a pivoted wide-format length-frequency data frame (one row per fishery
#' x timestep, bins as columns) into the set of model-ready arrays and vectors
#' that \code{opal_model} and \code{get_length_like} expect.  The function
#' validates bin alignment against the length structure in \code{data}, optionally
#' filters to a subset of fisheries, and appends every derived object directly to
#' \code{data}.
#'
#' @param data a \code{list} containing at minimum the length structure scalars
#'   \code{len_bin_start}, \code{len_bin_width}, \code{n_len}, and
#'   \code{n_fishery}.
#' @param lf_wide a \code{data.frame} in wide format: columns \code{fishery},
#'   \code{year}, \code{month}, \code{ts} plus one numeric column per length
#'   bin (named by the bin lower-bound value).  Typically produced by
#'   \code{tidyr::pivot_wider()} on the long-format \code{wcpo_bet_lf} object.
#' @param lf_keep_fisheries integer vector of fishery indices to retain.  Pass
#'   \code{NULL} (the default) to keep all fisheries present in \code{lf_wide}.
#' @param lf_switch integer likelihood type selector passed straight through to
#'   \code{data$lf_switch} (1 = multinomial, 2 = Dirichlet, 3
#'   = Dirichlet-multinomial).  Default \code{1L}.
#' @param lf_minbin integer vector length \code{n_fishery} giving the minimum
#'   bin index (1-based) used for each fishery.  Defaults to
#'   \code{rep(1L, n_fishery)}.
#' @param lf_maxbin integer vector length \code{n_fishery} giving the maximum
#'   bin index (1-based) used for each fishery.  Defaults to
#'   \code{rep(n_len, n_fishery)}.
#'
#' @return The input \code{data} list with the following elements appended or
#'   updated:
#'   \describe{
#'     \item{\code{len_lower}, \code{len_upper}, \code{len_mid}}{Length bin
#'       boundary and midpoint vectors derived from the scalar inputs.}
#'     \item{\code{lf_switch}}{Passed through from the argument.}
#'     \item{\code{lf_obs_in}}{Matrix of observed proportions,
#'       \code{[n_obs x n_len]}.}
#'     \item{\code{lf_n}}{Numeric vector of sample sizes per observation row.}
#'     \item{\code{lf_fishery}}{Integer vector of fishery index per observation
#'       row.}
#'     \item{\code{lf_fishery_f}}{Integer vector of unique fishery indices with
#'       LF data.}
#'     \item{\code{lf_n_f}}{Integer vector of observation counts per fishery.}
#'     \item{\code{lf_year}}{Integer vector of model timestep index (1-based)
#'       per observation row.}
#'     \item{\code{lf_season}}{Integer vector of season index (all 1).}
#'     \item{\code{lf_minbin}, \code{lf_maxbin}}{Passed through from arguments.}
#'     \item{\code{removal_switch_f}}{Integer vector \code{[n_fishery]} of
#'       removal flags (all 0).}
#'     \item{\code{n_lf}}{Total number of observation rows.}
#'     \item{\code{lf_obs_data}}{Named list \code{list(obs = lf_obs_list)} for
#'       plain-R access in \code{get_length_like}.}
#'     \item{\code{lf_obs_flat}}{Flattened numeric vector of counts (for
#'       multinomial, \code{lf_switch = 1}).}
#'     \item{\code{lf_obs_ints}}{Flattened integer vector of rounded counts
#'       (for Dirichlet-multinomial, \code{lf_switch = 3}).}
#'     \item{\code{lf_obs_prop}}{Flattened numeric vector of normalised
#'       proportions (for Dirichlet, \code{lf_switch = 2}).}
#'     \item{\code{lf_nbins}}{Number of bins used in each observation (scalar,
#'       derived from the first fishery's min/max bin setting).}
#'   }
#'
#' @details
#' Rows in \code{lf_wide} with a total sample size of zero are silently
#' removed before any other processing.  Bin alignment between the data frame
#' columns and the model's length structure is checked with
#' \code{stopifnot()}.
#'
#' @seealso \code{\link{get_length_like}}, \code{\link{opal_model}}
#' @export
prep_lf_data <- function(data,
                         lf_wide,
                         lf_keep_fisheries = NULL,
                         lf_switch         = 1L,
                         lf_minbin         = NULL,
                         lf_maxbin         = NULL) {

  # ---- 1. Extract bin columns and build obs-count matrix ----
  meta_cols <- c("fishery", "year", "month", "ts")
  bin_cols  <- sort(as.numeric(setdiff(names(lf_wide), meta_cols)))
  lf_obs_counts <- as.matrix(lf_wide[, as.character(bin_cols)])

  # ---- 2. Drop rows with zero total sample size ----
  lf_n <- rowSums(lf_obs_counts)
  keep <- lf_n > 0
  lf_wide       <- lf_wide[keep, ]
  lf_obs_counts <- lf_obs_counts[keep, ]
  lf_n          <- lf_n[keep]

  # ---- 3. Convert counts to proportions ----
  lf_obs <- lf_obs_counts / lf_n

  # ---- 4. Validate bin alignment against model length structure ----
  expected_bins <- seq(data$len_bin_start, by = data$len_bin_width,
                       length.out = data$n_len)
  stopifnot("Length bins in lf_wide do not match data length structure" =
              all(bin_cols == expected_bins))

  # ---- 5. Derive length bin boundaries (single source of truth) ----
  data$len_lower <- seq(from = data$len_bin_start, by = data$len_bin_width,
                        length.out = data$n_len)
  data$len_upper <- data$len_lower + data$len_bin_width
  data$len_mid   <- data$len_lower + data$len_bin_width / 2

  # ---- 6. Optional fishery filter ----
  if (!is.null(lf_keep_fisheries)) {
    lf_use        <- lf_wide$fishery %in% lf_keep_fisheries
    lf_wide       <- lf_wide[lf_use, ]
    lf_obs_counts <- lf_obs_counts[lf_use, ]
    lf_obs        <- lf_obs[lf_use, ]
    lf_n          <- lf_n[lf_use]
  }

  # ---- 7. Default min/max bin arguments ----
  if (is.null(lf_minbin)) lf_minbin <- rep(1L, data$n_fishery)
  if (is.null(lf_maxbin)) lf_maxbin <- rep(as.integer(data$n_len), data$n_fishery)

  # ---- 8. Attach index vectors and scalars to data ----
  data$lf_switch       <- as.integer(lf_switch)
  data$lf_obs_in       <- lf_obs
  data$lf_n            <- lf_n
  data$lf_fishery      <- as.integer(lf_wide$fishery)
  data$lf_fishery_f    <- unique(data$lf_fishery)
  data$lf_n_f          <- as.integer(table(data$lf_fishery))
  data$lf_year         <- as.integer(lf_wide$ts)  # model timestep (1-based)
  data$lf_season       <- rep(1L, nrow(lf_wide))
  data$lf_minbin       <- lf_minbin
  data$lf_maxbin       <- lf_maxbin
  data$removal_switch_f <- rep(0L, data$n_fishery)
  data$n_lf            <- nrow(lf_obs)

  # ---- 9. Build lf_obs_list (list of per-fishery obs matrices) ----
  lf_fishery_f  <- data$lf_fishery_f
  lf_n_f        <- data$lf_n_f
  n_f           <- length(lf_fishery_f)
  n_len_local   <- ncol(data$lf_obs_in)
  lf_n_fi       <- split(data$lf_n, data$lf_fishery)
  lf_row_fi     <- split(seq_len(nrow(data$lf_obs_in)), data$lf_fishery)

  lf_obs_list <- vector("list", n_f)
  for (j in seq_len(n_f)) {
    f    <- lf_fishery_f[j]
    bmin <- lf_minbin[f]
    bmax <- lf_maxbin[f]
    rows <- lf_row_fi[[j]]
    m    <- matrix(0, lf_n_f[j], bmax - bmin + 1L)
    for (i in seq_len(lf_n_f[j])) {
      obs <- data$lf_obs_in[rows[i], ]
      if (bmin > 1)          obs[bmin] <- sum(obs[1:bmin])
      if (bmax < n_len_local) obs[bmax] <- sum(obs[bmax:n_len_local])
      obs     <- obs[bmin:bmax]
      m[i, ]  <- obs * lf_n_fi[[j]][i]
    }
    lf_obs_list[[j]] <- m
  }

  # ---- 10. Attach flattened observation vectors ----
  # Plain-R nested list (for get_length_like reporting)
  data$lf_obs_data <- list(obs = lf_obs_list)

  # For multinomial (lf_switch = 1): unrounded counts
  data$lf_obs_flat <- unlist(lapply(lf_obs_list,
                                    function(m) as.numeric(t(m))))

  # For Dirichlet-multinomial (lf_switch = 3): rounded integer counts
  data$lf_obs_ints <- unlist(lapply(lf_obs_list,
                                    function(m) as.numeric(t(round(m)))))

  # For Dirichlet (lf_switch = 2): row-normalised proportions
  data$lf_obs_prop <- unlist(lapply(lf_obs_list, function(m) {
    props <- t(apply(m, 1, function(row) {
      p <- row / sum(row)
      p <- p + 1e-8
      p / sum(p)
    }))
    as.numeric(t(props))
  }))

  # Number of bins per observation (scalar; same for all obs when bmin/bmax are uniform)
  data$lf_nbins <- ncol(lf_obs_list[[1]])

  # Clean up any stale fields that earlier code may have left
  data$lf_obs     <- NULL
  data$lf_obs_vec <- NULL
  data$lf_nbins_f <- NULL

  return(data)
}
