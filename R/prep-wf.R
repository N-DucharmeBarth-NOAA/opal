#' Prepare weight composition data for model input
#'
#' Transforms a pivoted wide-format weight-frequency data frame into
#' model-ready arrays and vectors that get_weight_like() expects.
#' Also precomputes the rebinning matrix for converting predicted
#' length compositions to weight compositions.
#'
#' @param data list containing at minimum: lw_a, lw_b, len_bin_start,
#'   len_bin_width, n_len, n_fishery, and the weight bin scalars
#'   wt_bin_start, wt_bin_width, n_wt.
#' @param wf_wide data.frame in wide format: columns fishery, year, month,
#'   ts plus one numeric column per weight bin (named by bin value).
#' @param wf_keep_fisheries integer vector of fishery indices to retain
#'   (NULL = keep all).
#' @param wf_switch integer likelihood type (1=multinomial, 2=Dirichlet,
#'   3=Dirichlet-multinomial). Default 1L.
#' @param wf_minbin integer vector (length n_fishery) min weight bin index.
#'   Defaults to rep(1L, n_fishery).
#' @param wf_maxbin integer vector (length n_fishery) max weight bin index.
#'   Defaults to rep(n_wt, n_fishery).
#' @param wf_var_adjust numeric vector (length n_fishery) variance adjustment
#'   divisors. Defaults to rep(1, n_fishery).
#' @return data list with the following weight composition elements appended:
#'   \describe{
#'     \item{\code{wf_switch}}{Passed through from the argument.}
#'     \item{\code{wt_lower}, \code{wt_upper}, \code{wt_mid}}{Weight bin
#'       boundary and midpoint vectors derived from the scalar inputs.}
#'     \item{\code{wt_bin_edges}}{Weight bin boundary vector (length n_wt + 1).}
#'     \item{\code{wf_rebin_matrix}}{Precomputed rebinning matrix (n_wt x n_len)
#'       for converting predicted length compositions to weight compositions.}
#'     \item{\code{n_wf}}{Total number of WF observation rows.}
#'     \item{\code{wf_obs_in}}{Matrix of observed proportions (n_wf x n_wt).}
#'     \item{\code{wf_obs_flat}}{Flattened numeric vector of counts (for
#'       multinomial, \code{wf_switch = 1}).}
#'     \item{\code{wf_obs_ints}}{Flattened integer vector of rounded counts
#'       (for Dirichlet-multinomial, \code{wf_switch = 3}).}
#'     \item{\code{wf_obs_prop}}{Flattened numeric vector of normalised
#'       proportions (for Dirichlet, \code{wf_switch = 2}).}
#'     \item{\code{wf_n}}{Numeric vector of sample sizes per observation row.}
#'     \item{\code{wf_fishery}}{Integer vector of fishery index per observation
#'       row.}
#'     \item{\code{wf_fishery_f}}{Integer vector of unique fishery indices with
#'       WF data.}
#'     \item{\code{wf_n_f}}{Integer vector of observation counts per fishery.}
#'     \item{\code{wf_year}}{Integer vector of model timestep per observation
#'       row.}
#'     \item{\code{wf_minbin}, \code{wf_maxbin}}{Passed through from arguments.}
#'     \item{\code{wf_var_adjust}}{Passed through from argument; numeric vector
#'       \code{[n_fishery]} of variance-adjustment divisors applied to
#'       \code{wf_n}.}
#'   }
#' @seealso \code{\link{prep_lf_data}}, \code{\link{opal_model}}
#' @export
prep_wf_data <- function(data, wf_wide, wf_keep_fisheries = NULL,
                         wf_switch = 1L, wf_minbin = NULL,
                         wf_maxbin = NULL, wf_var_adjust = NULL) {

  # ---- 1. Extract bin columns and build obs-count matrix ----
  meta_cols <- c("fishery", "year", "month", "ts")
  bin_cols  <- sort(as.numeric(setdiff(names(wf_wide), meta_cols)))
  wf_obs_counts <- as.matrix(wf_wide[, as.character(bin_cols)])

  # ---- 2. Drop rows with zero total sample size ----
  wf_n <- rowSums(wf_obs_counts)
  keep <- wf_n > 0
  wf_wide       <- wf_wide[keep, ]
  wf_obs_counts <- wf_obs_counts[keep, ]
  wf_n          <- wf_n[keep]

  # ---- 3. Convert counts to proportions ----
  wf_obs <- wf_obs_counts / wf_n

  # ---- 4. Validate bin alignment against model weight structure ----
  expected_bins <- seq(data$wt_bin_start, by = data$wt_bin_width,
                       length.out = data$n_wt)
  stopifnot("Weight bins in wf_wide do not match data weight structure" =
              all(bin_cols == expected_bins))

  # ---- 5. Derive weight bin boundaries (single source of truth) ----
  wt_lower <- seq(from = data$wt_bin_start, by = data$wt_bin_width,
                  length.out = data$n_wt)
  wt_upper <- wt_lower + data$wt_bin_width
  wt_mid   <- wt_lower + data$wt_bin_width / 2
  wt_bin_edges <- c(wt_lower, wt_upper[data$n_wt])

  data$wt_lower    <- wt_lower
  data$wt_upper    <- wt_upper
  data$wt_mid      <- wt_mid
  data$wt_bin_edges <- wt_bin_edges

  # ---- 6. Compute the rebinning matrix (length -> weight) ----
  len_lower <- seq(from = data$len_bin_start, by = data$len_bin_width,
                   length.out = data$n_len)
  len_upper <- len_lower + data$len_bin_width
  len_edges <- c(len_lower, len_upper[data$n_len])

  # Convert length edges to weight edges via L-W relationship
  wt_edges_from_len <- data$lw_a * len_edges^data$lw_b

  # W maps predicted counts at length bins to weight bins [n_wt x n_len]
  data$wf_rebin_matrix <- rebin_matrix(wt_edges_from_len, wt_bin_edges)

  # ---- 7. Optional fishery filter ----
  if (!is.null(wf_keep_fisheries)) {
    wf_use        <- wf_wide$fishery %in% wf_keep_fisheries
    wf_wide       <- wf_wide[wf_use, ]
    wf_obs_counts <- wf_obs_counts[wf_use, ]
    wf_obs        <- wf_obs[wf_use, ]
    wf_n          <- wf_n[wf_use]
  }

  # ---- 8. Default min/max bin arguments and variance adjustment ----
  if (is.null(wf_minbin))     wf_minbin     <- rep(1L, data$n_fishery)
  if (is.null(wf_maxbin))     wf_maxbin     <- rep(as.integer(data$n_wt), data$n_fishery)
  if (is.null(wf_var_adjust)) wf_var_adjust <- rep(1, data$n_fishery)

  # Apply per-fishery variance adjustment to effective sample sizes
  wf_n <- wf_n / wf_var_adjust[wf_wide$fishery]

  # ---- 9. Attach index vectors and scalars to data ----
  data$wf_switch    <- as.integer(wf_switch)
  data$wf_obs_in    <- wf_obs
  data$wf_n         <- wf_n
  data$wf_fishery   <- as.integer(wf_wide$fishery)
  data$wf_fishery_f <- unique(data$wf_fishery)
  data$wf_n_f       <- as.integer(table(data$wf_fishery))
  data$wf_year      <- as.integer(wf_wide$ts)
  data$wf_minbin    <- wf_minbin
  data$wf_maxbin    <- wf_maxbin
  data$wf_var_adjust <- wf_var_adjust
  data$n_wf         <- nrow(wf_obs)

  # ---- 10. Build wf_obs_list (list of per-fishery obs matrices) ----
  wf_fishery_f  <- data$wf_fishery_f
  wf_n_f        <- data$wf_n_f
  n_f           <- length(wf_fishery_f)
  n_wt_local    <- ncol(data$wf_obs_in)
  wf_n_fi       <- split(data$wf_n, data$wf_fishery)
  wf_row_fi     <- split(seq_len(nrow(data$wf_obs_in)), data$wf_fishery)

  wf_obs_list <- vector("list", n_f)
  for (j in seq_len(n_f)) {
    f    <- wf_fishery_f[j]
    bmin <- wf_minbin[f]
    bmax <- wf_maxbin[f]
    rows <- wf_row_fi[[j]]
    m    <- matrix(0, wf_n_f[j], bmax - bmin + 1L)
    for (i in seq_len(wf_n_f[j])) {
      obs <- data$wf_obs_in[rows[i], ]
      if (bmin > 1)           obs[bmin] <- sum(obs[1:bmin])
      if (bmax < n_wt_local)  obs[bmax] <- sum(obs[bmax:n_wt_local])
      obs     <- obs[bmin:bmax]
      m[i, ]  <- obs * wf_n_fi[[j]][i]
    }
    wf_obs_list[[j]] <- m
  }

  # ---- 11. Attach flattened observation vectors ----
  # For multinomial (wf_switch = 1): unrounded counts
  data$wf_obs_flat <- unlist(lapply(wf_obs_list,
                                    function(m) as.numeric(t(m))))

  # For Dirichlet-multinomial (wf_switch = 3): rounded integer counts
  data$wf_obs_ints <- unlist(lapply(wf_obs_list,
                                    function(m) as.integer(t(round(m)))))

  # For Dirichlet (wf_switch = 2): row-normalised proportions
  data$wf_obs_prop <- unlist(lapply(wf_obs_list, function(m) {
    props <- t(apply(m, 1, function(row) {
      p <- row / sum(row)
      p <- p + 1e-8
      p / sum(p)
    }))
    as.numeric(t(props))
  }))

  return(data)
}
