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
#' @param wf_cap positive integer (or \code{NULL}).  When supplied, each
#'   effective sample size is capped: \eqn{n_i = \min(n_i, \text{wf\_cap})}.
#'   Applied after proportions are computed so that observed compositions are
#'   unaffected.  Default \code{NULL} (no cap).  This is a Multifan-CL legacy
#'   feature.
#' @param wf_addtocomp small non-negative numeric constant added to observed
#'   composition proportions after tail compression and before renormalisation.
#'   This robustifies zero bins. Default \code{1e-08}.
#' @param wf_upper_edge Final upper edge, or one per element of a list of wide
#'   data frames. Defaults to the native weight-grid upper edge. Long input supplies
#'   \code{bin_lower}, \code{bin_upper}, and \code{count}; extra columns such as
#'   \code{source} identify separate observations.
#' @param wf_tails Prediction tail handling: \code{"accumulate"} (default)
#'   or \code{"truncate"}.
#' @param wf_edge_tol Absolute tolerance for deduplicating scheme edges;
#'   defaults to \code{1e-8 * data$wt_bin_width}.
#' @details Input may be one wide data frame, a list of wide frames, or long
#'   data. \code{wf_scheme_edges} records observed edges, including native
#'   min/max pooling; \code{wf_rebin_s} contains observed-bin by model-length-bin
#'   prediction matrices. Length edges are converted using \code{lw_a * L^lw_b}.
#'   Native schemes preserve the legacy two-stage weight mapping and pooling.
#'   \code{wf_scheme} and \code{wf_nbins_i} follow fishery-grouped flat order;
#'   \code{wf_row_fi} indexes retained input order. Heterogeneous observations
#'   use lists where a rectangular matrix is not possible.
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
#'     \item{\code{wf_obs_ints}}{Flattened integer counts that sum to
#'       \code{wf_n_int} per observation (for Dirichlet-multinomial,
#'       \code{wf_switch = 3}).}
#'     \item{\code{wf_obs_prop}}{Flattened numeric vector of normalised
#'       proportions (for Dirichlet, \code{wf_switch = 2}).}
#'     \item{\code{wf_addtocomp}}{Stored value of the add-to-composition
#'       constant used to robustify observed composition bins.}
#'     \item{\code{wf_n}}{Numeric vector of sample sizes per observation row.}
#'     \item{\code{wf_n_int}}{Integer Dirichlet-multinomial sample sizes
#'       obtained by half-up rounding of \code{wf_n}.}
#'     \item{\code{wf_fishery}}{Integer vector of fishery index per observation
#'       row.}
#'     \item{\code{wf_fishery_f}}{Integer vector of unique fishery indices with
#'       WF data.}
#'     \item{\code{wf_n_f}}{Integer vector of observation counts per fishery.}
#'     \item{\code{wf_year}}{Integer vector of model timestep per observation
#'       row.}
#'     \item{\code{wf_year_fi}, \code{wf_n_fi}, \code{wf_n_int_fi}, \code{wf_row_fi}}{Lists split
#'       by fishery containing model timesteps, effective sample sizes, and
#'       row indices. Precomputed so \code{opal_model()} does not rebuild them
#'       on every objective evaluation.}
#'     \item{\code{wf_minbin}, \code{wf_maxbin}}{Passed through from arguments.}
#'     \item{\code{wf_var_adjust}}{Passed through from argument; numeric vector
#'       \code{[n_fishery]} of variance-adjustment divisors applied to
#'       \code{wf_n}.}
#'   }
#' @details Dirichlet-multinomial counts use half-up rounding of effective
#' sample sizes and largest-remainder allocation. Observations that round to
#' zero are skipped by the Dirichlet-multinomial likelihood.
#' @seealso \code{\link{prep_lf_data}}, \code{\link{opal_model}}
#' @export
prep_wf_data <- function(data, wf_wide, wf_keep_fisheries = NULL,
                         wf_switch = 1L, wf_minbin = NULL,
                         wf_maxbin = NULL, wf_var_adjust = NULL,
                         wf_cap = NULL,
                         wf_addtocomp = 1e-08,
                         wf_upper_edge = NULL,
                         wf_tails = c("accumulate", "truncate"),
                         wf_edge_tol = 1e-8 * data$wt_bin_width) {

  wf_tails <- match.arg(wf_tails)
  if (length(wf_edge_tol) != 1L || !is.finite(wf_edge_tol) || wf_edge_tol < 0)
    stop("wf_edge_tol must be finite and nonnegative.", call. = FALSE)
  edges <- data$wt_bin_start + (0:data$n_wt) * data$wt_bin_width
  input <- .opal_composition_input(wf_wide, wf_upper_edge, edges, wf_edge_tol)
  .opal_prepare_compositions(data, input, "wf", wf_keep_fisheries,
    wf_switch, wf_minbin, wf_maxbin, wf_var_adjust, wf_cap, wf_addtocomp,
    wf_tails, wf_edge_tol)
}
