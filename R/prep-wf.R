#' Prepare weight composition data for model input
#'
#' Transforms wide or long weight-frequency observations into
#' model-ready arrays and vectors that get_weight_like() expects.
#' Also precomputes the rebinning matrix for converting predicted
#' length compositions to weight compositions.
#'
#' @param data list containing at minimum: lw_a, lw_b, len_bin_start,
#'   len_bin_width, n_len, n_fishery, and the weight bin scalars
#'   wt_bin_start, wt_bin_width, n_wt.
#' @param wf_wide A wide data frame, list of wide data frames, or long data
#'   frame. Wide frames contain \code{fishery}, \code{year}, \code{month},
#'   \code{ts}, and numeric count columns named by bin lower edges. Long input
#'   contains the same metadata plus \code{bin_lower}, \code{bin_upper}, and
#'   \code{count}. Long observations are grouped by \code{fishery}, \code{ts},
#'   and optionally \code{source}; calendar metadata is taken from the first row.
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
#'   data frames. Defaults to the largest observed bin lower edge plus
#'   \code{data$wt_bin_width}. Ignored for long input, which supplies both
#'   boundaries explicitly.
#' @param wf_tails Prediction tail handling: \code{"accumulate"} (default)
#'   extends custom scheme matrix edges to include the converted length range;
#'   \code{"truncate"} excludes outside prediction mass. Native schemes retain
#'   the legacy finite intermediate weight range; see Details.
#' @param wf_edge_tol Absolute tolerance for deduplicating scheme edges;
#'   defaults to \code{1e-8 * data$wt_bin_width}.
#' @details Input may be one wide data frame, a list of wide frames, or long
#'   data. \code{wf_scheme_edges} records observed edges, including native
#'   min/max pooling; \code{wf_rebin_s} contains observed-bin by model-length-bin
#'   prediction matrices. Length edges are converted using \code{lw_a * L^lw_b}.
#'   Native schemes preserve the legacy two-stage weight mapping and, by default,
#'   tail pooling. Even under \code{"accumulate"}, native predictions outside
#'   the full native weight grid remain excluded, preserving the old numerical
#'   mapping rather than accumulating the entire converted length range.
#'   With \code{wf_tails = "truncate"}, their prediction matrices
#'   select native min/max rows without accumulating excluded prediction mass.
#'   Prediction tail policy does not change observed native min/max count pooling.
#'   Native and custom weight operators remain distinct even when their displayed
#'   edges coincide, because their prediction tail handling can differ.
#'   Within each mapping type, schemes deduplicate observed edges within
#'   \code{wf_edge_tol}. Registry edges are never expanded for prediction tails.
#'   \code{wf_scheme} and \code{wf_nbins_i} follow fishery-grouped flat order;
#'   \code{wf_row_fi} indexes retained input order. Heterogeneous observations
#'   use lists where a rectangular matrix is not possible.
#'   The logical \code{wf_native_i}, in flat observation order, identifies
#'   observations originally on the full native grid before min/max pooling.
#'   Observed bins narrower than overlapping converted-length source bins
#'   trigger a warning, as do weight schemes entirely outside the source range.
#' @return data list with the following weight composition elements appended:
#'   \describe{
#'     \item{\code{wf_switch}}{Passed through from the argument.}
#'     \item{\code{wt_lower}, \code{wt_upper}, \code{wt_mid}}{Weight bin
#'       boundary and midpoint vectors derived from the scalar inputs.}
#'     \item{\code{wt_bin_edges}}{Weight bin boundary vector (length n_wt + 1).}
#'     \item{\code{wf_rebin_matrix}}{Precomputed rebinning matrix (n_wt x n_len)
#'       for converting predicted length compositions to weight compositions.}
#'     \item{\code{n_wf}}{Total number of WF observation rows.}
#'     \item{\code{wf_obs_in}}{Observed proportions in retained input order:
#'       a matrix with one row per observation when bin counts agree, otherwise
#'       a list of row vectors. Native input has \code{n_wt} columns.}
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
  input <- .opal_composition_input(wf_wide, wf_upper_edge, edges, wf_edge_tol,
                                   default_width = data$wt_bin_width)
  .opal_prepare_compositions(data, input, "wf", wf_keep_fisheries,
    wf_switch, wf_minbin, wf_maxbin, wf_var_adjust, wf_cap, wf_addtocomp,
    wf_tails, wf_edge_tol)
}
