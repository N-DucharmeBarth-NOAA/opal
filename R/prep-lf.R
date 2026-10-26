# Round composition counts while preserving the integer sample size.
.opal_integerize_composition <- function(composition, sample_size) {
  stopifnot(
    is.numeric(composition),
    all(is.finite(composition)),
    all(composition >= 0),
    length(sample_size) == 1L,
    is.finite(sample_size),
    sample_size >= 0
  )
  size <- as.integer(floor(sample_size + 0.5))
  if (size == 0L) {
    return(list(counts = integer(length(composition)), size = size))
  }
  composition_sum <- sum(composition)
  if (composition_sum == 0) {
    stop("`composition` must contain positive mass when `sample_size` rounds above zero.",
         call. = FALSE)
  }
  scaled <- composition / composition_sum * size
  counts <- as.integer(floor(scaled))
  remainder <- size - sum(counts)
  if (remainder > 0L) {
    order_index <- order(-(scaled - counts), seq_along(counts))
    counts[order_index[seq_len(remainder)]] <- counts[order_index[seq_len(remainder)]] + 1L
  }
  list(counts = counts, size = size)
}

#' Prepare length composition data for model input
#'
#' Transforms wide or long length-frequency observations (one composition per
#' fishery x timestep and optional source) into model-ready arrays and vectors
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
#' @param lf_var_adjust numeric vector of length \code{n_fishery} used to
#'   scale the effective sample size for each fishery.  For observation row
#'   \eqn{i} belonging to fishery \eqn{f}, the sample size is adjusted as
#'   \eqn{n_i / \text{lf\_var\_adjust}_f}.  Values greater than 1 downweight
#'   (reduce effective \eqn{N}) and values less than 1 upweight the fishery.
#'   Defaults to \code{rep(1, n_fishery)} (no adjustment).
#' @param lf_cap positive integer (or \code{NULL}).  When supplied, each
#'   effective sample size is capped: \eqn{n_i = \min(n_i, \text{lf\_cap})}.
#'   Applied after proportions are computed so that observed compositions are
#'   unaffected.  Default \code{NULL} (no cap).  This is a Multifan-CL legacy
#'   feature.
#' @param lf_addtocomp small non-negative numeric constant added to observed
#'   composition proportions after tail compression and before renormalisation.
#'   This robustifies zero bins. Default \code{1e-08}.
#' @param lf_upper_edge Final upper edge, or one per element of a list of wide
#'   data frames. Defaults to the model length-grid upper edge.
#'   Long input instead supplies
#'   \code{bin_lower}, \code{bin_upper}, and \code{count}; additional columns
#'   (for example \code{source}) identify separate observations.
#' @param lf_tails Prediction tail handling: \code{"accumulate"} (default)
#'   extends the outer matrix edges, or \code{"truncate"} excludes outside mass.
#' @param lf_edge_tol Absolute edge tolerance; defaults to
#'   \code{1e-8 * data$len_bin_width}. Nearby model edges snap silently.
#' @param lf_allow_misaligned Allow observed edges to split model bins, with a
#'   diagnostic warning. The default rejects all edges not matching a model
#'   boundary, including edges beyond the model range.
#'
#' @details Inputs may also be lists of wide data frames or long data frames.
#'   Observed counts remain on their original grid. \code{lf_scheme_edges} is
#'   a deduplicated list of observed edges, including native min/max pooling.
#'   \code{lf_rebin_s} contains prediction matrices (observed bins by model
#'   length bins). \code{lf_scheme} and \code{lf_nbins_i} follow fishery-grouped
#'   flat observation order, while \code{lf_row_fi} indexes retained input order.
#'   Heterogeneous observations use lists instead of rectangular matrices.
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
#'     \item{\code{lf_n_int}}{Integer Dirichlet-multinomial sample sizes
#'       obtained by half-up rounding of \code{lf_n}.}
#'     \item{\code{lf_fishery}}{Integer vector of fishery index per observation
#'       row.}
#'     \item{\code{lf_fishery_f}}{Integer vector of unique fishery indices with
#'       LF data.}
#'     \item{\code{lf_n_f}}{Integer vector of observation counts per fishery.}
#'     \item{\code{lf_year}}{Integer vector of model timestep index (1-based)
#'       per observation row.}
#'     \item{\code{lf_year_fi}, \code{lf_n_fi}, \code{lf_n_int_fi}, \code{lf_row_fi}}{Lists split
#'       by fishery containing model timesteps, effective sample sizes, and
#'       row indices. Precomputed so \code{opal_model()} does not rebuild them
#'       on every objective evaluation.}
#'     \item{\code{lf_season}}{Integer vector of season index (all 1).}
#'     \item{\code{lf_minbin}, \code{lf_maxbin}}{Passed through from arguments.}
#'     \item{\code{lf_var_adjust}}{Passed through from argument; numeric vector
#'       \code{[n_fishery]} of variance-adjustment divisors applied to
#'       \code{lf_n}.}
#'     \item{\code{removal_switch_f}}{Integer vector \code{[n_fishery]} of
#'       removal flags (all 0).}
#'     \item{\code{n_lf}}{Total number of observation rows.}
#'     \item{\code{lf_obs_data}}{Named list \code{list(obs = lf_obs_list)} for
#'       plain-R access in \code{get_length_like}.}
#'     \item{\code{lf_obs_flat}}{Flattened numeric vector of counts (for
#'       multinomial, \code{lf_switch = 1}).}
#'     \item{\code{lf_obs_ints}}{Flattened integer counts that sum to
#'       \code{lf_n_int} per observation (for Dirichlet-multinomial,
#'       \code{lf_switch = 3}).}
#'     \item{\code{lf_obs_prop}}{Flattened numeric vector of normalised
#'       proportions (for Dirichlet, \code{lf_switch = 2}).}
#'     \item{\code{lf_addtocomp}}{Stored value of the add-to-composition
#'       constant used to robustify observed composition bins.}
#'     \item{\code{lf_nbins}}{Number of bins used in each observation (scalar,
#'       derived from the first fishery's min/max bin setting).}
#'   }
#'
#' @details
#' Rows in \code{lf_wide} with a total sample size of zero are silently
#' removed before any other processing.  Bin alignment between the data frame
#' columns and the model's length structure is checked with
#'   the edge tolerance. Dirichlet-multinomial counts use half-up rounding of
#' effective sample sizes and largest-remainder allocation. Observations that
#' round to zero are skipped by the Dirichlet-multinomial likelihood.
#'
#' @seealso \code{\link{get_length_like}}, \code{\link{opal_model}}
#' @export
prep_lf_data <- function(data,
                         lf_wide,
                         lf_keep_fisheries = NULL,
                         lf_switch         = 1L,
                         lf_minbin         = NULL,
                         lf_maxbin         = NULL,
                         lf_var_adjust     = NULL,
                         lf_cap            = NULL,
                         lf_addtocomp      = 1e-08,
                         lf_upper_edge     = NULL,
                         lf_tails          = c("accumulate", "truncate"),
                         lf_edge_tol       = 1e-8 * data$len_bin_width,
                         lf_allow_misaligned = FALSE) {

  lf_tails <- match.arg(lf_tails)
  if (length(lf_edge_tol) != 1L || !is.finite(lf_edge_tol) || lf_edge_tol < 0)
    stop("lf_edge_tol must be finite and nonnegative.", call. = FALSE)
  if (!is.logical(lf_allow_misaligned) || length(lf_allow_misaligned) != 1L ||
      is.na(lf_allow_misaligned))
    stop("lf_allow_misaligned must be TRUE or FALSE.", call. = FALSE)
  edges <- data$len_bin_start + (0:data$n_len) * data$len_bin_width
  input <- .opal_composition_input(lf_wide, lf_upper_edge, edges, lf_edge_tol)
  .opal_prepare_compositions(data, input, "lf", lf_keep_fisheries,
    lf_switch, lf_minbin, lf_maxbin, lf_var_adjust, lf_cap, lf_addtocomp,
    lf_tails, lf_edge_tol, lf_allow_misaligned)
}
