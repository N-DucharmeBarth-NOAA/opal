#' Linear area rebinning of frequency data
#'
#' Redistributes counts from source bins to destination bins using
#' proportional overlap. For each destination bin, the function finds
#' all source bins that overlap it and assigns a fraction of each
#' source bin's count proportional to the overlap width relative to
#' the source bin width.
#'
#' AD-safe: src_edges and dest_edges are data; only src_counts may
#' carry AD values. The operation is linear in src_counts.
#'
#' @param src_edges numeric vector (length N+1). Boundaries of source bins.
#'   Must be monotonically increasing.
#' @param src_counts numeric vector (length N). Counts in each source bin.
#'   May be an AD vector from RTMB.
#' @param dest_edges numeric vector (length M+1). Boundaries of destination bins.
#'   Must be monotonically increasing.
#' @return numeric vector (length M). Rebinned counts in destination bins.
#' @importFrom RTMB ADoverload
#' @examples
#' # Identity rebin
#' edges <- 0:5
#' counts <- c(10, 20, 15, 5, 30)
#' rebin_counts(edges, counts, edges)
#'
#' # Coarsen: merge pairs of bins
#' src_edges  <- 0:4
#' src_counts <- c(10, 20, 15, 5)
#' dest_edges <- c(0, 2, 4)
#' rebin_counts(src_edges, src_counts, dest_edges)
#' @export
rebin_counts <- function(src_edges, src_counts, dest_edges) {
  if (any(diff(src_edges) <= 0)) {
    stop("src_edges must be monotonically increasing")
  }
  if (any(diff(dest_edges) <= 0)) {
    stop("dest_edges must be monotonically increasing")
  }
  N <- length(src_edges) - 1L
  M <- length(dest_edges) - 1L
  if (length(src_counts) != N) {
    stop("length(src_counts) must equal length(src_edges) - 1")
  }

  "[<-" <- RTMB::ADoverload("[<-")

  dest_counts <- numeric(M)
  for (i in seq_len(M)) {
    d_low  <- dest_edges[i]
    d_high <- dest_edges[i + 1L]
    for (j in seq_len(N)) {
      s_low  <- src_edges[j]
      s_high <- src_edges[j + 1L]
      overlap_low  <- max(d_low, s_low)
      overlap_high <- min(d_high, s_high)
      if (overlap_low < overlap_high) {
        overlap_width <- overlap_high - overlap_low
        src_bin_width <- s_high - s_low
        dest_counts[i] <- dest_counts[i] + src_counts[j] * (overlap_width / src_bin_width)
      }
    }
  }
  dest_counts
}

#' Compute rebinning weight matrix
#'
#' Returns an M x N matrix W where `W[i,j]` is the fraction of source bin j
#' that falls within destination bin i. Applying \code{W \%*\% src_counts} gives
#' the same result as \code{rebin_counts(src_edges, src_counts, dest_edges)}.
#'
#' Precomputing W is efficient when the same bin edges are reused for
#' many observations (e.g., all weight composition observations share the
#' same length-to-weight bin mapping).
#'
#' @param src_edges numeric vector (length N+1). Source bin boundaries.
#'   Must be monotonically increasing.
#' @param dest_edges numeric vector (length M+1). Destination bin boundaries.
#'   Must be monotonically increasing.
#' @return numeric matrix (M x N) of overlap fractions.
#' @examples
#' src_edges  <- 0:4
#' dest_edges <- c(0, 2, 4)
#' W <- rebin_matrix(src_edges, dest_edges)
#' # W %*% c(10, 20, 15, 5) gives the coarsened counts
#' @export
rebin_matrix <- function(src_edges, dest_edges) {
  if (any(diff(src_edges) <= 0)) {
    stop("src_edges must be monotonically increasing")
  }
  if (any(diff(dest_edges) <= 0)) {
    stop("dest_edges must be monotonically increasing")
  }
  N <- length(src_edges) - 1L
  M <- length(dest_edges) - 1L

  W <- matrix(0.0, nrow = M, ncol = N)
  for (i in seq_len(M)) {
    d_low  <- dest_edges[i]
    d_high <- dest_edges[i + 1L]
    for (j in seq_len(N)) {
      s_low  <- src_edges[j]
      s_high <- src_edges[j + 1L]
      overlap_low  <- max(d_low, s_low)
      overlap_high <- min(d_high, s_high)
      if (overlap_low < overlap_high) {
        W[i, j] <- (overlap_high - overlap_low) / (s_high - s_low)
      }
    }
  }
  W
}
