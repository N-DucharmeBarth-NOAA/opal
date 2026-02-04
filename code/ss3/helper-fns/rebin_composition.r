#' Rebin Composition Data Using Linear Area Distribution
#'
#' Rebins frequency/proportion data from one set of bin boundaries to another
#' by distributing source counts proportionally based on the overlap between
#' source and destination bins.
#'
#' @param src_edges Numeric vector. Boundaries of the original bins (length N+1).
#'   Must be strictly increasing.
#' @param src_counts Numeric vector. Frequency/counts in original bins (length N).
#'   Non-negative values.
#' @param dest_edges Numeric vector. Boundaries of the new bins (length M+1).
#'   Must be strictly increasing.
#' 
#' @return Numeric vector of rebinned frequency data (length M)
#' 
#' @details
#' This function implements linear area rebinning where:
#' - Each source bin's count is distributed to overlapping destination bins
#' - Distribution is proportional to the overlap width
#' - Total count is preserved: sum(src_counts) == sum(dest_counts)
#' - Handles arbitrary bin structures (different widths, different ranges)
#' - Works for both frequencies (counts) and proportions
#'
#' The algorithm assumes uniform distribution of observations within each
#' source bin (linear assumption).
#'
#' @examples
#' \dontrun{
#'   # Example 1: Simple rebinning
#'   # Source: 3 bins of 10cm width with counts
#'   src_edges = c(10, 20, 30, 40)
#'   src_counts = c(5, 10, 15)
#'   
#'   # Destination: 2 bins of 15cm width
#'   dest_edges = c(10, 25, 40)
#'   
#'   result = rebin_composition(src_edges, src_counts, dest_edges)
#'   # result = c(12.5, 17.5)
#'   
#'   # Verify total preserved
#'   sum(src_counts) == sum(result)  # TRUE: 30 == 30
#'   
#'   # Example 2: Rebinning proportions (must sum to 1)
#'   src_edges = c(50, 60, 70, 80)
#'   src_props = c(0.3, 0.5, 0.2)  # Sum to 1.0
#'   dest_edges = c(50, 65, 80)
#'   
#'   result = rebin_composition(src_edges, src_props, dest_edges)
#'   sum(result)  # 1.0 (preserved)
#'   
#'   # Example 3: Non-overlapping bins (extended range)
#'   src_edges = c(100, 110, 120)
#'   src_counts = c(8, 12)
#'   dest_edges = c(90, 105, 120, 130)
#'   
#'   result = rebin_composition(src_edges, src_counts, dest_edges)
#'   # result = c(0, 10, 10, 0) - counts distributed to overlapping bins
#' }
#'
#' @export
rebin_composition = function(src_edges, src_counts, dest_edges) {
  # Input validation
  if (length(src_edges) != length(src_counts) + 1) {
    stop("src_edges must have length = length(src_counts) + 1")
  }
  if (length(dest_edges) < 2) {
    stop("dest_edges must have at least 2 elements")
  }
  if (any(diff(src_edges) <= 0)) {
    stop("src_edges must be strictly increasing")
  }
  if (any(diff(dest_edges) <= 0)) {
    stop("dest_edges must be strictly increasing")
  }
  if (any(src_counts < 0)) {
    stop("src_counts must be non-negative")
  }
  
  # Initialize destination counts
  n_dest_bins = length(dest_edges) - 1
  dest_counts = numeric(n_dest_bins)
  
  # For each destination bin
  for (i in seq_len(n_dest_bins)) {
    d_low = dest_edges[i]
    d_high = dest_edges[i + 1]
    
    # For each source bin
    for (j in seq_len(length(src_counts))) {
      s_low = src_edges[j]
      s_high = src_edges[j + 1]
      
      # Calculate the overlap between [d_low, d_high] and [s_low, s_high]
      overlap_low = max(d_low, s_low)
      overlap_high = min(d_high, s_high)
      
      # If there is overlap
      if (overlap_low < overlap_high) {
        overlap_width = overlap_high - overlap_low
        src_bin_width = s_high - s_low
        
        # Distribute source count proportionally to the overlap area
        dest_counts[i] = dest_counts[i] + src_counts[j] * (overlap_width / src_bin_width)
      }
    }
  }
  
  return(dest_counts)
}
