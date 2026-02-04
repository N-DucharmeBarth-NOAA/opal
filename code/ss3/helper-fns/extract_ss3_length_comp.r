#' Extract Length Composition Fit from SS3 Model
#'
#' Extracts observed and expected length composition data from Stock Synthesis
#' output, aggregates across time, and optionally rebins to a target bin structure.
#'
#' @param model_dir Character. Path to SS3 model directory containing Report.sso
#' @param model_id Character. Model identifier for output
#' @param harmonize_bins Logical. Rebin to target structure? Default FALSE
#' @param target_bins Numeric vector. Target bin edges (lower bounds + upper edge).
#'   Only used if harmonize_bins = TRUE
#' @param verbose Logical. Print progress messages? Default TRUE
#' 
#' @return data.table with columns: id, Fleet, Fleet_name, Used, Kind, Sex, Bin,
#'   Obs, Exp, Dev, effN, Nsamp_in, Nsamp_adj
#' 
#' @details
#' Extracts length composition data from SS3 lendbase component and aggregates
#' across all years. If harmonize_bins = TRUE, applies rebin_composition() to
#' convert data to target bin structure for cross-model comparison.
#'
#' @examples
#' \dontrun{
#'   # Without bin harmonization
#'   len_comp = extract_ss3_length_comp(
#'     model_dir = "model-files/ss3/01-bet-base",
#'     model_id = "01-bet-base"
#'   )
#'   
#'   # With bin harmonization to 5cm bins
#'   target_bins = seq(10, 200, by = 5)
#'   len_comp = extract_ss3_length_comp(
#'     model_dir = "model-files/ss3/01-bet-base",
#'     model_id = "01-bet-base",
#'     harmonize_bins = TRUE,
#'     target_bins = target_bins
#'   )
#' }
#'
#' @export
extract_ss3_length_comp = function(model_dir, model_id,
                                   harmonize_bins = FALSE,
                                   target_bins = NULL,
                                   verbose = TRUE) {
  
  # Check if Report.sso exists
  if(!file.exists(file.path(model_dir, "Report.sso"))) {
    stop(sprintf("Report.sso not found in %s", model_dir))
  }
  
  # Load required packages
  if(!requireNamespace("r4ss", quietly = TRUE)) {
    stop("Package 'r4ss' is required but not installed")
  }
  if(!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required but not installed")
  }
  if(!requireNamespace("magrittr", quietly = TRUE)) {
    stop("Package 'magrittr' is required but not installed")
  }
  
  if(verbose) message("Reading SS3 output from ", model_dir)
  
  # Step 1: Read SS3 output
  tmp_report = r4ss::SS_output(model_dir, verbose = FALSE, printstats = FALSE)
  tmp_flt = tmp_report$definitions
  
  # Step 2: Extract length composition data
  if(is.null(tmp_report$lendbase)) {
    stop("No length composition data found in Report.sso")
  }
  
  len_dt = data.table::as.data.table(tmp_report$lendbase)
  
  # Step 3: Select relevant columns
  required_cols = c("Fleet", "Used", "Kind", "Sex", "Bin", "Obs", "Exp", "effN", "Nsamp_in", "Nsamp_adj")
  missing_cols = setdiff(required_cols, names(len_dt))
  if(length(missing_cols) > 0) {
    stop(sprintf("Missing required columns in lendbase: %s", paste(missing_cols, collapse = ", ")))
  }
  
  len_dt = len_dt[, .SD, .SDcols = required_cols]
  
  # Step 4: Aggregate across time (sum by Fleet, Used, Kind, Sex, Bin)
  if(verbose) message("Aggregating length composition data across time")
  
  len_agg = len_dt[, .(
    Obs = sum(Obs, na.rm = TRUE),
    Exp = sum(Exp, na.rm = TRUE),
    effN = sum(effN, na.rm = TRUE),
    Nsamp_in = sum(Nsamp_in, na.rm = TRUE),
    Nsamp_adj = sum(Nsamp_adj, na.rm = TRUE)
  ), by = .(Fleet, Used, Kind, Sex, Bin)]
  
  # Step 5: Calculate deviation
  len_agg[, Dev := Obs - Exp]
  
  # Step 6: Apply bin harmonization if requested
  if(harmonize_bins) {
    if(is.null(target_bins)) {
      stop("target_bins must be provided when harmonize_bins = TRUE")
    }
    if(length(target_bins) < 2) {
      stop("target_bins must have at least 2 elements (bin edges)")
    }
    
    if(verbose) message("Applying bin harmonization")
    
    # Source rebin_composition if not already loaded
    if(!exists("rebin_composition")) {
      # Try multiple possible paths
      possible_paths = c(
        file.path(dirname(sys.frame(1)$ofile), "rebin_composition.r"),
        file.path(model_dir, "..", "..", "code", "ss3", "helper-fns", "rebin_composition.r"),
        "code/ss3/helper-fns/rebin_composition.r"
      )
      
      rebin_file = NULL
      for(path in possible_paths) {
        if(file.exists(path)) {
          rebin_file = path
          break
        }
      }
      
      if(!is.null(rebin_file)) {
        source(rebin_file)
      } else {
        stop("rebin_composition.r not found. Please ensure it is in code/ss3/helper-fns/")
      }
    }
    
    if(!exists("rebin_composition")) {
      stop("rebin_composition function not found. Please source rebin_composition.r")
    }
    
    # Rebin by Fleet, Used, Kind, Sex
    len_agg = len_agg[, {
      # Construct source bin edges
      bins_sorted = sort(unique(Bin))
      bin_width = unique(diff(bins_sorted))
      if(length(bin_width) > 1) {
        # Variable bin widths - use minimum
        bin_width = min(bin_width)
      }
      src_edges = c(bins_sorted, max(bins_sorted) + bin_width)
      
      # Match bins to data
      obs_vec = Obs[match(bins_sorted, Bin)]
      exp_vec = Exp[match(bins_sorted, Bin)]
      
      # Rebin observed and expected
      obs_rebinned = rebin_composition(src_edges, obs_vec, target_bins)
      exp_rebinned = rebin_composition(src_edges, exp_vec, target_bins)
      
      # Get sample sizes (use first value as they're aggregated)
      effN_val = effN[1]
      Nsamp_in_val = Nsamp_in[1]
      Nsamp_adj_val = Nsamp_adj[1]
      
      # Return rebinned data
      data.table::data.table(
        Bin = target_bins[-length(target_bins)],
        Obs = obs_rebinned,
        Exp = exp_rebinned,
        effN = effN_val,
        Nsamp_in = Nsamp_in_val,
        Nsamp_adj = Nsamp_adj_val
      )
    }, by = .(Fleet, Used, Kind, Sex)]
    
    # Recalculate deviation
    len_agg[, Dev := Obs - Exp]
  }
  
  # Step 7: Add fleet names
  if(!is.null(tmp_flt) && "Fleet_name" %in% names(tmp_flt)) {
    len_agg[, Fleet_name := tmp_flt$Fleet_name[Fleet]]
  } else {
    len_agg[, Fleet_name := paste0("Fleet_", Fleet)]
  }
  
  # Step 8: Add model id
  len_agg[, id := model_id]
  
  # Step 9: Reorder columns
  len_agg = len_agg[, .(id, Fleet, Fleet_name, Used, Kind, Sex, Bin,
                        Obs, Exp, Dev, effN, Nsamp_in, Nsamp_adj)]
  
  # Step 10: Write CSV
  output_file = file.path(model_dir, "comp_len.csv")
  if(verbose) message("Writing output to ", output_file)
  data.table::fwrite(len_agg, output_file)
  
  # Step 11: Return data.table
  return(len_agg)
}
