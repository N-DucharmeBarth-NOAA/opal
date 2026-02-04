#' Extract Weight Composition Fit from SS3 Model
#'
#' Extracts observed and expected weight composition data from Stock Synthesis
#' generalized size composition output, aggregates across time, and optionally
#' rebins to a target bin structure.
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
#' Extracts weight composition data from SS3 sizedbase component (generalized
#' size composition) and aggregates across all years. If harmonize_bins = TRUE,
#' applies rebin_composition() to convert data to target bin structure.
#'
#' @examples
#' \dontrun{
#'   # Without bin harmonization
#'   wt_comp = extract_ss3_weight_comp(
#'     model_dir = "model-files/ss3/01-bet-base",
#'     model_id = "01-bet-base"
#'   )
#'   
#'   # With bin harmonization to 2kg bins
#'   target_bins = seq(0, 140, by = 2)
#'   wt_comp = extract_ss3_weight_comp(
#'     model_dir = "model-files/ss3/01-bet-base",
#'     model_id = "01-bet-base",
#'     harmonize_bins = TRUE,
#'     target_bins = target_bins
#'   )
#' }
#'
#' @export
extract_ss3_weight_comp = function(model_dir, model_id,
                                   harmonize_bins = FALSE,
                                   target_bins = NULL,
                                   verbose = TRUE) {
  
  require(r4ss)
  require(data.table)
  require(magrittr)
  
  # Check if Report.sso exists
  if(!file.exists(file.path(model_dir, "Report.sso"))) {
    stop(sprintf("Report.sso not found in %s", model_dir))
  }
  
  # Read SS3 output
  if(verbose) cat("Reading SS3 model output...\n")
  tmp_report = r4ss::SS_output(model_dir, verbose = FALSE, printstats = FALSE)
  tmp_flt = tmp_report$definitions
  
  # Extract weight composition data (generalized size comp)
  if(is.null(tmp_report$sizedbase) || nrow(tmp_report$sizedbase) == 0) {
    warning("No weight composition data found in sizedbase")
    # Return empty data.table with correct structure
    return(data.table(
      id = character(),
      Fleet = integer(),
      Fleet_name = character(),
      Used = integer(),
      Kind = integer(),
      Sex = integer(),
      Bin = numeric(),
      Obs = numeric(),
      Exp = numeric(),
      Dev = numeric(),
      effN = numeric(),
      Nsamp_in = numeric(),
      Nsamp_adj = numeric()
    ))
  }
  
  if(verbose) cat("Extracting weight composition data...\n")
  wt_dt = as.data.table(tmp_report$sizedbase) %>%
    .[,.(Yr.S, Fleet, Used, Kind, Sex, Bin, Obs, Exp, effN, Nsamp_in, Nsamp_adj)]
  
  # Aggregate across time (sum by Fleet, Used, Kind, Sex, Bin)
  if(verbose) cat("Aggregating across time...\n")
  wt_agg = wt_dt[,.(Obs = sum(Obs), Exp = sum(Exp), effN = sum(effN),
                    Nsamp_in = sum(Nsamp_in), Nsamp_adj = sum(Nsamp_adj)),
                 by = .(Fleet, Used, Kind, Sex, Bin)]
  
  # Calculate deviation
  wt_agg[,Dev := Obs - Exp]
  
  # Apply bin harmonization if requested
  if (harmonize_bins && !is.null(target_bins)) {
    if(verbose) cat("Harmonizing bins...\n")
    
    # Source rebin_composition function if not already loaded
    if(!exists("rebin_composition")) {
      rebin_path = file.path(dirname(model_dir), "..", "..", "code", "ss3", "helper-fns", "rebin_composition.r")
      if(!file.exists(rebin_path)) {
        rebin_path = file.path(dirname(dirname(dirname(model_dir))), "code", "ss3", "helper-fns", "rebin_composition.r")
      }
      if(file.exists(rebin_path)) {
        source(rebin_path)
      } else {
        stop("Cannot find rebin_composition.r")
      }
    }
    
    wt_agg = wt_agg[, {
      # Sort data by Bin to ensure proper ordering
      setorder(.SD, Bin)
      
      # Construct source bin edges from sorted bins
      bin_width = unique(diff(Bin))[1]
      src_edges = c(Bin, max(Bin) + bin_width)
      
      # Rebin observed and expected
      obs_rebinned = rebin_composition(src_edges, Obs, target_bins)
      exp_rebinned = rebin_composition(src_edges, Exp, target_bins)
      
      # Return rebinned data
      data.table(
        Bin = target_bins[-length(target_bins)],
        Obs = obs_rebinned,
        Exp = exp_rebinned,
        effN = sum(effN),          # Sum sample sizes
        Nsamp_in = sum(Nsamp_in),  # Sum sample sizes
        Nsamp_adj = sum(Nsamp_adj) # Sum sample sizes
      )
    }, by = .(Fleet, Used, Kind, Sex)]
    
    # Recalculate deviation
    wt_agg[,Dev := Obs - Exp]
  }
  
  # Add fleet names
  wt_agg[,Fleet_name := tmp_flt$Fleet_name[Fleet]]
  
  # Add model id
  wt_agg[,id := model_id]
  
  # Reorder columns
  wt_agg = wt_agg[,.(id, Fleet, Fleet_name, Used, Kind, Sex, Bin,
                     Obs, Exp, Dev, effN, Nsamp_in, Nsamp_adj)]
  
  # Write CSV
  if(verbose) cat("Writing comp_size.csv...\n")
  fwrite(wt_agg, file.path(model_dir, "comp_size.csv"))
  
  if(verbose) cat("Done!\n")
  
  # Return data.table
  return(wt_agg)
}
