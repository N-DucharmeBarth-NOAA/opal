#' Extract Weight Composition Fit from MFCL Model
#'
#' Extracts observed and expected weight composition data from MFCL weight.fit
#' file, aggregates across time, and optionally rebins to a target bin structure.
#'
#' @param weight_fit_file Character. Path to MFCL weight.fit file
#' @param frq_file Character. Path to MFCL .frq file (for bin structure)
#' @param model_id Character. Model identifier for output
#' @param fishery_names Character vector. Optional fleet names
#' @param harmonize_bins Logical. Rebin to target structure? Default FALSE
#' @param target_bins Numeric vector. Target bin edges. Only used if harmonize_bins = TRUE
#' @param output_dir Character. Directory for output CSV
#' @param verbose Logical. Print progress messages? Default TRUE
#' 
#' @return data.table with columns: id, Fleet, Fleet_name, Used, Kind, Sex, Bin,
#'   Obs, Exp, Dev, effN, Nsamp_in, Nsamp_adj
#' 
#' @details
#' Reads MFCL weight.fit file using either FLR4MFCL::read.MFCLWgtFit() or
#' the parse_fit_file() function. Aggregates across time and optionally rebins.
#' 
#' Note: Not all MFCL models have weight composition data. Function handles
#' gracefully when weight.fit file doesn't exist.
#'
#' @examples
#' \dontrun{
#'   wt_comp = extract_mfcl_weight_comp(
#'     weight_fit_file = "model-files/mfcl/v11/weight.fit",
#'     frq_file = "model-files/mfcl/v11/bet.frq",
#'     model_id = "mfcl-v11",
#'     output_dir = "model-files/mfcl/v11",
#'     harmonize_bins = TRUE,
#'     target_bins = seq(0, 140, by = 2)
#'   )
#' }
#'
#' @export
extract_mfcl_weight_comp = function(weight_fit_file, frq_file, model_id,
                                    fishery_names = NULL,
                                    harmonize_bins = FALSE,
                                    target_bins = NULL,
                                    output_dir,
                                    verbose = TRUE) {
  
  require(data.table)
  require(magrittr)
  
  # Check if weight.fit file exists
  if (!file.exists(weight_fit_file)) {
    warning("Weight composition file not found: ", weight_fit_file)
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
  
  # Check if FLR4MFCL is available
  if (requireNamespace("FLR4MFCL", quietly = TRUE)) {
    if(verbose) cat("Reading MFCL weight.fit file using FLR4MFCL...\n")
    
    # Use FLR4MFCL to read weight.fit
    wt_fit = FLR4MFCL::read.MFCLWgtFit(weight_fit_file)
    
    # Convert to data.table format
    # The structure from FLR4MFCL needs to be converted
    # This is a placeholder - actual implementation would depend on FLR4MFCL structure
    wt_dt = as.data.table(wt_fit)
    
  } else {
    # Fallback: Try to parse manually or use parse_fit_file if available
    if(verbose) cat("FLR4MFCL not available, attempting alternative parsing...\n")
    
    # Check if parse_fit_file exists
    parse_fit_file_path = file.path(dirname(dirname(output_dir)), "code", "ss3", "helper-fns", "parse_fit_file.r")
    
    if (file.exists(parse_fit_file_path)) {
      source(parse_fit_file_path)
      frq_stem = sub("\\.frq$", "", basename(frq_file))
      wt_dt = parse_fit_file(dirname(weight_fit_file), "weight", frq_stem)
    } else {
      stop("Cannot parse weight.fit file: FLR4MFCL not available and parse_fit_file.r not found")
    }
  }
  
  # Check if we got any data
  if(is.null(wt_dt) || nrow(wt_dt) == 0) {
    warning("No weight composition data found in file: ", weight_fit_file)
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
  
  # Extract relevant columns
  # Column names depend on the parsing method used
  # Typical MFCL columns: fishery, bin_mid, prop_obs, prop_pred, nsamp
  if(verbose) cat("Processing weight composition data...\n")
  
  # Standardize column names based on what's available
  if("fishery" %in% names(wt_dt)) {
    # Already has correct column names
  } else if("Fleet" %in% names(wt_dt)) {
    setnames(wt_dt, "Fleet", "fishery", skip_absent = TRUE)
  }
  
  # Map other column names
  if("Obs" %in% names(wt_dt)) setnames(wt_dt, "Obs", "prop_obs", skip_absent = TRUE)
  if("Exp" %in% names(wt_dt)) setnames(wt_dt, "Exp", "prop_pred", skip_absent = TRUE)
  if("Bin" %in% names(wt_dt)) setnames(wt_dt, "Bin", "bin_mid", skip_absent = TRUE)
  if("Nsamp_in" %in% names(wt_dt)) setnames(wt_dt, "Nsamp_in", "nsamp", skip_absent = TRUE)
  
  # Aggregate across time (sum proportions, sum sample sizes)
  wt_agg = wt_dt[,.(Obs = sum(prop_obs, na.rm = TRUE), 
                    Exp = sum(prop_pred, na.rm = TRUE),
                    Nsamp_in = sum(nsamp, na.rm = TRUE)),
                 by = .(fishery, bin_mid)]
  
  # Calculate deviation
  wt_agg[,Dev := Obs - Exp]
  
  # Apply bin harmonization if requested
  if (harmonize_bins && !is.null(target_bins)) {
    if(verbose) cat("Harmonizing bins...\n")
    
    # Get MFCL bin structure from frq file
    parse_frq_path = file.path(dirname(dirname(output_dir)), "code", "mfcl", "helper-fns", "parse-frq.r")
    if(!file.exists(parse_frq_path)) {
      parse_frq_path = file.path(dirname(dirname(dirname(output_dir))), "code", "mfcl", "helper-fns", "parse-frq.r")
    }
    source(parse_frq_path)
    
    frq = parse_frq(frq_file)
    
    # Get weight bin structure (indices 6, 7, 8 for WFIntervals, WFFirst, WFWidth)
    bin_lower = seq(from = lf_range(frq)[7], by = lf_range(frq)[8],
                   length.out = lf_range(frq)[6])
    src_edges = c(bin_lower, max(bin_lower) + lf_range(frq)[8])
    
    # Source rebin_composition function
    rebin_path = file.path(dirname(dirname(output_dir)), "code", "ss3", "helper-fns", "rebin_composition.r")
    if(!file.exists(rebin_path)) {
      rebin_path = file.path(dirname(dirname(dirname(output_dir))), "code", "ss3", "helper-fns", "rebin_composition.r")
    }
    source(rebin_path)
    
    wt_agg = wt_agg[, {
      # Rebin
      obs_rebinned = rebin_composition(src_edges, Obs, target_bins)
      exp_rebinned = rebin_composition(src_edges, Exp, target_bins)
      
      data.table(
        Bin = target_bins[-length(target_bins)],
        Obs = obs_rebinned,
        Exp = exp_rebinned,
        Nsamp_in = Nsamp_in[1]  # Keep first value (not rebinned)
      )
    }, by = .(fishery)]
    
    wt_agg[,Dev := Obs - Exp]
  } else {
    # Rename bin_mid to Bin
    setnames(wt_agg, "bin_mid", "Bin", skip_absent = TRUE)
  }
  
  # Add standard columns
  wt_agg[,Used := 1]  # All MFCL data is used
  wt_agg[,Kind := 3]  # Weight composition type
  wt_agg[,Sex := 0]   # Aggregated
  wt_agg[,effN := Nsamp_in]
  wt_agg[,Nsamp_adj := Nsamp_in]
  
  # Add fleet names
  if (!is.null(fishery_names) && length(fishery_names) >= max(wt_agg$fishery)) {
    wt_agg[,Fleet_name := fishery_names[fishery]]
  } else {
    wt_agg[,Fleet_name := paste0("Fishery_", fishery)]
  }
  
  # Rename columns
  setnames(wt_agg, "fishery", "Fleet", skip_absent = TRUE)
  
  # Add model id
  wt_agg[,id := model_id]
  
  # Reorder columns
  wt_agg = wt_agg[,.(id, Fleet, Fleet_name, Used, Kind, Sex, Bin,
                     Obs, Exp, Dev, effN, Nsamp_in, Nsamp_adj)]
  
  # Write CSV
  if(verbose) cat("Writing comp_size.csv...\n")
  fwrite(wt_agg, file.path(output_dir, "comp_size.csv"))
  
  if(verbose) cat("Done!\n")
  
  # Return data.table
  return(wt_agg)
}
