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
#' @param output_dir Character. Directory for output CSV. Required if save_csv = TRUE
#' @param save_csv Logical. Save output as CSV file? Default TRUE
#' @param verbose Logical. Print progress messages? Default TRUE
#' @param aggregate Logical. Aggregate across time? Default TRUE. When FALSE,
#'   returns per-observation data with year, month, and ts columns and saves
#'   as comp_size_time.csv instead of comp_size.csv.
#' @param par_file Character. Path to MFCL .par file (e.g., 10.par). Used to
#'   read fish-flag 50 values per fishery for computing Nsamp_adj. When NULL
#'   (default), Nsamp_adj equals Nsamp_in (backward compatible).
#' @param first_yr Integer. First year of the model, passed to
#'   FLR4MFCL::read.MFCLPar(). Default 1952.
#' 
#' @return data.table with columns: id, Fleet, Fleet_name, Used, Kind, Sex, Bin,
#'   Obs, Exp, Dev, effN, Nsamp_in, Nsamp_adj. When aggregate = FALSE, also
#'   includes year, month, and ts (quarterly time step from 1952 Q1).
#' 
#' @details
#' Reads MFCL weight.fit file using FLR4MFCL::read.MFCLWgtFit() if available,
#' or can use alternative parsing methods. Aggregates across time and optionally
#' rebins to target bin structure for cross-model comparison.
#' 
#' Note: Not all MFCL models have weight composition data. Function returns
#' error if weight.fit file doesn't exist.
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
#'   
#'   # Return data.table without saving CSV
#'   wt_comp = extract_mfcl_weight_comp(
#'     weight_fit_file = "model-files/mfcl/v11/weight.fit",
#'     frq_file = "model-files/mfcl/v11/bet.frq",
#'     model_id = "mfcl-v11",
#'     save_csv = FALSE
#'   )
#' }
#'
#' @export
extract_mfcl_weight_comp = function(weight_fit_file, frq_file, model_id,
                                    fishery_names = NULL,
                                    harmonize_bins = FALSE,
                                    target_bins = NULL,
                                    output_dir = NULL,
                                    save_csv = TRUE,
                                    verbose = TRUE,
                                    zero_replace = NULL,
                                    aggregate = TRUE,
                                    par_file = NULL,
                                    first_yr = 1952L) {
  
  # Check if weight.fit file exists
  if(!file.exists(weight_fit_file)) {
    stop(sprintf("weight.fit file not found: %s", weight_fit_file))
  }
  
  # Check if frq file exists
  if(!file.exists(frq_file)) {
    stop(sprintf(".frq file not found: %s", frq_file))
  }
  
  # Check if output_dir is provided when save_csv = TRUE
  if(save_csv && is.null(output_dir)) {
    stop("output_dir must be provided when save_csv = TRUE")
  }
  
  # Step 0: Read flag 50 values from par file for Nsamp_adj calculation
  # Fish flags have flagtype == -n where n is the fleet number.
  # Nsamp_adj = min(Nsamp_in, 1000) / flag50_value_for_fleet
  flag50_lookup = NULL
  if(!is.null(par_file)) {
    if(!file.exists(par_file)) {
      stop(sprintf("par file not found: %s", par_file))
    }
    if(!requireNamespace("FLR4MFCL", quietly = TRUE)) {
      stop("Package 'FLR4MFCL' is required to read par file for flag 50 values")
    }
    if(verbose) message("Reading MFCL par file for flag 50 values: ", par_file)
    base_par = FLR4MFCL::read.MFCLPar(par_file, first.yr = first_yr)
    flag_dt = data.table::as.data.table(FLR4MFCL::flags(base_par))
    # Fish flags: flagtype is negative, abs(flagtype) = fleet number
    flag50 = flag_dt[flag == 50 & flagtype < 0]
    flag50[, fleet := abs(flagtype)]
    # Build named lookup vector: fleet number -> flag 50 value
    flag50_lookup = stats::setNames(flag50$value, flag50$fleet)
    if(verbose) message("  Flag 50 values per fishery: ",
                        paste(paste0("F", names(flag50_lookup), "=", flag50_lookup),
                              collapse = ", "))
  }
  
  # Load required packages
  if(!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required but not installed")
  }
  if(!requireNamespace("magrittr", quietly = TRUE)) {
    stop("Package 'magrittr' is required but not installed")
  }
  
  if(verbose) message("Reading MFCL weight.fit from ", weight_fit_file)
  
  # Step 1: Parse weight.fit file
  # Try FLR4MFCL first
  wt_dt = NULL
  
  if(requireNamespace("FLR4MFCL", quietly = TRUE)) {
    if(verbose) message("Using FLR4MFCL::read.MFCLWgtFit()")
    tryCatch({
      wt_fit = FLR4MFCL::read.MFCLWgtFit(weight_fit_file)
      # Use wgtfits() accessor if available to get properly structured data
      if(exists("wgtfits", where = asNamespace("FLR4MFCL"))) {
        wt_dt = data.table::as.data.table(FLR4MFCL::wgtfits(wt_fit))
      } else {
        wt_dt = data.table::as.data.table(wt_fit)
      }
    }, error = function(e) {
      warning("FLR4MFCL::read.MFCLWgtFit() failed: ", e$message)
    })
  }
  
  # If FLR4MFCL didn't work, try custom parsing
  if(is.null(wt_dt)) {
    if(verbose) message("Using custom weight.fit parser")
    
    # Read weight.fit file
    lines = readLines(weight_fit_file)
    
    # Find data lines (skip comments)
    data_lines = lines[!grepl("^#", lines) & nchar(trimws(lines)) > 0]
    
    if(length(data_lines) == 0) {
      stop("No data found in weight.fit file")
    }
    
    # Parse data (space or tab separated)
    data_list = lapply(data_lines, function(x) {
      as.numeric(unlist(strsplit(trimws(x), "\\s+")))
    })
    
    # Convert to matrix
    data_mat = do.call(rbind, data_list)
    
    # Get bin structure from frq file
    if(!exists("parse_frq")) {
      frq_parser_paths = c(
        file.path(dirname(dirname(dirname(weight_fit_file))), "code", "mfcl", "helper-fns", "parse-frq.r"),
        "code/mfcl/helper-fns/parse-frq.r"
      )
      
      frq_parser_file = NULL
      for(path in frq_parser_paths) {
        if(file.exists(path)) {
          frq_parser_file = path
          break
        }
      }
      
      if(!is.null(frq_parser_file)) {
        source(frq_parser_file)
      } else {
        stop("Could not find parse-frq.r to determine bin structure")
      }
    }
    
    if(exists("parse_frq")) {
      frq = parse_frq(frq_file)
      # Weight bin structure from indices 6, 7, 8: WFIntervals, WFFirst, WFWidth
      n_bins = lf_range(frq)[6]
      bin_first = lf_range(frq)[7]
      bin_width = lf_range(frq)[8]
      
      # Create bin structure
      bin_lower = seq(from = bin_first, by = bin_width, length.out = n_bins)
      
      n_cols = ncol(data_mat)
      
      # Typical format: fishery, year, month, week, nsamp, then weight bins
      if(n_cols >= 5 + 2 * n_bins) {
        obs_cols = 6:(5 + n_bins)
        pred_cols = (6 + n_bins):(5 + 2 * n_bins)
        
        # Reshape to long format
        wt_list = list()
        for(i in 1:nrow(data_mat)) {
          for(j in 1:n_bins) {
            wt_list[[length(wt_list) + 1]] = list(
              fishery = data_mat[i, 1],
              year = data_mat[i, 2],
              month = data_mat[i, 3],
              bin_mid = bin_lower[j] + bin_width / 2,
              bin_lower = bin_lower[j],
              prop_obs = data_mat[i, obs_cols[j]],
              prop_pred = data_mat[i, pred_cols[j]],
              nsamp = data_mat[i, 5]
            )
          }
        }
        
        wt_dt = data.table::rbindlist(wt_list)
      } else {
        stop("Unexpected weight.fit file format")
      }
    } else {
      stop("Could not find parse-frq.r to determine bin structure")
    }
  }
  
  # Step 2: Ensure required columns exist and standardize names
  if("fishery" %in% names(wt_dt)) {
    data.table::setnames(wt_dt, "fishery", "Fleet", skip_absent = TRUE)
  }
  if("weight" %in% names(wt_dt)) {
    data.table::setnames(wt_dt, "weight", "Bin", skip_absent = TRUE)
  }
  if("obs" %in% names(wt_dt)) {
    data.table::setnames(wt_dt, "obs", "Obs", skip_absent = TRUE)
  }
  if("pred" %in% names(wt_dt)) {
    data.table::setnames(wt_dt, "pred", "Exp", skip_absent = TRUE)
  }
  if("sample_size" %in% names(wt_dt)) {
    data.table::setnames(wt_dt, "sample_size", "Nsamp_in", skip_absent = TRUE)
  }
  # Handle custom parser column names
  if("bin_lower" %in% names(wt_dt)) {
    data.table::setnames(wt_dt, "bin_lower", "Bin", skip_absent = TRUE)
  }
  if("prop_obs" %in% names(wt_dt)) {
    data.table::setnames(wt_dt, "prop_obs", "Obs", skip_absent = TRUE)
  }
  if("prop_pred" %in% names(wt_dt)) {
    data.table::setnames(wt_dt, "prop_pred", "Exp", skip_absent = TRUE)
  }
  if("nsamp" %in% names(wt_dt)) {
    data.table::setnames(wt_dt, "nsamp", "Nsamp_in", skip_absent = TRUE)
  }
  
  # Step 2b: Return pre-aggregated (time-resolved) data if requested
  if(!aggregate) {
    if(harmonize_bins) {
      warning("harmonize_bins is ignored when aggregate = FALSE")
    }
    
    if(verbose) message("Returning pre-aggregated (time-resolved) weight composition data")
    
    # Ensure year and month columns exist
    if(!"year" %in% names(wt_dt)) {
      stop("year column not found in parsed weight.fit data")
    }
    if(!"month" %in% names(wt_dt)) {
      stop("month column not found in parsed weight.fit data")
    }
    
    # Ensure year and month are integer type
    wt_dt[, year := as.integer(year)]
    wt_dt[, month := as.integer(month)]
    
    # Calculate ts (quarterly time step index, 1-based starting 1952 Q1)
    wt_dt[, ts := (year - 1952L) * 4L + match(month, c(2L, 5L, 8L, 11L))]
    
    # Apply zero replacement if requested
    if(!is.null(zero_replace)) {
      if(!is.numeric(zero_replace) || length(zero_replace) != 1) {
        stop("zero_replace must be a single numeric value or NULL")
      }
      if("Obs" %in% names(wt_dt)) {
        wt_dt[is.na(Obs), Obs := 0]
        wt_dt[Obs == 0, Obs := zero_replace]
      }
    }
    
    # Ensure Bin is numeric
    wt_dt[, Bin := as.numeric(Bin)]
    
    # Calculate deviation
    wt_dt[, Dev := Obs - Exp]
    
    # Add standard columns
    wt_dt[, Used := "yes"]
    wt_dt[, Kind := "SIZE"]
    wt_dt[, Sex := 1L]
    if(!"Nsamp_in" %in% names(wt_dt)) {
      wt_dt[, Nsamp_in := NA_real_]
    }
    wt_dt[, effN := Nsamp_in]
    # Compute Nsamp_adj using flag 50 values if available
    if(!is.null(flag50_lookup)) {
      wt_dt[, Nsamp_adj := {
        f50 = flag50_lookup[as.character(Fleet)]
        pmin(Nsamp_in, 1000) / f50
      }, by = Fleet]
    } else {
      wt_dt[, Nsamp_adj := Nsamp_in]
    }
    
    # Add fleet names
    if(!is.null(fishery_names) && length(fishery_names) >= max(wt_dt$Fleet)) {
      wt_dt[, Fleet_name := fishery_names[Fleet]]
    } else {
      wt_dt[, Fleet_name := paste0("Fishery_", Fleet)]
    }
    
    # Add model id
    wt_dt[, id := model_id]
    
    # Reorder columns
    wt_dt = wt_dt[, .(id, Fleet, Fleet_name, year, month, ts, Used, Kind, Sex, Bin,
                      Obs, Exp, Dev, effN, Nsamp_in, Nsamp_adj)]
    
    # Write CSV (optional)
    if(save_csv) {
      output_file = file.path(output_dir, "comp_size_time.csv")
      if(verbose) message("Writing output to ", output_file)
      data.table::fwrite(wt_dt, output_file)
    }
    
    return(wt_dt)
  }
  
  # Step 3: Aggregate across time (convert to counts, sum, then back to proportions)
  if(verbose) message("Aggregating weight composition data across time")
  
  # Apply zero replacement if requested
  if(!is.null(zero_replace)) {
    if(!is.numeric(zero_replace) || length(zero_replace) != 1) {
      stop("zero_replace must be a single numeric value or NULL")
    }
    if("Obs" %in% names(wt_dt)) {
      wt_dt[is.na(Obs), Obs := 0]
      wt_dt[Obs == 0, Obs := zero_replace]
    }
  }
  
  # Pre-compute Nsamp_adj on wt_dt before aggregation so it can be summed
  if(!is.null(flag50_lookup)) {
    wt_dt[, Nsamp_adj := {
      f50 = flag50_lookup[as.character(Fleet)]
      pmin(Nsamp_in, 1000) / f50
    }, by = Fleet]
  } else {
    wt_dt[, Nsamp_adj := Nsamp_in]
  }
  
  wt_agg = wt_dt[, .(
    Obs = sum(Obs * Nsamp_in, na.rm = TRUE) / sum(Nsamp_in, na.rm = TRUE),
    Exp = sum(Exp * Nsamp_in, na.rm = TRUE) / sum(Nsamp_in, na.rm = TRUE),
    Nsamp_in = sum(Obs * Nsamp_in, na.rm = TRUE),
    Nsamp_adj = sum(Nsamp_adj, na.rm = TRUE)
  ), by = .(Fleet, Bin)]
  
  # Ensure Bin is numeric
  wt_agg[, Bin := as.numeric(Bin)]
  
  # Step 4: Calculate deviation
  wt_agg[, Dev := Obs - Exp]
  
  # Step 5: Apply bin harmonization if requested
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
      possible_paths = c(
        file.path(dirname(dirname(dirname(weight_fit_file))), "code", "ss3", "helper-fns", "rebin_composition.r"),
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
      stop("rebin_composition function not found")
    }
    
    # Get MFCL bin structure from frq file
    if(!exists("parse_frq")) {
      frq_parser_paths = c(
        file.path(dirname(dirname(dirname(weight_fit_file))), "code", "mfcl", "helper-fns", "parse-frq.r"),
        "code/mfcl/helper-fns/parse-frq.r"
      )
      
      frq_parser_file = NULL
      for(path in frq_parser_paths) {
        if(file.exists(path)) {
          frq_parser_file = path
          break
        }
      }
      
      if(!is.null(frq_parser_file)) {
        source(frq_parser_file)
      } else {
        stop("parse-frq.r not found. Please ensure it is in code/mfcl/helper-fns/")
      }
    }
    
    if(exists("parse_frq")) {
      frq = parse_frq(frq_file)
      # Weight bins: indices 6, 7, 8 for WFIntervals, WFFirst, WFWidth
      bin_lower = seq(from = lf_range(frq)[7], by = lf_range(frq)[8],
                     length.out = lf_range(frq)[6])
      src_edges = c(bin_lower, max(bin_lower) + lf_range(frq)[8])
      
      # Rebin by Fleet
      wt_agg = wt_agg[, {
        # Match bins to source edges
        obs_vec = Obs[match(bin_lower, Bin)]
        obs_vec[is.na(obs_vec)] = 0
        exp_vec = Exp[match(bin_lower, Bin)]
        exp_vec[is.na(exp_vec)] = 0
        
        # Rebin
        obs_rebinned = rebin_composition(src_edges, obs_vec, target_bins)
        exp_rebinned = rebin_composition(src_edges, exp_vec, target_bins)
        
        # Get sample sizes (use sum)
        Nsamp_in_val = sum(Nsamp_in, na.rm = TRUE)
        Nsamp_adj_val = sum(Nsamp_adj, na.rm = TRUE)
        
        # Return rebinned data
        data.table::data.table(
          Bin = as.numeric(target_bins[-length(target_bins)]),
          Obs = obs_rebinned,
          Exp = exp_rebinned,
          Nsamp_in = Nsamp_in_val,
          Nsamp_adj = Nsamp_adj_val
        )
      }, by = .(Fleet)]
      
      # Recalculate deviation
      wt_agg[, Dev := Obs - Exp]
    }
  }
  
  # Step 6: Add standard columns
  wt_agg[, Used := "yes"]  # All MFCL data is used (character, matching SS3)
  wt_agg[, Kind := "SIZE"]  # Weight composition type (character)
  wt_agg[, Sex := 1L]      # Aggregated (integer)
  wt_agg[, effN := Nsamp_in]
  # Nsamp_adj already computed during aggregation
  
  # Step 7: Add fleet names
  if(!is.null(fishery_names) && length(fishery_names) >= max(wt_agg$Fleet)) {
    wt_agg[, Fleet_name := fishery_names[Fleet]]
  } else {
    wt_agg[, Fleet_name := paste0("Fishery_", Fleet)]
  }
  
  # Step 8: Add model id
  wt_agg[, id := model_id]
  
  # Step 9: Reorder columns
  wt_agg = wt_agg[, .(id, Fleet, Fleet_name, Used, Kind, Sex, Bin,
                      Obs, Exp, Dev, effN, Nsamp_in, Nsamp_adj)]
  
  # Step 10: Write CSV (optional)
  if(save_csv) {
    output_file = file.path(output_dir, "comp_size.csv")
    if(verbose) message("Writing output to ", output_file)
    data.table::fwrite(wt_agg, output_file)
  }
  
  # Step 11: Return data.table
  return(wt_agg)
}
