#' Extract Length Composition Fit from MFCL Model
#'
#' Extracts observed and expected length composition data from MFCL length.fit
#' file, aggregates across time, and optionally rebins to a target bin structure.
#'
#' @param length_fit_file Character. Path to MFCL length.fit file
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
#'   as comp_len_time.csv instead of comp_len.csv.
#' 
#' @return data.table with columns: id, Fleet, Fleet_name, Used, Kind, Sex, Bin,
#'   Obs, Exp, Dev, effN, Nsamp_in, Nsamp_adj. When aggregate = FALSE, also
#'   includes year, month, and ts (quarterly time step from 1952 Q1).
#' 
#' @details
#' Reads MFCL length.fit file using FLR4MFCL::read.MFCLLenFit() if available,
#' or can use alternative parsing methods. Aggregates across time and optionally
#' rebins to target bin structure for cross-model comparison.
#'
#' @examples
#' \dontrun{
#'   len_comp = extract_mfcl_length_comp(
#'     length_fit_file = "model-files/mfcl/v11/length.fit",
#'     frq_file = "model-files/mfcl/v11/bet.frq",
#'     model_id = "mfcl-v11",
#'     output_dir = "model-files/mfcl/v11",
#'     harmonize_bins = TRUE,
#'     target_bins = seq(10, 200, by = 5)
#'   )
#'   
#'   # Return data.table without saving CSV
#'   len_comp = extract_mfcl_length_comp(
#'     length_fit_file = "model-files/mfcl/v11/length.fit",
#'     frq_file = "model-files/mfcl/v11/bet.frq",
#'     model_id = "mfcl-v11",
#'     save_csv = FALSE
#'   )
#' }
#'
#' @export
extract_mfcl_length_comp = function(length_fit_file, frq_file, model_id,
                                    fishery_names = NULL,
                                    harmonize_bins = FALSE,
                                    target_bins = NULL,
                                    output_dir = NULL,
                                    save_csv = TRUE,
                                    verbose = TRUE,
                                    zero_replace = NULL,
                                    aggregate = TRUE) {
  
  # Check if length.fit file exists
  if(!file.exists(length_fit_file)) {
    stop(sprintf("length.fit file not found: %s", length_fit_file))
  }
  
  # Check if frq file exists
  if(!file.exists(frq_file)) {
    stop(sprintf(".frq file not found: %s", frq_file))
  }
  
  # Check if output_dir is provided when save_csv = TRUE
  if(save_csv && is.null(output_dir)) {
    stop("output_dir must be provided when save_csv = TRUE")
  }
  
  # Load required packages
  if(!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required but not installed")
  }
  if(!requireNamespace("magrittr", quietly = TRUE)) {
    stop("Package 'magrittr' is required but not installed")
  }
  
  if(verbose) message("Reading MFCL length.fit from ", length_fit_file)
  
  # Step 1: Parse length.fit file
  # Try FLR4MFCL first, fall back to custom parsing
  len_dt = NULL
  
  if(requireNamespace("FLR4MFCL", quietly = TRUE)) {
    if(verbose) message("Using FLR4MFCL::read.MFCLLenFit() and lenfits() accessor")
    tryCatch({
      len_fit = FLR4MFCL::read.MFCLLenFit(length_fit_file)
      # Use lenfits() accessor to get properly structured data
      len_dt = data.table::as.data.table(FLR4MFCL::lenfits(len_fit))
    }, error = function(e) {
      warning("FLR4MFCL::read.MFCLLenFit() or lenfits() failed: ", e$message)
    })
  }
  
  # If FLR4MFCL didn't work or isn't available, use simple parsing
  if(is.null(len_dt)) {
    if(verbose) message("Using custom length.fit parser")
    
    # Read length.fit file - simple format parser
    lines = readLines(length_fit_file)
    
    # Parse header to get structure
    # Typical format:
    # # fishery year month week ...bins...
    # data rows
    
    # Find data lines (skip comments)
    data_lines = lines[!grepl("^#", lines) & nchar(trimws(lines)) > 0]
    
    if(length(data_lines) == 0) {
      stop("No data found in length.fit file")
    }
    
    # Parse data (assuming space or tab separated)
    data_list = lapply(data_lines, function(x) {
      as.numeric(unlist(strsplit(trimws(x), "\\s+")))
    })
    
    # Convert to matrix
    data_mat = do.call(rbind, data_list)
    
    # Determine structure from frq file
    if(!exists("parse_frq")) {
      frq_parser_paths = c(
        file.path(dirname(dirname(dirname(length_fit_file))), "code", "mfcl", "helper-fns", "parse-frq.r"),
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
      source(frq_parser_file)
      frq = parse_frq(frq_file)
      n_bins = lf_range(frq)[2]  # LFIntervals
      bin_first = lf_range(frq)[3]  # LFFirst
      bin_width = lf_range(frq)[4]  # LFWidth
      
      # Create bin structure
      bin_lower = seq(from = bin_first, by = bin_width, length.out = n_bins)
      
      # Assuming columns: fishery, year, month, week, nsamp, then bins (obs and pred alternating or separate)
      # This is a simplified parser - actual format may vary
      n_cols = ncol(data_mat)
      
      # Typical format has: fishery, year, month, week, nsamp, then length bins
      # Or: fishery, year, then obs bins, then pred bins
      
      # Simple assumption: first 4 cols are fishery/year/month/week, 5th is nsamp
      # Next n_bins are observed proportions, next n_bins are predicted
      if(n_cols >= 5 + 2 * n_bins) {
        obs_cols = 6:(5 + n_bins)
        pred_cols = (6 + n_bins):(5 + 2 * n_bins)
        
        # Reshape to long format
        len_list = list()
        for(i in 1:nrow(data_mat)) {
          for(j in 1:n_bins) {
            len_list[[length(len_list) + 1]] = list(
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
        
        len_dt = data.table::rbindlist(len_list)
      } else {
        stop("Unexpected length.fit file format")
      }
    } else {
      stop("Could not find parse-frq.r to determine bin structure")
    }
  }
  
  # Step 2: Ensure required columns exist
  # Map MFCL column names to standard names
  if("fishery" %in% names(len_dt)) {
    data.table::setnames(len_dt, "fishery", "Fleet", skip_absent = TRUE)
  }
  if("length" %in% names(len_dt)) {
    data.table::setnames(len_dt, "length", "Bin", skip_absent = TRUE)
  }
  if("obs" %in% names(len_dt)) {
    data.table::setnames(len_dt, "obs", "Obs", skip_absent = TRUE)
  }
  if("pred" %in% names(len_dt)) {
    data.table::setnames(len_dt, "pred", "Exp", skip_absent = TRUE)
  }
  if("sample_size" %in% names(len_dt)) {
    data.table::setnames(len_dt, "sample_size", "Nsamp_in", skip_absent = TRUE)
  }
  # Also handle the old column names for backward compatibility with custom parser
  if("bin_lower" %in% names(len_dt)) {
    data.table::setnames(len_dt, "bin_lower", "Bin", skip_absent = TRUE)
  }
  if("prop_obs" %in% names(len_dt)) {
    data.table::setnames(len_dt, "prop_obs", "Obs", skip_absent = TRUE)
  }
  if("prop_pred" %in% names(len_dt)) {
    data.table::setnames(len_dt, "prop_pred", "Exp", skip_absent = TRUE)
  }
  if("nsamp" %in% names(len_dt)) {
    data.table::setnames(len_dt, "nsamp", "Nsamp_in", skip_absent = TRUE)
  }
  
  # Step 2b: Return pre-aggregated (time-resolved) data if requested
  if(!aggregate) {
    if(harmonize_bins) {
      warning("harmonize_bins is ignored when aggregate = FALSE")
    }
    
    if(verbose) message("Returning pre-aggregated (time-resolved) length composition data")
    
    # Ensure year and month columns exist
    if(!"year" %in% names(len_dt)) {
      stop("year column not found in parsed length.fit data")
    }
    if(!"month" %in% names(len_dt)) {
      stop("month column not found in parsed length.fit data")
    }
    
    # Ensure year and month are integer type
    len_dt[, year := as.integer(year)]
    len_dt[, month := as.integer(month)]
    
    # Calculate ts (quarterly time step index, 1-based starting 1952 Q1)
    len_dt[, ts := (year - 1952L) * 4L + match(month, c(2L, 5L, 8L, 11L))]
    
    # Apply zero replacement if requested
    if(!is.null(zero_replace)) {
      if(!is.numeric(zero_replace) || length(zero_replace) != 1) {
        stop("zero_replace must be a single numeric value or NULL")
      }
      if("Obs" %in% names(len_dt)) {
        len_dt[is.na(Obs), Obs := 0]
        len_dt[Obs == 0, Obs := zero_replace]
      }
    }
    
    # Ensure Bin is numeric
    len_dt[, Bin := as.numeric(Bin)]
    
    # Calculate deviation
    len_dt[, Dev := Obs - Exp]
    
    # Add standard columns
    len_dt[, Used := "yes"]
    len_dt[, Kind := "LEN"]
    len_dt[, Sex := 1L]
    if(!"Nsamp_in" %in% names(len_dt)) {
      len_dt[, Nsamp_in := NA_real_]
    }
    len_dt[, effN := Nsamp_in]
    len_dt[, Nsamp_adj := Nsamp_in]
    
    # Add fleet names
    if(!is.null(fishery_names) && length(fishery_names) >= max(len_dt$Fleet)) {
      len_dt[, Fleet_name := fishery_names[Fleet]]
    } else {
      len_dt[, Fleet_name := paste0("Fishery_", Fleet)]
    }
    
    # Add model id
    len_dt[, id := model_id]
    
    # Reorder columns
    len_dt = len_dt[, .(id, Fleet, Fleet_name, year, month, ts, Used, Kind, Sex, Bin,
                        Obs, Exp, Dev, effN, Nsamp_in, Nsamp_adj)]
    
    # Write CSV (optional)
    if(save_csv) {
      output_file = file.path(output_dir, "comp_len_time.csv")
      if(verbose) message("Writing output to ", output_file)
      data.table::fwrite(len_dt, output_file)
    }
    
    return(len_dt)
  }
  
  # Step 3: Aggregate across time (convert to counts, sum, then back to proportions)
  if(verbose) message("Aggregating length composition data across time")

  # Optionally replace zero observed proportions with a small positive value
  # This helps match SS3 behaviour where tiny non-zero fill-values are used
  # for empty bins and prevents systematic aggregation differences.
  if(!is.null(zero_replace)) {
    if(!is.numeric(zero_replace) || length(zero_replace) != 1) {
      stop("zero_replace must be a single numeric value or NULL")
    }
    if("Obs" %in% names(len_dt)) {
      len_dt[is.na(Obs), Obs := 0]
      len_dt[Obs == 0, Obs := zero_replace]
    }
  }
  
  len_agg = len_dt[, .(
    Obs = sum(Obs * Nsamp_in, na.rm = TRUE) / sum(Nsamp_in, na.rm = TRUE),
    Exp = sum(Exp * Nsamp_in, na.rm = TRUE) / sum(Nsamp_in, na.rm = TRUE),
    Nsamp_in = sum(Obs * Nsamp_in, na.rm = TRUE)
  ), by = .(Fleet, Bin)]
  
  # Ensure Bin is numeric (consistent with SS3)
  len_agg[, Bin := as.numeric(Bin)]
  
  # Step 4: Calculate deviation
  len_agg[, Dev := Obs - Exp]
  
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
        file.path(dirname(dirname(dirname(length_fit_file))), "code", "ss3", "helper-fns", "rebin_composition.r"),
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
        file.path(dirname(dirname(dirname(length_fit_file))), "code", "mfcl", "helper-fns", "parse-frq.r"),
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
      source(frq_parser_file)
      frq = parse_frq(frq_file)
      bin_lower = seq(from = lf_range(frq)[3], by = lf_range(frq)[4],
                     length.out = lf_range(frq)[2])
      src_edges = c(bin_lower, max(bin_lower) + lf_range(frq)[4])
      
      # Rebin by Fleet
      len_agg = len_agg[, {
        # Match bins to source edges
        obs_vec = Obs[match(bin_lower, Bin)]
        obs_vec[is.na(obs_vec)] = 0
        exp_vec = Exp[match(bin_lower, Bin)]
        exp_vec[is.na(exp_vec)] = 0
        
        # Rebin
        obs_rebinned = rebin_composition(src_edges, obs_vec, target_bins)
        exp_rebinned = rebin_composition(src_edges, exp_vec, target_bins)
        
        # Get sample size (use sum)
        Nsamp_in_val = sum(Nsamp_in, na.rm = TRUE)
        
        # Return rebinned data
        data.table::data.table(
          Bin = as.numeric(target_bins[-length(target_bins)]),
          Obs = obs_rebinned,
          Exp = exp_rebinned,
          Nsamp_in = Nsamp_in_val
        )
      }, by = .(Fleet)]
      
      # Recalculate deviation
      len_agg[, Dev := Obs - Exp]
    }
  }
  
  # Step 6: Add standard columns
  len_agg[, Used := "yes"]  # All MFCL data is used (character, matching SS3)
  len_agg[, Kind := "LEN"]  # Length composition type (character, matching SS3)
  len_agg[, Sex := 1L]      # Aggregated (MFCL typically doesn't separate) (integer)
  len_agg[, effN := Nsamp_in]      # Use input sample size as effective
  len_agg[, Nsamp_adj := Nsamp_in]
  
  # Step 7: Add fleet names
  if(!is.null(fishery_names) && length(fishery_names) >= max(len_agg$Fleet)) {
    len_agg[, Fleet_name := fishery_names[Fleet]]
  } else {
    len_agg[, Fleet_name := paste0("Fishery_", Fleet)]
  }
  
  # Step 8: Add model id
  len_agg[, id := model_id]
  
  # Step 9: Reorder columns
  len_agg = len_agg[, .(id, Fleet, Fleet_name, Used, Kind, Sex, Bin,
                        Obs, Exp, Dev, effN, Nsamp_in, Nsamp_adj)]
  
  # Step 10: Write CSV (optional)
  if(save_csv) {
    output_file = file.path(output_dir, "comp_len.csv")
    if(verbose) message("Writing output to ", output_file)
    data.table::fwrite(len_agg, output_file)
  }
  
  # Step 11: Return data.table
  return(len_agg)
}
