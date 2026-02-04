#' Extract Selectivity-at-Length from MFCL Model
#'
#' Extracts selectivity curves from MFCL report file, converts from age-based
#' to length-based using the growth curve, and writes to standardized CSV format.
#'
#' @param rep_file Character. Path to MFCL plot-*.rep file
#' @param par_file Character. Path to MFCL *.par file
#' @param model_id Character. Model identifier for output
#' @param first_year Numeric. First year of model. Default 1952
#' @param fishery_names Character vector. Optional fleet names. If NULL, uses default naming
#' @param output_dir Character. Directory for output CSV
#' @param verbose Logical. Print progress messages? Default TRUE
#' 
#' @return data.table with columns: id, Fleet, Fleet_name, Yr, Sex, variable, value
#' 
#' @details
#' MFCL stores selectivity by age, which must be converted to length-based
#' selectivity by:
#' 1. Extracting selectivity-at-age using FLR4MFCL::sel()
#' 2. Extracting mean length-at-age using FLR4MFCL::mean_laa()
#' 3. Joining these to get selectivity-at-length
#' 
#' The final year is determined from the report file.
#'
#' @examples
#' \dontrun{
#'   selex_dt = extract_mfcl_selectivity(
#'     rep_file = "model-files/mfcl/v11/plot-10.par.rep",
#'     par_file = "model-files/mfcl/v11/10.par",
#'     model_id = "mfcl-v11",
#'     output_dir = "model-files/mfcl/v11"
#'   )
#' }
#'
#' @export
extract_mfcl_selectivity = function(rep_file, par_file, model_id, 
                                    first_year = 1952,
                                    fishery_names = NULL,
                                    output_dir,
                                    verbose = TRUE) {
  # Validate inputs
  if(!file.exists(rep_file)) {
    stop(sprintf("Report file not found: %s", rep_file))
  }
  if(!file.exists(par_file)) {
    stop(sprintf("Parameter file not found: %s", par_file))
  }
  
  if(verbose) {
    message(sprintf("Reading MFCL files from %s", dirname(rep_file)))
  }
  
  # Read MFCL files
  tmp.rep = FLR4MFCL::read.MFCLRep(rep_file)
  tmp.par = FLR4MFCL::read.MFCLPar(par_file, first.yr = first_year)
  
  # Extract growth curve (mean length-at-age)
  tmp_growth = as.data.table(FLR4MFCL::mean_laa(tmp.rep)) %>%
    .[, .(value)] %>%
    .[order(value)] %>%
    .[, age := 1:.N] %>%
    setnames(., "value", "len") %>%
    .[, .(age, len)]
  
  if(verbose) {
    message(sprintf("Extracted growth curve: %d age classes", nrow(tmp_growth)))
  }
  
  # Extract selectivity-at-age
  selex_dt = as.data.table(FLR4MFCL::sel(tmp.rep)) %>%
    .[, model_name := model_id] %>%
    .[, .(model_name, unit, age, value)] %>%
    setnames(., c("unit", "value"), c("fishery", "selex")) %>%
    .[, fishery := as.numeric(as.character(fishery))] %>%
    .[, age := as.numeric(as.character(age))]
  
  if(verbose) {
    message(sprintf("Extracted selectivity for %d fisheries", length(unique(selex_dt$fishery))))
  }
  
  # Merge selectivity with growth to get selectivity-at-length
  selex_dt = merge(selex_dt, tmp_growth, by = "age")
  
  # Get final year from report
  # Extract year range from report file
  tmp_years = as.data.table(FLR4MFCL::sel(tmp.rep))
  if("year" %in% names(tmp_years)) {
    final_year = max(as.numeric(as.character(tmp_years$year)), na.rm = TRUE)
  } else {
    # If year not available, use dimensions or assume end year
    final_year = first_year + length(unique(tmp_years$season)) / 4 - 1
  }
  
  # Create fleet names
  n_fisheries = length(unique(selex_dt$fishery))
  if(is.null(fishery_names)) {
    Fleet_names = paste0("Fishery_", sprintf("%02d", 1:n_fisheries))
  } else {
    if(length(fishery_names) != n_fisheries) {
      warning(sprintf("Number of fishery_names (%d) does not match number of fisheries (%d). Using default names.",
                      length(fishery_names), n_fisheries))
      Fleet_names = paste0("Fishery_", sprintf("%02d", 1:n_fisheries))
    } else {
      Fleet_names = fishery_names
    }
  }
  
  # Create fleet name lookup
  fleet_lookup = data.table(
    fishery = 1:n_fisheries,
    Fleet_name = Fleet_names
  )
  
  # Standardize column names and format
  result = selex_dt %>%
    .[, .(fishery, len, selex)] %>%
    merge(., fleet_lookup, by = "fishery") %>%
    .[, .(
      id = model_id,
      Fleet = as.integer(fishery),
      Fleet_name = Fleet_name,
      Yr = as.integer(final_year),
      Sex = 0L,  # Aggregated
      variable = len,
      value = selex
    )] %>%
    .[order(Fleet, variable)]
  
  # Write CSV output
  output_file = file.path(output_dir, "selex_l.csv")
  data.table::fwrite(result, output_file)
  
  if(verbose) {
    message(sprintf("Written %d rows to %s", nrow(result), output_file))
    message(sprintf("  Fleets: %s", paste(unique(result$Fleet), collapse = ", ")))
    message(sprintf("  Length range: %.1f - %.1f cm", 
                    min(result$variable), 
                    max(result$variable)))
  }
  
  return(result)
}
