#' Extract CPUE Index Fit from MFCL Model
#'
#' Extracts observed CPUE values from MFCL frequency file (frq),
#' calculates standard errors from penalty values, and optionally 
#' writes to standardized CSV format compatible with plot_index_comparison().
#'
#' @param frq_file Character. Path to MFCL .frq file
#' @param rep_file Character. Path to MFCL plot-*.par.rep file
#' @param model_id Character. Model identifier for output
#' @param fishery_names Character vector. Optional fleet names. If NULL, uses default naming
#' @param output_dir Character. Directory for output CSV. Required if save_csv = TRUE
#' @param save_csv Logical. Save output as CSV file? Default TRUE
#' @param verbose Logical. Print progress messages? Default TRUE
#' 
#' @return data.table with columns: id, Fleet, Fleet_name, Time, Obs, Exp, SE, Dev, Use
#' 
#' @details
#' MFCL catch and effort data extraction:
#' 1. Parse frq file using parse_frq() to get simpleFrq object
#' 2. Extract cateffpen data: catch, effort, penalty by fishery/year/month
#' 3. Calculate observed CPUE: cpue = catch / effort
#' 4. Extract predicted biomass from report file and normalize
#' 5. Normalize CPUE: obs = cpue / mean(cpue)
#' 6. Calculate CV from penalty: CV = 1/sqrt(2*penalty)
#' 7. Calculate SE: SE = sqrt(log(1 + CV^2))
#' 8. Filter to fisheries with penalty > 0 (CPUE data exists)
#' 
#' Time is calculated as: year + (season-1)/4 for quarterly models
#'
#' @examples
#' \dontrun{
#'   # Extract and save CSV
#'   cpue_dt = extract_mfcl_cpue(
#'     frq_file = "model-files/mfcl/v11/bet.frq",
#'     rep_file = "model-files/mfcl/v11/plot-10.par.rep",
#'     model_id = "mfcl-v11",
#'     output_dir = "model-files/mfcl/v11"
#'   )
#'   
#'   # Extract without saving CSV
#'   cpue_dt = extract_mfcl_cpue(
#'     frq_file = "model-files/mfcl/v11/bet.frq",
#'     rep_file = "model-files/mfcl/v11/plot-10.par.rep",
#'     model_id = "mfcl-v11",
#'     save_csv = FALSE
#'   )
#' }
#'
#' @export
extract_mfcl_cpue = function(frq_file, rep_file, model_id,
                             fishery_names = NULL,
                             output_dir = NULL,
                             save_csv = TRUE,
                             verbose = TRUE) {
  # Validate inputs
  if(!file.exists(frq_file)) {
    stop(sprintf("Frequency file not found: %s", frq_file))
  }
  
  if(!file.exists(rep_file)) {
    stop(sprintf("Report file not found: %s", rep_file))
  }
  
  if(save_csv && is.null(output_dir)) {
    stop("output_dir must be provided when save_csv = TRUE")
  }
  
  if(verbose) {
    message(sprintf("Parsing MFCL frequency file: %s...", basename(frq_file)))
  }
  
  # 1. Parse frq file
  base_frq = parse_frq(frq_file)
  
  if(verbose) {
    message(sprintf("Reading MFCL report file: %s...", basename(rep_file)))
  }
  
  tmp.rep = FLR4MFCL::read.MFCLRep(rep_file)
  
  if(verbose) {
    message("Extracting catch and effort data...")
  }
  
  # 2-3. Extract cateffpen data (catch, effort, penalty) and calculate CPUE
  cateffpen_dt = data.table::as.data.table(cateffpen(base_frq)) %>%
                 .[, .(cpue = catch/effort, penalty = penalty), 
                   by = .(fishery, year, month)]
  
  if(nrow(cateffpen_dt) == 0) {
    stop("No catch/effort data found")
  }
  
  if(verbose) {
    message("Filtering to fisheries with CPUE data (penalty > 0)...")
  }
  
  # 3. Filter early to fisheries with penalty > 0
  cateffpen_dt = cateffpen_dt[penalty > 0]
  
  if(nrow(cateffpen_dt) == 0) {
    stop("No CPUE data found (no fisheries with penalty > 0)")
  }
  
  # Identify fisheries with CPUE data for filtering predicted biomass
  fisheries_with_cpue = unique(cateffpen_dt$fishery)
  
  if(verbose) {
    message("Extracting predicted CPUE from report...")
  }
  
  # 4. Extract predicted biomass by unit (fishery) and aggregate to quarterly level
  pred_biomass_dt = data.table::as.data.table(FLR4MFCL::vulnBiomass(tmp.rep)) %>%
                    .[, .(unit, year, season, value)] %>%
                    .[, month := (as.numeric(season) * 3)-1] %>%
                    setnames(., "value", "pred_cpue") %>%
                    .[, unit := as.numeric(unit)] %>%
                    .[, year := as.numeric(year)] %>%
                    .[unit %in% fisheries_with_cpue]
  
  # Normalize predicted CPUE by mean for each fishery
  pred_biomass_dt = pred_biomass_dt[, pred_cpue := pred_cpue / mean(pred_cpue), by = unit]
  pred_biomass_dt = pred_biomass_dt[, .(unit, year, month, pred_cpue)] %>%
                    setnames(.,c("unit"), c("fishery"))
  
  if(verbose) {
    message("Merging observed and predicted CPUE...")
  }
  
  # 5. Merge observed and predicted by fishery/unit, year and month
  cateffpen_dt = merge(cateffpen_dt, pred_biomass_dt, 
                       by.x = c("fishery", "year", "month"),
                       by.y = c("fishery", "year", "month"),
                       all = TRUE)
  
  if(verbose) {
    message("Calculating normalized CPUE and standard errors...")
  }
  
  # 6. Normalize CPUE by mean (for each fishery separately)
  cateffpen_dt = cateffpen_dt[, obs := cpue / mean(cpue), by = fishery]
  
  # 7-8. Calculate CV and SE
  cateffpen_dt = cateffpen_dt[, cv := 1/sqrt(2*penalty)]
  cateffpen_dt = cateffpen_dt[, se_log := sqrt(log(1 + cv^2))]
  
  # 9. Select columns for output
  cpue_dt = cateffpen_dt[, .(fishery, year, month, obs, pred_cpue, se_log)]
  
  # 10. Add fleet names
  if(is.null(fishery_names)) {
    # Use default naming: S{fishery_number}_INDEX
    cpue_dt[, Fleet_name := sprintf("S%02d_INDEX", fishery)]
  } else {
    # Map fishery numbers to provided names
    fishery_map = data.table::data.table(
      fishery = seq_along(fishery_names),
      Fleet_name = fishery_names
    )
    cpue_dt = merge(cpue_dt, fishery_map, by = "fishery", all.x = TRUE)
    
    # For any unmapped fisheries, use default naming
    cpue_dt[is.na(Fleet_name), Fleet_name := sprintf("S%02d_INDEX", fishery)]
  }
  
  # 11. Add Use = 1 (all MFCL CPUE is used if penalty > 0)
  cpue_dt[, Use := 1L]
  
  # 12. Convert to standard output format
  cpue_dt[order(fishery,year,month)]
  cpue_dt[, Time := 1:.N, by = fishery]  # Create Time as sequential index within each fishery
  cpue_dt[, Dev := obs - pred_cpue]  # Deviation: observed - predicted
  cpue_dt[, Obs := obs]
  cpue_dt[, Exp := pred_cpue]
  cpue_dt[, SE := se_log]
  cpue_dt[, Fleet := fishery]
  
  # 13. Add id column
  cpue_dt[, id := model_id]
  
  # 14. Select and reorder columns to standard format
  cpue_dt = cpue_dt[, .(id, Fleet, Fleet_name, Time, Obs, Exp, SE, Dev, Use)]
  
  # Sort by Fleet and Time
  data.table::setorderv(cpue_dt, c("Fleet", "Time"))
  
  if(verbose) {
    message(sprintf("Extracted %d CPUE observations from %d fisheries", 
                    nrow(cpue_dt), 
                    length(unique(cpue_dt$Fleet))))
  }
  
  # 15. Write CSV (optional)
  if(save_csv) {
    output_file = file.path(output_dir, "cpue.csv")
    data.table::fwrite(cpue_dt, output_file)
    
    if(verbose) {
      message(sprintf("CPUE data written to %s", output_file))
    }
  }
  
  # 16. Return data.table
  return(cpue_dt)
}
