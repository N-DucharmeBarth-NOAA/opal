#' Extract CPUE Index Fit from MFCL Model
#'
#' Extracts observed and predicted CPUE values from MFCL report file,
#' calculates standard errors from penalty values in par file, and writes
#' to standardized CSV format.
#'
#' @param rep_file Character. Path to MFCL plot-*.rep file
#' @param par_file Character. Path to MFCL *.par file
#' @param model_id Character. Model identifier for output
#' @param first_year Numeric. First year of model. Default 1952
#' @param fishery_names Character vector. Optional fleet names. If NULL, uses default naming
#' @param output_dir Character. Directory for output CSV
#' @param verbose Logical. Print progress messages? Default TRUE
#' 
#' @return data.table with columns: id, Fleet, Fleet_name, Time, Obs, Exp, SE, Dev, Use
#' 
#' @details
#' MFCL stores CPUE in log space and penalties separately:
#' 1. Extract observed CPUE: FLR4MFCL::cpue_obs() - in log space, need exp()
#' 2. Extract predicted CPUE: FLR4MFCL::cpue_pred() - in log space, need exp()
#' 3. Extract penalties from par file: flag==92, flagtype<0
#' 4. Calculate CV from penalty: CV = 1/sqrt(2*penalty/10)
#' 5. Calculate SE: SE = sqrt(log(1 + CV^2))
#' 6. Merge and calculate deviations
#' 
#' Time is calculated as: year + (season-1)/4 for quarterly models
#'
#' @examples
#' \dontrun{
#'   cpue_dt = extract_mfcl_cpue(
#'     rep_file = "model-files/mfcl/v11/plot-10.par.rep",
#'     par_file = "model-files/mfcl/v11/10.par",
#'     model_id = "mfcl-v11",
#'     output_dir = "model-files/mfcl/v11"
#'   )
#' }
#'
#' @export
extract_mfcl_cpue = function(rep_file, par_file, model_id,
                             first_year = 1952,
                             fishery_names = NULL,
                             output_dir,
                             verbose = TRUE) {
  # Validate inputs
  if(!file.exists(rep_file)) {
    stop(sprintf("Report file not found: %s", rep_file))
  }
  
  if(!file.exists(par_file)) {
    stop(sprintf("Par file not found: %s", par_file))
  }
  
  if(verbose) {
    message(sprintf("Reading MFCL files: %s and %s...", 
                    basename(rep_file), basename(par_file)))
  }
  
  # 1. Read MFCL files
  tmp.rep = FLR4MFCL::read.MFCLRep(rep_file)
  tmp.par = FLR4MFCL::read.MFCLPar(par_file, first.yr = first_year)
  
  if(verbose) {
    message("Extracting observed CPUE...")
  }
  
  # 2. Extract observed CPUE (in log space)
  tmp_cpue_obs = data.table::as.data.table(FLR4MFCL::cpue_obs(tmp.rep))
  tmp_cpue_obs = tmp_cpue_obs[, ts := as.numeric(year) + (as.numeric(season)-1)/4]
  tmp_cpue_obs = tmp_cpue_obs[, .(unit, ts, value = exp(value))]  # Convert from log space
  data.table::setnames(tmp_cpue_obs, c("unit", "ts", "value"), c("fishery", "ts", "cpue_obs"))
  
  if(verbose) {
    message("Extracting predicted CPUE...")
  }
  
  # 3. Extract predicted CPUE (in log space)
  tmp_cpue_pred = data.table::as.data.table(FLR4MFCL::cpue_pred(tmp.rep))
  tmp_cpue_pred = tmp_cpue_pred[, ts := as.numeric(year) + (as.numeric(season)-1)/4]
  tmp_cpue_pred = tmp_cpue_pred[, .(unit, ts, value = exp(value))]  # Convert from log space
  data.table::setnames(tmp_cpue_pred, c("unit", "ts", "value"), c("fishery", "ts", "cpue_pred"))
  
  if(verbose) {
    message("Extracting penalties...")
  }
  
  # 4. Extract penalty values from par file
  tmp_penalty = data.table::as.data.table(FLR4MFCL::flags(tmp.par))
  tmp_penalty = tmp_penalty[flag == 92 & flagtype < 0]
  tmp_penalty = tmp_penalty[, flagtype := -flagtype]
  data.table::setnames(tmp_penalty, c("flagtype", "value"), c("fishery", "penalty"))
  tmp_penalty = tmp_penalty[, .(fishery, penalty)]
  tmp_penalty = tmp_penalty[, penalty := penalty/10]
  
  if(verbose) {
    message("Merging and calculating standard errors...")
  }
  
  # 5. Merge observed and predicted
  tmp_cpue = merge(tmp_cpue_obs, tmp_cpue_pred, by = c("fishery", "ts"))
  tmp_cpue = tmp_cpue[, fishery := as.numeric(fishery)]
  
  # 6. Merge with penalties and calculate SE
  cpue_dt = merge(tmp_cpue, tmp_penalty, by = "fishery")
  cpue_dt = cpue_dt[, input_cv := 1/(sqrt(2*penalty))]
  cpue_dt = cpue_dt[, se_log := sqrt(log(1 + input_cv^2))]
  cpue_dt = cpue_dt[, Dev := cpue_obs - cpue_pred]
  
  # 7. Filter to fisheries with penalties > 0 (i.e., CPUE data exists)
  cpue_dt = cpue_dt[penalty > 0]
  
  if(nrow(cpue_dt) == 0) {
    stop("No CPUE data found (no fisheries with penalty > 0)")
  }
  
  # 8. Add fleet names
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
  
  # 9. Add Use = 1 (all MFCL CPUE is used if penalty > 0)
  cpue_dt[, Use := 1L]
  
  # 10. Rename columns to standard format
  data.table::setnames(cpue_dt, 
                       c("fishery", "ts", "cpue_obs", "cpue_pred", "se_log"),
                       c("Fleet", "Time", "Obs", "Exp", "SE"))
  
  # 11. Add id column
  cpue_dt[, id := model_id]
  
  # 12. Select and reorder columns to standard format
  cpue_dt = cpue_dt[, .(id, Fleet, Fleet_name, Time, Obs, Exp, SE, Dev, Use)]
  
  # Sort by Fleet and Time
  data.table::setorderv(cpue_dt, c("Fleet", "Time"))
  
  if(verbose) {
    message(sprintf("Extracted %d CPUE observations from %d fleets", 
                    nrow(cpue_dt), 
                    length(unique(cpue_dt$Fleet))))
  }
  
  # 13. Write CSV
  output_file = file.path(output_dir, "cpue.csv")
  data.table::fwrite(cpue_dt, output_file)
  
  if(verbose) {
    message(sprintf("CPUE data written to %s", output_file))
  }
  
  # 14. Return data.table
  return(cpue_dt)
}
