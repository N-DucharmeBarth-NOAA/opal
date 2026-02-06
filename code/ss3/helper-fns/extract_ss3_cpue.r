#' Extract CPUE Index Fit from SS3 Model
#'
#' Extracts observed and predicted CPUE values with standard errors from
#' Stock Synthesis output and optionally writes to standardized CSV format 
#' compatible with plot_index_comparison().
#'
#' @param model_dir Character. Path to SS3 model directory containing Report.sso
#' @param model_id Character. Model identifier for output
#' @param save_csv Logical. Save output as CSV file? Default TRUE
#' @param verbose Logical. Print progress messages? Default TRUE
#' 
#' @return data.table with columns: id, Fleet, Fleet_name, Time, Obs, Exp, SE, Dev, Use
#' 
#' @details
#' Reads SS3 output using r4ss::SS_output() and extracts CPUE fit data.
#' The cpue component already contains observed, expected, SE, and deviation values.
#' Time is in decimal years (year + (season-1)/nseasons).
#'
#' @examples
#' \dontrun{
#'   # Extract and save CSV
#'   cpue_dt = extract_ss3_cpue(
#'     model_dir = "model-files/ss3/01-bet-base",
#'     model_id = "01-bet-base"
#'   )
#'   
#'   # Extract without saving CSV
#'   cpue_dt = extract_ss3_cpue(
#'     model_dir = "model-files/ss3/01-bet-base",
#'     model_id = "01-bet-base",
#'     save_csv = FALSE
#'   )
#'   
#'   # Check output
#'   head(cpue_dt)
#' }
#'
#' @export
extract_ss3_cpue = function(model_dir, model_id, save_csv = TRUE, verbose = TRUE) {
  # Validate inputs
  if(!file.exists(file.path(model_dir, "Report.sso"))) {
    stop(sprintf("Report.sso not found in %s", model_dir))
  }
  
  if(verbose) {
    message(sprintf("Reading SS3 model from %s...", model_dir))
  }
  
  # 1. Read SS3 output
  tmp_report = r4ss::SS_output(model_dir, verbose = FALSE, printstats = FALSE)
  
  # Check if cpue component exists
  if(is.null(tmp_report$cpue)) {
    stop("No CPUE data found in SS3 output")
  }
  
  if(verbose) {
    message("Extracting CPUE data...")
  }
  
  # 2. Extract CPUE data
  cpue_dt = data.table::as.data.table(tmp_report$cpue)
  
  # 3. Select and rename columns to standard format, convert Fleet and Time to numeric
  cpue_dt = cpue_dt[, .(Fleet = as.numeric(Fleet), Fleet_name, Time = as.numeric(Time), Obs, Exp, SE, Dev, Use)]
  
  # 4. Add model id
  cpue_dt[, id := model_id]
  
  # 5. Reorder columns to standard format
  cpue_dt = cpue_dt[, .(id, Fleet, Fleet_name, Time, Obs, Exp, SE, Dev, Use)]
  
  if(verbose) {
    message(sprintf("Extracted %d CPUE observations from %d fleets", 
                    nrow(cpue_dt), 
                    length(unique(cpue_dt$Fleet))))
  }
  
  # 6. Write CSV (optional)
  if(save_csv) {
    output_file = file.path(model_dir, "cpue.csv")
    data.table::fwrite(cpue_dt, output_file)
    
    if(verbose) {
      message(sprintf("CPUE data written to %s", output_file))
    }
  }
  
  # 7. Return data.table
  return(cpue_dt)
}
