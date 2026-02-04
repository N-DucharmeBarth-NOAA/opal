#' Extract Selectivity-at-Length from SS3 Model
#'
#' Extracts final year selectivity curves for all fleets from Stock Synthesis
#' output and writes to standardized CSV format.
#'
#' @param model_dir Character. Path to SS3 model directory containing Report.sso
#' @param model_id Character. Model identifier for output
#' @param verbose Logical. Print progress messages? Default TRUE
#' 
#' @return data.table with columns: id, Fleet, Fleet_name, Yr, Sex, variable, value
#' 
#' @details
#' Reads SS3 output using r4ss::SS_output() and extracts length-based selectivity
#' from the sizeselex component for the final model year. Converts wide format
#' (one column per length bin) to long format for plotting compatibility.
#'
#' @examples
#' \dontrun{
#'   selex_dt = extract_ss3_selectivity(
#'     model_dir = "model-files/ss3/01-bet-base",
#'     model_id = "01-bet-base"
#'   )
#'   
#'   # Check output
#'   head(selex_dt)
#'   
#'   # Verify CSV written
#'   file.exists("model-files/ss3/01-bet-base/selex_l.csv")
#' }
#'
#' @export
extract_ss3_selectivity = function(model_dir, model_id, verbose = TRUE) {
  # Validate inputs
  if(!file.exists(file.path(model_dir, "Report.sso"))) {
    stop(sprintf("Report.sso not found in %s", model_dir))
  }
  
  if(verbose) {
    message(sprintf("Reading SS3 model from %s", model_dir))
  }
  
  # Read SS3 output
  tmp_report = r4ss::SS_output(model_dir, verbose = FALSE, printstats = FALSE)
  
  # Extract fleet definitions
  tmp_flt = tmp_report$definitions
  
  if(verbose) {
    message(sprintf("Extracting selectivity for final year %d", tmp_report$endyr))
  }
  
  # Extract and transform selectivity data
  # Filter for length selectivity (Factor=="Lsel") at final year
  tmp_len_selex = as.data.table(tmp_report$sizeselex) %>%
    .[Factor == "Lsel"] %>%
    .[Yr == tmp_report$endyr] %>%
    .[, Factor := NULL] %>%
    .[, Label := NULL] %>%
    # Convert from wide to long format
    melt(., id.vars = c("Fleet", "Yr", "Sex")) %>%
    # Add model identifier
    .[, id := model_id] %>%
    # Add fleet names from definitions
    .[, Fleet_name := tmp_flt$Fleet_name[Fleet]] %>%
    # Select and reorder columns to match standard format
    .[, .(id, Fleet, Fleet_name, Yr, Sex, variable, value)]
  
  # Convert variable (length bins) to numeric
  tmp_len_selex[, variable := as.numeric(as.character(variable))]
  
  # Write CSV output
  output_file = file.path(model_dir, "selex_l.csv")
  data.table::fwrite(tmp_len_selex, output_file)
  
  if(verbose) {
    message(sprintf("Written %d rows to %s", nrow(tmp_len_selex), output_file))
    message(sprintf("  Fleets: %s", paste(unique(tmp_len_selex$Fleet), collapse = ", ")))
    message(sprintf("  Length range: %.1f - %.1f cm", 
                    min(tmp_len_selex$variable), 
                    max(tmp_len_selex$variable)))
  }
  
  return(tmp_len_selex)
}
