#' Modular CPUE comparison plotter (project-style)
#'
#' Create a CPUE comparison plot from pre-processed CPUE data without requiring
#' model directory scanning or variance adjustments. This follows the style of
#' `plot_selectivity_modular()` and extracts core plotting logic from
#' `plot_index_comparison()`.
#'
#' @param cpue_input Either a `data.table` already read from `cpue.csv`, or
#'   a character path to a `cpue.csv` file.
#' @param fleet_names Optional character vector of fleet names to assign to
#'   unique fleet IDs (length must match number of unique Fleet ids when used).
#' @param model_id Optional model identifier used to label `model_label` column
#'   if not present in the input data.
#' @param show_se Logical; show error bars for observed CPUE values.
#' @param use_log_scale Logical; apply log transformation to y-axis.
#' @param n_col Number of columns for `facet_wrap`.
#' @param verbose Logical; print progress messages.
#'
#' @return A `ggplot` object (CPUE comparison plot facetted by `Fleet_name`).
#' @export
#' @examples
#' \dontrun{
#' # Read CPUE data from two models
#' mfcl_cpue = data.table::fread("model-files/mfcl/v11/cpue.csv")
#' mfcl_cpue[, model_label := "MFCL"]
#' 
#' ss3_cpue = data.table::fread("model-files/ss3/01-bet-base/cpue.csv")
#' ss3_cpue[, model_label := "SS3"]
#' 
#' # Combine and plot
#' cpue_combined = data.table::rbindlist(list(mfcl_cpue, ss3_cpue), 
#'                                        use.names = TRUE, fill = TRUE)
#' p = plot_cpue_modular(cpue_combined, show_se = TRUE, use_log_scale = FALSE)
#' print(p)
#' 
#' # Or read from file path
#' p2 = plot_cpue_modular("model-files/ss3/01-bet-base/cpue.csv")
#' print(p2)
#' }
plot_cpue_modular = function(cpue_input,
                             fleet_names = NULL,
                             model_id = NULL,
                             show_se = TRUE,
                             use_log_scale = FALSE,
                             n_col = 2,
                             verbose = TRUE) {
  # Read or coerce input
  if (is.character(cpue_input) && length(cpue_input) == 1) {
    if (!file.exists(cpue_input)) stop(sprintf("CPUE file not found: %s", cpue_input))
    cpue_dt = data.table::fread(cpue_input)
  } else if (data.table::is.data.table(cpue_input) || is.data.frame(cpue_input)) {
    cpue_dt = data.table::as.data.table(cpue_input)
  } else {
    stop("cpue_input must be a file path or a data.table/data.frame")
  }
  
  # Normalize column names (handle both upper and lower case)
  col_map = list(
    fleet = c("fleet", "Fleet"),
    time = c("time", "Time"),
    obs = c("obs", "Obs"),
    exp = c("exp", "Exp"),
    se = c("se", "SE"),
    Fleet_name = c("Fleet_name", "fleetname")
  )
  
  for (target_col in names(col_map)) {
    possible_cols = col_map[[target_col]]
    found_col = intersect(possible_cols, colnames(cpue_dt))
    if (length(found_col) > 0) {
      # Use the first matching column and rename it
      if (found_col[1] != target_col) {
        data.table::setnames(cpue_dt, found_col[1], target_col, skip_absent = TRUE)
      }
    }
  }
  
  # Validate required columns
  required_cols = c("fleet", "time", "obs", "exp", "se")
  missing_cols = setdiff(required_cols, colnames(cpue_dt))
  if (length(missing_cols) > 0) {
    stop(sprintf("cpue_input missing required columns: %s", paste(missing_cols, collapse = ", ")))
  }
  
  # Ensure fleet is integer-like
  cpue_dt[, fleet := as.integer(fleet)]
  
  # Handle Fleet_name column
  unique_fleets = sort(unique(cpue_dt$fleet))
  if (!is.null(fleet_names)) {
    if (length(fleet_names) != length(unique_fleets)) {
      stop("Length of fleet_names must match number of unique fleet IDs in cpue input")
    }
    fleet_map = data.table::data.table(fleet = unique_fleets, Fleet_name = fleet_names)
  } else {
    if (!"Fleet_name" %in% colnames(cpue_dt)) {
      fleet_map = data.table::data.table(fleet = unique_fleets, 
                                        Fleet_name = paste0("Fleet_", unique_fleets))
    } else {
      fleet_map = unique(cpue_dt[, .(fleet, Fleet_name)])
    }
  }
  
  # Add model_label if missing
  if (!"model_label" %in% colnames(cpue_dt)) {
    cpue_dt[, model_label := ifelse(is.null(model_id), "model", as.character(model_id))]
  }
  
  # Merge fleet names (select only needed columns first to avoid conflicts)
  cpue_dt = merge(cpue_dt[, .SD, .SDcols = c("model_label", "fleet", "time", "obs", "exp", "se")], 
                  fleet_map, by = "fleet", all.x = TRUE, sort = FALSE)
  
  if (verbose) {
    message(sprintf("Preparing CPUE plot for %d fleets", data.table::uniqueN(cpue_dt$fleet)))
  }
  
  # Aggregate over any additional grouping variables not specified
  # For observed data, take first occurrence to avoid duplication
  obs_dt = cpue_dt[, .SD[1], by = .(Fleet_name, time)]
  
  # For fitted values, aggregate by model_label
  fit_dt = cpue_dt[, .(exp = mean(exp, na.rm = TRUE)), 
                   by = .(model_label, Fleet_name, fleet, time)]
  
  # Apply transformations and calculate error bars
  if (use_log_scale) {
    # Apply log transformation first
    obs_dt[, obs := log(obs)]
    fit_dt[, exp := log(exp)]
    # For log scale, use additive SE
    obs_dt[, lse := obs - se]
    obs_dt[, use := obs + se]
    ylab_txt = "Index (log-scale)"
    yint = 0
  } else {
    # For normal scale, use multiplicative SE
    obs_dt[, lse := exp(log(obs) - se)]
    obs_dt[, use := exp(log(obs) + se)]
    ylab_txt = "Index"
    yint = 1
  }
  
  # Create plot
  p = ggplot2::ggplot()
  
  # Add error bars if requested
  if (show_se) {
    p = p + ggplot2::geom_errorbar(data = obs_dt, 
                                   ggplot2::aes(x = time, ymin = lse, ymax = use),
                                   width = 0.2, linewidth = 0.7, color = "gray60")
  }
  
  # Add observed points
  p = p + ggplot2::geom_point(data = obs_dt, 
                              ggplot2::aes(x = time, y = obs),
                              shape = 21, size = 3, stroke = 0.5, fill = "gray90")
  
  # Add fitted lines
  p = p + ggplot2::geom_line(data = fit_dt, 
                             ggplot2::aes(x = time, y = exp, color = model_label, group = model_label),
                             linewidth = 1)
  
  # Add reference line
  p = p + ggplot2::geom_hline(yintercept = yint, linetype = "dashed", linewidth = 1)
  
  # Add faceting and labels
  p = p + ggplot2::facet_wrap(~ Fleet_name, ncol = n_col)
  p = p + ggplot2::xlab("Time")
  p = p + ggplot2::ylab(ylab_txt)
  
  # Set y-axis limits
  if (!use_log_scale) {
    p = p + ggplot2::coord_cartesian(ylim = c(0, NA))
  }
  
  # Apply color scale
  p = p + viridis::scale_color_viridis(name = "Model", begin = 0.1, end = 0.9,
                                       direction = -1, option = "H", discrete = TRUE)
  
  # Apply theme
  p = p + ggplot2::theme(
    text = ggplot2::element_text(size = 14),
    panel.background = ggplot2::element_rect(fill = "white", color = "black", linewidth = 0.8),
    panel.grid.major = ggplot2::element_line(color = "gray70", linetype = "dotted", linewidth = 0.3),
    panel.grid.minor = ggplot2::element_line(color = "gray85", linetype = "dotted", linewidth = 0.2),
    strip.background = ggplot2::element_rect(fill = "white"),
    legend.key = ggplot2::element_rect(fill = "white"),
    legend.position = "right"
  )
  
  return(p)
}
