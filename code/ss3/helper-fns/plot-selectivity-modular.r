library(data.table)
library(ggplot2)

#' Modular selectivity plotter (project-style)
#'
#' Create a selectivity plot from a `selex_l`-style input without requiring
#' `fleet_summary.csv`. This follows the style of the `extract_*_selectivity` helpers
#' used elsewhere in the project: validate inputs, support either a path or a
#' pre-read `data.table`, allow optional fleet name mapping, aggregate over
#' `Yr` and `Sex`, and return a `ggplot` object.
#'
#' @param selex_input Either a `data.table` already read from `selex_l.csv`, or
#'   a character path to a `selex_l.csv` file.
#' @param fleet_names Optional character vector of fleet names to assign to
#'   unique `Fleet` IDs (length must match number of unique Fleet ids when used).
#' @param model_id Optional model identifier used to label `model_label` column.
#' @param n_col Number of columns for `facet_wrap`.
#' @param verbose Logical; print progress messages.
#'
#' @return A `ggplot` object (selectivity curves facetted by `Fleet_name`).
#' @export
## Example (don't run)
#' \dontrun{
#' # Read selectivity CSVs from two model folders (MFCL and SS3), combine, and plot
#' mfcl_selex = file.path(this.path::this.proj(), "model-files", "mfcl", "v11", "selex_l.csv")
#' ss3_selex = file.path(this.path::this.proj(), "model-files", "ss3", "01-bet-base", "selex_l.csv")
#'
#' # Read both files (if they exist), add model labels, and combine
#' selex_list = list()
#' if (file.exists(mfcl_selex)) selex_list[["mfcl"]] = data.table::fread(mfcl_selex)[, model_label := "mfcl_v11"]
#' if (file.exists(ss3_selex)) selex_list[["ss3"]] = data.table::fread(ss3_selex)[, model_label := "ss3_01-bet-base"]
#'
#' if (length(selex_list) > 0) {
#'   selex_dt = data.table::rbindlist(selex_list, use.names = TRUE, fill = TRUE)
#'   # Plot (ignores Yr and Sex by design)
#'   p = plot_selectivity_modular(selex_dt)
#'   print(p)
#' }
#' }
plot_selectivity_modular = function(selex_input,
                                     fleet_names = NULL,
                                     model_id = NULL,
                                     n_col = 4,
                                     verbose = TRUE) {
  # Read or coerce input
  if (is.character(selex_input) && length(selex_input) == 1) {
    if (!file.exists(selex_input)) stop(sprintf("selex file not found: %s", selex_input))
    selex_dt = data.table::fread(selex_input)
  } else if (is.data.table(selex_input) || is.data.frame(selex_input)) {
    selex_dt = as.data.table(selex_input)
  } else {
    stop("selex_input must be a file path or a data.table/data.frame")
  }

  # Basic validation
  required_cols = c("Fleet", "variable", "value")
  missing_cols = setdiff(required_cols, colnames(selex_dt))
  if (length(missing_cols) > 0) {
    stop(sprintf("selex input missing required columns: %s", paste(missing_cols, collapse = ", ")))
  }

  # Ensure Fleet is integer-like
  selex_dt[, Fleet := as.integer(Fleet)]

  # If Fleet_name absent, create default mapping from Fleet (or use provided fleet_names)
  unique_fleets = sort(unique(selex_dt$Fleet))
  if (!is.null(fleet_names)) {
    if (length(fleet_names) != length(unique_fleets)) {
      stop("Length of fleet_names must match number of unique Fleet IDs in selex input")
    }
    fleet_map = data.table(Fleet = unique_fleets, Fleet_name = fleet_names)
  } else {
    if (!"Fleet_name" %in% colnames(selex_dt)) {
      fleet_map = data.table(Fleet = unique_fleets, Fleet_name = paste0("Fleet_", unique_fleets))
    } else {
      fleet_map = unique(selex_dt[, .(Fleet, Fleet_name)])
    }
  }

  # Add model_label if missing
  if (!"model_label" %in% colnames(selex_dt)) {
    selex_dt[, model_label := ifelse(is.null(model_id), "model", as.character(model_id))]
  }

  # Merge fleet names
  selex_dt = merge(selex_dt[,.(model_label, Fleet, variable, value)], fleet_map, by = "Fleet", all.x = TRUE, sort = FALSE)

  if (verbose) {
    message(sprintf("Preparing selectivity plot for %d fleets", uniqueN(selex_dt$Fleet)))
  }

  # Aggregate over Yr and Sex (ignore them) by averaging values
  agg_dt = selex_dt[, .(value = mean(value, na.rm = TRUE)), by = .(model_label, Fleet, Fleet_name, variable)]

  # Create plot (styled to match project example)
  p = ggplot(agg_dt, aes(x = variable, y = value, color = model_label, group = model_label)) +
    geom_line(linewidth = 1.2) +
    facet_wrap(~ Fleet_name, scales = "free_x", ncol = n_col) +
    coord_cartesian(ylim = c(0, 1), expand = FALSE) +
    xlab("Length") + ylab("Selectivity") +
    viridis::scale_color_viridis(name = "Model", begin = 0.1, end = 0.9, direction = -1, option = "H", discrete = TRUE) +
    theme(
      text = element_text(size = 14),
      panel.background = element_rect(fill = "white", color = "black", linewidth = 0.8),
      panel.grid.major = element_line(color = "gray70", linetype = "dotted", linewidth = 0.3),
      panel.grid.minor = element_line(color = "gray85", linetype = "dotted", linewidth = 0.2),
      strip.background = element_rect(fill = "white"),
      legend.key = element_rect(fill = "white"),
      legend.position = "right"
    )

  return(p)
}



