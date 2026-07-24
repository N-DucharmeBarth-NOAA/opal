#' Plot catch
#'
#' Plot catch by year and fishery.
#'
#' @param data A model data list passed to \code{MakeADFun}.
#' @param obj The AD object created using \code{MakeADFun}.
#' @param plot_resid Logical; plot catch residuals instead of observed and
#'   predicted catch.
#' @return A \code{ggplot2} object.
#' @import ggplot2
#' @import dplyr
#' @importFrom reshape2 melt
#' @importFrom rlang .data
#' @export
#'
plot_catch <- function(data, obj, plot_resid = FALSE) {
  yrs <- data$years
  fisheries <- paste0("Fishery: ", seq_len(data$n_fishery))

  df_obs <- melt(data$catch_obs_ysf, value.name = "obs") %>%
    filter(.data$obs > 0) %>%
    mutate(
      season = paste("Season:", 1),
      fishery = fisheries[.data$fishery],
      Type = "Observed",
      fishery = factor(.data$fishery, levels = fisheries)
    )

  df_pred <- obj$report()$catch_pred_ysf %>%
    melt(value.name = "pred") %>%
    mutate(
      year = yrs[.data$Var1],
      season = paste("Season:", .data$Var2),
      fishery = fisheries[.data$Var3]
    ) %>%
    right_join(df_obs, by = join_by("year", "season", "fishery")) %>%
    mutate(
      resid = .data$obs - .data$pred,
      fishery = factor(.data$fishery, levels = fisheries)
    )

  message("The maximum catch difference was: ", max(df_pred$resid))

  if (plot_resid) {
    p <- ggplot(df_pred, aes(x = .data$year, y = .data$resid)) +
      geom_point(color = "red") +
      labs(x = "Year", y = "Catch residual (tonnes)")
  } else {
    p <- ggplot(df_pred, aes(x = .data$year, y = .data$obs)) +
      geom_point(aes(color = .data$Type)) +
      geom_line(aes(y = .data$pred), group = 1) +
      labs(x = "Year", y = "Catch (tonnes)", color = NULL) +
      scale_y_continuous(
        limits = c(0, NA),
        expand = expansion(mult = c(0, 0.05))
      )
  }

  p + facet_wrap(fishery ~ .)
}

#' Plot CPUE
#'
#' Plot observed and predicted CPUE by season and fishery.
#'
#' @param data A model data list passed to \code{MakeADFun}.
#' @param object The AD object created using \code{MakeADFun}.
#' @return A \code{ggplot2} object.
#' @import ggplot2
#' @import dplyr
#' @importFrom scales pretty_breaks
#' @export
#'
plot_cpue <- function(data, object) {
  report <- object$report(object$env$last.par.best)

  df <- data$cpue_data %>%
    mutate(pred = report$cpue_pred, sigma = report$cpue_sigma)

  ggplot(df, aes(x = .data$year, y = .data$value)) +
    geom_point(aes(color = "Observed")) +
    geom_linerange(aes(
      ymin = exp(log(.data$value) - .data$sigma),
      ymax = exp(log(.data$value) + .data$sigma),
      color = "Observed"
    )) +
    geom_line(aes(y = .data$pred, color = "Predicted")) +
    labs(x = "Year", y = "CPUE", color = NULL) +
    scale_x_continuous(breaks = pretty_breaks()) +
    scale_y_continuous(
      limits = c(0, NA),
      expand = expansion(mult = c(0, 0.05))
    ) +
    facet_wrap(season ~ fishery)
}

#' Plot spawning biomass
#'
#' Plot spawning biomass or relative spawning biomass by year for one or more
#' model runs.
#'
#' @param data_list A list of model data lists passed to \code{MakeADFun}.
#' @param object_list A list of AD objects created using \code{MakeADFun}.
#' @param relative Logical; plot spawning biomass relative to unfished biomass.
#' @param labels Optional labels for the model runs.
#' @return A \code{ggplot2} object.
#' @import ggplot2
#' @import dplyr
#' @importFrom scales pretty_breaks
#' @export
#'
plot_biomass_spawning <- function(data_list, object_list, relative = TRUE,
                                  labels = NULL) {
  n_model <- length(data_list)
  if (is.null(labels)) labels <- seq_len(n_model)

  model_data <- vector("list", n_model)
  for (j in seq_len(n_model)) {
    years <- data_list[[j]]$first_yr:(data_list[[j]]$last_yr + 1)
    report <- object_list[[j]]$report()
    model_data[[j]] <- data.frame(
      Model = labels[j],
      year = years,
      B0 = report$B0,
      value = report$spawning_biomass_y
    )
  }
  df <- bind_rows(model_data) %>%
    mutate(Model = factor(.data$Model, levels = labels))

  if (relative) {
    df <- df %>% mutate(value = .data$value / .data$B0)
    ylab <- "Relative spawning biomass"
  } else {
    df <- df %>% mutate(value = .data$value / 1e6)
    ylab <- "Spawning biomass (millions of tonnes)"
  }

  p <- ggplot(df, aes(
    x = .data$year,
    y = .data$value,
    color = .data$Model
  )) +
    geom_line() +
    scale_x_continuous(breaks = pretty_breaks()) +
    scale_y_continuous(
      limits = c(0, NA),
      expand = expansion(mult = c(0, 0.05))
    ) +
    labs(x = "Year", y = ylab)

  if (n_model == 1) p <- p + theme(legend.position = "none")
  p
}

#' Plot initial numbers
#'
#' Plot initial population numbers by age.
#'
#' @param data A model data list passed to \code{MakeADFun}.
#' @param object The AD object created using \code{MakeADFun}.
#' @return A \code{ggplot2} object.
#' @import ggplot2
#' @export
#'
plot_initial_numbers <- function(data, object) {
  ages <- data$min_age:data$max_age
  number_ysa <- object$report()$number_ysa
  initial_numbers <- data.frame(age = ages, value = number_ysa[1, 1, ])

  ggplot(initial_numbers, aes(x = .data$age, y = .data$value / 1e6)) +
    geom_line(linetype = "dashed") +
    labs(x = "Age", y = "Initial numbers (millions)") +
    scale_x_continuous(
      limits = c(0, NA),
      expand = expansion(mult = c(0, 0.05))
    ) +
    scale_y_continuous(
      limits = c(0, NA),
      expand = expansion(mult = c(0, 0.05))
    )
}

#' Plot harvest rate
#'
#' Plot harvest rate by year, season, and age.
#'
#' @param data A model data list passed to \code{MakeADFun}.
#' @param object The AD object created using \code{MakeADFun}.
#' @param years Optional years to show. The default plots every model year from
#'   the first catch year onwards.
#' @param ... Options passed to \code{geom_density_ridges}.
#' @return A \code{ggplot2} object.
#' @import ggplot2
#' @import dplyr
#' @importFrom reshape2 melt
#' @importFrom ggridges geom_density_ridges
#' @importFrom scales pretty_breaks
#' @export
#'
plot_hrate <- function(data, object, years = NULL, ...) {
  model_years <- data$first_yr:data$last_yr
  ages <- data$min_age:data$max_age
  if (is.null(years)) years <- model_years

  df <- object$report()$hrate_ysa %>%
    melt() %>%
    mutate(
      year = model_years[.data$Var1],
      season = .data$Var2,
      age = ages[.data$Var3]
    ) %>%
    filter(.data$year %in% years, .data$year >= data$first_yr_catch)

  ggplot(df, aes(
    x = .data$age,
    y = .data$year,
    height = .data$value,
    group = .data$year
  )) +
    geom_density_ridges(
      stat = "identity",
      alpha = 0.75,
      rel_min_height = 0,
      color = NA,
      ...
    ) +
    facet_wrap(season ~ .) +
    labs(x = "Age", y = "Year") +
    scale_x_continuous(
      limits = c(0, NA),
      expand = expansion(mult = c(0, 0.05))
    ) +
    scale_y_reverse(breaks = pretty_breaks())
}

#' Plot recruitment
#'
#' Plot recruitment by year.
#'
#' @param data A model data list passed to \code{MakeADFun}.
#' @param object The AD object created using \code{MakeADFun}.
#' @return A \code{ggplot2} object.
#' @import ggplot2
#' @export
#'
plot_recruitment <- function(data, object) {
  report <- object$report(object$env$last.par.best)
  years <- data$first_yr:(data$last_yr + 1)
  recruitment <- data.frame(year = years, value = report$number_ysa[, 1, 1])

  ggplot(recruitment, aes(x = .data$year, y = .data$value / 1e6)) +
    geom_hline(yintercept = report$R0 / 1e6, linetype = "dashed") +
    geom_line(linetype = "dashed") +
    labs(x = "Year", y = "Recruitment (millions)") +
    scale_y_continuous(
      limits = c(0, NA),
      expand = expansion(mult = c(0, 0.05))
    )
}
