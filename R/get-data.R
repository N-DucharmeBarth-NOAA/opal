#' Get bundled model data
#'
#' Loads one of the packaged opal model data objects. This replaces the legacy
#' data-construction helper, which depended on historical raw inputs that are no
#' longer bundled with the package.
#'
#' @param model Character model identifier. Supported values are
#'   \code{"opal_baseline"}, \code{"opakapaka"}, and \code{"wcpo_bet"}.
#'   Aliases \code{"baseline"}, \code{"opaka"}, and \code{"bet"} are also
#'   accepted.
#' @param include_parameters Logical; if \code{TRUE}, return a list with both
#'   \code{data} and matching initial \code{parameters}.
#'
#' @return A data list ready for \code{\link{opal_model}}, or a list with
#'   elements \code{data} and \code{parameters} when
#'   \code{include_parameters = TRUE}.
#' @export
get_data <- function(model = c("opal_baseline", "opakapaka", "wcpo_bet"),
                     include_parameters = FALSE) {
  if (is.list(model)) {
    stop(
      "`get_data()` now loads bundled model data by name; pass ",
      "`model = \"opal_baseline\"`, `\"opakapaka\"`, or `\"wcpo_bet\"`.",
      call. = FALSE
    )
  }
  if (!is.logical(include_parameters) || length(include_parameters) != 1L) {
    stop("`include_parameters` must be a single logical value.", call. = FALSE)
  }

  model <- .resolve_bundled_model(model = model)
  data <- .load_bundled_data(model)

  if (isTRUE(include_parameters)) {
    return(list(data = data, parameters = .load_bundled_parameters(model)))
  }

  data
}

.resolve_bundled_model <- function(model = NULL, data = NULL) {
  aliases <- c(
    opal_baseline = "opal_baseline",
    baseline = "opal_baseline",
    opal = "opal_baseline",
    opakapaka = "opakapaka",
    opaka = "opakapaka",
    wcpo_bet = "wcpo_bet",
    bet = "wcpo_bet"
  )

  if (!is.null(model)) {
    if (length(model) > 1L) {
      model <- match.arg(model, choices = c("opal_baseline", "opakapaka", "wcpo_bet"))
    }
    model <- tolower(gsub("-", "_", as.character(model[[1L]]), fixed = TRUE))
    if (!model %in% names(aliases)) {
      stop(
        "`model` must be one of: opal_baseline, opakapaka, wcpo_bet.",
        call. = FALSE
      )
    }
    return(unname(aliases[[model]]))
  }

  if (is.null(data)) {
    return("opal_baseline")
  }
  if (!is.list(data)) {
    stop("`data` must be a model data list when `model` is not supplied.", call. = FALSE)
  }

  if (identical(data, .load_bundled_data("opal_baseline"))) {
    return("opal_baseline")
  }
  if (identical(data, .load_bundled_data("opakapaka"))) {
    return("opakapaka")
  }
  if (identical(data, .load_bundled_data("wcpo_bet"))) {
    return("wcpo_bet")
  }

  n_fishery <- as.integer(data$n_fishery %||% NA_integer_)
  n_year <- as.integer(data$n_year %||% NA_integer_)
  n_len <- as.integer(data$n_len %||% NA_integer_)
  n_age <- as.integer(data$n_age %||% NA_integer_)

  if (identical(c(n_fishery, n_year, n_len, n_age), c(3L, 75L, 17L, 44L))) {
    return("opakapaka")
  }
  if (identical(c(n_fishery, n_year, n_len, n_age), c(15L, 268L, 95L, 40L))) {
    return("wcpo_bet")
  }

  stop(
    "Could not infer the bundled parameter set from `data`; pass `model` ",
    "explicitly.",
    call. = FALSE
  )
}

.load_bundled_data <- function(model) {
  object <- switch(model,
    opal_baseline = "opal_baseline_data",
    opakapaka = "opaka_data",
    wcpo_bet = "wcpo_bet_data"
  )
  env <- new.env(parent = emptyenv())
  utils::data(list = object, package = "opal", envir = env)
  env[[object]]
}

.load_bundled_parameters <- function(model) {
  object <- switch(model,
    opal_baseline = "opal_baseline_parameters",
    opakapaka = "opaka_parameters",
    wcpo_bet = "wcpo_bet_parameters"
  )
  env <- new.env(parent = emptyenv())
  utils::data(list = object, package = "opal", envir = env)
  env[[object]]
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
