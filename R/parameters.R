#' Get bundled initial parameter values
#'
#' Loads the packaged initial parameter list that matches one of the bundled
#' opal model data sets.
#'
#' @param data Optional model data list used to infer the bundled parameter set.
#' @param model Optional character model identifier. Supported values are
#'   \code{"opal_baseline"}, \code{"opakapaka"}, and \code{"wcpo_bet"}.
#'   Aliases \code{"baseline"}, \code{"opaka"}, and \code{"bet"} are also
#'   accepted. Supplying \code{model} is preferred when \code{data} has been
#'   modified after loading.
#'
#' @return A \code{list} of initial parameter values.
#' @export
#' 
get_parameters <- function(data = NULL, model = NULL) {
  if (is.null(model) && is.character(data) && length(data) == 1L) {
    model <- data
    data <- NULL
  }

  model <- .resolve_bundled_model(model = model, data = data)
  .load_bundled_parameters(model)
}

.match_selectivity_age_indices <- function(target_ages, source_ages) {
  if (length(source_ages) == 0L) {
    stop("`source_ages` must contain at least one age.", call. = FALSE)
  }
  source_ages <- sort(unique(as.integer(source_ages)))
  target_ages <- as.integer(target_ages)
  pmin(pmax(target_ages, min(source_ages)), max(source_ages))
}

.match_selectivity_year_indices <- function(target_years, first_yr, n_year) {
  if (length(target_years) == 0L) {
    return(integer())
  }
  if (n_year < 1L) {
    stop("`n_year` must be positive.", call. = FALSE)
  }
  last_yr <- first_yr + n_year - 1L
  matched_years <- pmin(pmax(as.integer(target_years), first_yr), last_yr)
  matched_years - first_yr + 1L
}

#' Get default parameter mapping
#' 
#' Get a default parameter mapping. Parameter mapping is used by \code{MakeADFun} 
#' to turn parameters on/off.
#' 
#' @param parameters a \code{list} containing the initial parameter values to be passed to \code{MakeADFun}.
#' @return a named \code{list} of parameter mapping.
#' @export
#' 
get_map <- function(parameters) {
  fixed <- intersect(
    c(
      "log_h", "log_sigma_r", "cpue_creep", "log_cpue_tau",
      "log_cpue_omega", "log_lf_tau", "log_wf_tau", "log_L1", "log_L2",
      "log_k", "log_CV1", "log_CV2", "par_sel", "log_init_F_f",
      "init_rdev_a"
    ),
    names(parameters)
  )

  lapply(parameters[fixed], function(x) factor(rep(NA, length(x))))
}

#' Get default parameter bounds
#' 
#' Get \code{data.frame} of default parameter bounds.
#' 
#' @param obj a \code{list} specifying the AD object created using the \code{MakeADFun} function.
#' @param parameters The parameter list used to construct \code{obj}. Retained
#'   for API compatibility.
#' @return a \code{data.frame} of parameter bounds.
#' @importFrom RTMB qlogis
#' @export
#' 
get_bounds <- function(obj, parameters) {
  lower <- rep(-Inf, length(obj$par))
  upper <- rep(Inf, length(obj$par))
  parameter_names <- names(obj$par)

  set_bounds <- function(pattern, lower_value, upper_value) {
    index <- grep(pattern, parameter_names, fixed = TRUE)
    lower[index] <<- lower_value
    upper[index] <<- upper_value
  }

  set_bounds("log_B0", log(1), 22)
  set_bounds("log_h", log(0.21), log(1))
  set_bounds("log_sigma_r", log(0.1), log(2))
  set_bounds("log_cpue_q", log(0.001), log(10))
  set_bounds("log_lf_tau", -9, 9)
  set_bounds("log_wf_tau", -9, 9)
  set_bounds("log_init_F_f", log(1e-12), log(3))
  set_bounds("rdev_y", -5, 5)
  set_bounds("init_rdev_a", -5, 5)
  set_bounds("par_sel", -7, 7)

  check_bounds(opt = obj, lower = lower, upper = upper)

  data.frame(
    parameter = parameter_names,
    init = obj$par,
    lower = lower,
    upper = upper
  )
}

#' Check if parameters are up against the bounds
#' 
#' A \code{data.frame} containing the parameters that are against their lower or 
#' upper bound.
#' 
#' @param opt an optimized TMB object.
#' @param lower a vector of lower bounds.
#' @param upper a vector of upper bounds.
#' @return a \code{data.frame}.
#' @export
#' 
check_bounds <- function(opt, lower, upper) {
  df <- data.frame(par = names(opt$par), lb = lower, value = opt$par, ub = upper) %>%
    mutate(index = 1:n())
  rownames(df) <- NULL
  ilb <- which(df$value <= df$lower)
  iub <- which(df$value >= df$upper)
  return(df[c(ilb, iub), ])
}
