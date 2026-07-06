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
  map <- list()
  map[["par_log_psi"]] <- factor(NA)
  map[["par_log_m0"]] <- factor(NA)
  map[["par_log_m10"]] <- factor(NA)
  map[["par_log_h"]] <- factor(NA)
  map[["par_log_sigma_r"]] <- factor(NA)
  map[["par_log_cpue_sigma"]] <- factor(NA)
  map[["par_log_cpue_omega"]] <- factor(NA)
  map[["par_cpue_creep"]] <- factor(NA)
  map[["par_log_af_alpha"]] <- factor(rep(NA, 2))
  map[["par_log_lf_alpha"]] <- factor(rep(NA, 5))
  map[["par_sel_rho_y"]] <- factor(rep(NA, length(parameters$par_sel_rho_y)))
  map[["par_sel_rho_a"]] <- factor(rep(NA, length(parameters$par_sel_rho_a)))
  map[["par_log_sel_sigma"]] <- factor(rep(NA, length(parameters$par_log_sel_sigma)))
  map[["par_log_sel_4"]] <- factor(matrix(NA, nrow = nrow(parameters$par_log_sel_4), ncol = ncol(parameters$par_log_sel_4)))
  map[["log_init_F_f"]] <- factor(rep(NA, length(parameters$log_init_F_f)))
  map[["init_rdev_a"]] <- factor(rep(NA, length(parameters$init_rdev_a)))
  # map[["par_rec_dev_y"]] <- factor(rep(NA, length(parameters$par_rdev_y)))
  return(map)
}

#' Get default parameter bounds
#' 
#' Get \code{data.frame} of default parameter bounds.
#' 
#' @param obj a \code{list} specifying the AD object created using the \code{MakeADFun} function.
#' @param parameters a \code{list} specifying the AD object created using the \code{MakeADFun} function.
#' @return a \code{data.frame} of parameter bounds.
#' @importFrom RTMB qlogis
#' @export
#' 
get_bounds <- function(obj, parameters) {
  
  Lwr <- rep(-Inf, length(obj$par))
  Upr <- rep(Inf, length(obj$par))
  
  Lwr[grep("par_log_psi", names(obj$par))] <- log(0.5)
  Upr[grep("par_log_psi", names(obj$par))] <- log(3)
  # These were the old M bounds
  # Lwr[grep("par_log_m0", names(obj$par))] <- log(0.2)
  # Upr[grep("par_log_m0", names(obj$par))] <- log(0.55)
  # Lwr[grep("par_log_m4", names(obj$par))] <- parameters$par_log_m10
  # Upr[grep("par_log_m4", names(obj$par))] <- log(0.333 * exp(parameters$par_log_m10) + 0.667 * exp(parameters$par_log_m0))
  # Lwr[grep("par_log_m10", names(obj$par))] <- log(0.029)
  # Upr[grep("par_log_m10", names(obj$par))] <- log(0.21)
  # Lwr[grep("par_log_m30", names(obj$par))] <- log(0.2)
  # Upr[grep("par_log_m30", names(obj$par))] <- log(0.7)
  Lwr[grep("par_log_m0", names(obj$par))] <- log(1e-6)
  Upr[grep("par_log_m0", names(obj$par))] <- log(1)
  Lwr[grep("par_log_m4", names(obj$par))] <- log(1e-6)
  Upr[grep("par_log_m4", names(obj$par))] <- log(1)
  Lwr[grep("par_log_m10", names(obj$par))] <- log(1e-6)
  Upr[grep("par_log_m10", names(obj$par))] <- log(1)
  Lwr[grep("par_log_m30", names(obj$par))] <- log(1e-6)
  Upr[grep("par_log_m30", names(obj$par))] <- log(1)
  
  # Lwr[grep("par_log_cpue_tau", names(obj$par))] <- log(0.20)
  # Upr[grep("par_log_cpue_tau", names(obj$par))] <- log(0.20)
  Lwr[grep("par_log_sigma_r", names(obj$par))] <- log(0.1)
  Upr[grep("par_log_sigma_r", names(obj$par))] <- log(2.0)
  Lwr[grep("par_log_h", names(obj$par))] <- log(0.21)
  Upr[grep("par_log_h", names(obj$par))] <- log(1.0)
  # Keep initial F positive; upper cap F <= 3 follows the requested broad bound.
  Lwr[grep("log_init_F_f", names(obj$par))] <- log(1e-12)
  Upr[grep("log_init_F_f", names(obj$par))] <- log(3)
  Lwr[grep("par_rdev_y", names(obj$par))] <- rep(-5, length(parameters$par_rdev_y))
  Upr[grep("par_rdev_y", names(obj$par))] <- rep(5, length(parameters$par_rdev_y))
  Lwr[grep("init_rdev_a", names(obj$par))] <- rep(-5, length(parameters$init_rdev_a))
  Upr[grep("init_rdev_a", names(obj$par))] <- rep(5, length(parameters$init_rdev_a))
  
  check_bounds(opt = obj, lower = Lwr, upper = Upr)
  
  df <- data.frame(parameter = names(obj$par), init = obj$par, lower = Lwr, upper = Upr)
  
  return(df)
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
