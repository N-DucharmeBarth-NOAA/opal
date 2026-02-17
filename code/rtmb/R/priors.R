#' Get priors
#' 
#' Don't include priors for recruitment deviates (par_rdev_y) or selectivity 
#' (e.g., par_log_sel_1) here because they are dealt in the [get_recruitment_prior] 
#' and [get_selectivity_prior] functions.
#' 
#' @param parameters A \code{list} specifying the parameters to be passed to \code{MakeADFun}. Can be generated using the [get_parameters] function.
#' @return A \code{list} of priors.
#' @export
#' 
get_priors <- function(parameters) {
  priors <- list()
  priors[["log_B0"]] <- list(type = "normal", par1 = 0, par2 = 1.5, index = which("log_B0" == names(parameters)))
  # priors[["par_log_m0"]] <- list(type = "normal", par1 = 0, par2 = 1.5, index = which("par_log_m0" == names(parameters)))
  # priors[["par_log_m4"]] <- list(type = "normal", par1 = log(0.12), par2 = 0.4, index = which("par_log_m4" == names(parameters)))
  # priors[["par_log_m10"]] <- list(type = "normal", par1 = log(0.1), par2 = 0.06, index = which("par_log_m10" == names(parameters)))
  # priors[["par_log_m30"]] <- list(type = "normal", par1 = log(2), par2 = 1.655705, index = which("par_log_m30" == names(parameters)))
  # priors[["par_log_h"]] <- list(type = "normal", par1 = log(1), par2 = 1.5, index = which("par_log_h" == names(parameters)))
  # priors[["par_log_sigma_r"]] <- list(type = "normal", par1 = 0, par2 = 1.5, index = which("par_log_sigma_r" == names(parameters)))
  priors[["log_cpue_q"]] <- list(type = "normal", par1 = 0, par2 = 1.5, index = which("log_cpue_q" == names(parameters)))
  # priors[["par_cpue_creep"]] <- list(type = "normal", par1 = 0, par2 = 1.5, index = which("par_cpue_creep" == names(parameters)))
  # priors[["par_log_cpue_sigma"]] <- list(type = "normal", par1 = 0, par2 = 1.5, index = which("par_log_cpue_sigma" == names(parameters)))
  # priors[["par_log_cpue_omega"]] <- list(type = "normal", par1 = 0.875, par2 = 0.1, index = which("par_log_cpue_omega" == names(parameters)))
  # priors[["par_log_af_alpha"]] <- list(type = "normal", par1 = 0, par2 = 1.5, index = which("par_log_af_alpha" == names(parameters)))
  # priors[["par_log_lf_alpha"]] <- list(type = "normal", par1 = 0, par2 = 1.5, index = which("par_log_lf_alpha" == names(parameters)))
  # priors[["par_sel_rho_y"]] <- list(type = "normal", par1 = 0, par2 = 1.5, index = which("par_sel_rho_y" == names(parameters)))
  # priors[["par_sel_rho_a"]] <- list(type = "normal", par1 = 0, par2 = 1.5, index = which("par_sel_rho_a" == names(parameters)))
  # priors[["par_log_sel_sigma"]] <- list(type = "normal", par1 = 0, par2 = 1.5, index = which("par_log_sel_sigma" == names(parameters)))
  
  # Broad normal priors for selectivity parameters (all on real line)
  # par_sel is a matrix [n_fishery, 6] — treat as a single block
  # Normal(0, 2) is vague: allows peak to shift ~2 SD from mean length,
  # widths to vary by exp(±2) ≈ 0.14x to 7.4x the length SD, etc.
  if ("par_sel" %in% names(parameters)) {
    priors[["par_sel"]] <- list(type = "normal", par1 = 0, par2 = 2, index = which("par_sel" == names(parameters)))
  }
  
  evaluate_priors(parameters, priors)
  return(priors)
}

#' Evaluate priors
#' 
#' @param parameters A \code{list} specifying the parameters to be passed to \code{MakeADFun}. Can be generated using the [get_parameters] function.
#' @param priors A \code{list} of named \code{list}s specifying priors for the parameters. Can be generated using the [get_priors] function.
#' @return A \code{numeric} value.
#' @importFrom RTMB dnorm dlnorm
#' @importFrom RTMBdist dbeta2 dt2
#' @export
#' @examples
#' parameters <- list(par_log_m4 = log(0.167))
#' priors <- list(
#'   par_log_m4 = list(type = "normal", par1 = log(0.12), par2 = 0.4, 
#'                     index = which("par_log_m4" == names(parameters)))
#' )
#' evaluate_priors(parameters, priors)
#' 
evaluate_priors <- function(parameters, priors) {
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  n <- length(priors)
  lp <- numeric(n)
  for (i in 1:n) {
    type <- priors[[i]]$type
    par1 <- priors[[i]]$par1
    par2 <- priors[[i]]$par2
    idx <- priors[[i]]$index
    x <- parameters[[idx]]
    if (type == "normal") lp[i] <- sum(dnorm(x = x, mean = par1, sd = par2, log = TRUE))
    if (type == "student") lp[i] <- sum(dt2(x = x, mu = par1, sigma = par2, df = 3, log = TRUE))
    # if (type == "cauchy") lp[i] <- sum(dcauchy(x = x, location = par1, scale = par2, log = TRUE))
    if (type == "lognormal") lp[i] <- sum(dlnorm(x = x, meanlog = par1, sdlog = par2, log = TRUE))
    if (type == "beta") lp[i] <- sum(dbeta2(x = x, mu = par1, phi = par2, log = TRUE))
  }
  return(-sum(lp))
}
