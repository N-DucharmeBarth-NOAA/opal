#' Recruitment prior
#'
#' Calculates the prior for the recruitment deviations.
#'
#' @param rdev a \code{vector} of recruitment deviations.
#' @param sigma_r recruitment standard deviation.
#' @return negative log-prior (scalar).
#' @importFrom RTMB dnorm
#' @export
#'
get_recruitment_prior <- function(rdev, sigma_r) {
  "[<-" <- ADoverload("[<-")
  -sum(dnorm(x = rdev, mean = 0, sd = sigma_r, log = TRUE))
}

#' Calculate recruitment
#'
#' Computes recruitment based on Beverton-Holt with depensation and log-normal deviations.
#'
#' @param sbio Spawning biomass.
#' @param rdev Recruitment deviations.
#' @param B0 Unfished biomass.
#' @param alpha,beta Beverton-Holt stock recruitment parameters.
#' @param sigma_r Lognormal SD of recruitment deviations.
#' @param bias_adj Bias adjustment scalar (typically year-specific) applied to
#'   the lognormal correction term.
#' @return Recruitment value (numeric).
#' @export
#' 
get_recruitment <- function(sbio, rdev, B0, alpha, beta, sigma_r = 0.6, bias_adj = 1.0) {
  "[<-" <- ADoverload("[<-")
  rec <- (alpha * sbio) / (beta + sbio) * exp(rdev - bias_adj * 0.5 * sigma_r^2)
  return(rec)
}

#' Calculate Recruitment Bias Adjustment Ramp
#'
#' @param years Numeric vector of model years.
#' @param do_rec_bias_ramp Integer flag (0 = off, 1 = on).
#' @param bias_years Numeric vector of length 4 defining the ramp (ascend
#'   start, plateau start, plateau end, descend end).
#' @param max_bias_adj Numeric scalar for the maximum bias adjustment fraction.
#' @return A numeric vector of length \code{length(years)} containing the bias
#'   adjustment scalars.
#' @export
get_bias_adj_vector <- function(years, do_rec_bias_ramp, bias_years, max_bias_adj) {
  if (is.null(do_rec_bias_ramp) || do_rec_bias_ramp == 0) {
    return(rep(1.0, length(years)))
  }

  bias_adj_y <- rep(0.0, length(years))

  idx_asc <- which(years >= bias_years[1] & years < bias_years[2])
  if (length(idx_asc) > 0) {
    bias_adj_y[idx_asc] <- (years[idx_asc] - bias_years[1]) / (bias_years[2] - bias_years[1])
  }

  idx_full <- which(years >= bias_years[2] & years <= bias_years[3])
  if (length(idx_full) > 0) {
    bias_adj_y[idx_full] <- 1.0
  }

  idx_desc <- which(years > bias_years[3] & years <= bias_years[4])
  if (length(idx_desc) > 0) {
    bias_adj_y[idx_desc] <- 1.0 - ((years[idx_desc] - bias_years[3]) / (bias_years[4] - bias_years[3]))
  }

  bias_adj_y <- bias_adj_y * max_bias_adj

  return(bias_adj_y)
}
