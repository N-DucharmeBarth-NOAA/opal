#' Estimate temporal autocorrelation in recruitment deviations
#'
#' Calculates the AR1 autocorrelation coefficient (phi) for recruitment deviations.
#'
#' @param first_yr First model year.
#' @param last_yr Last model year.
#' @param rdev Vector of recruitment deviations.
#' @return Estimated autocorrelation.
#' @export
#' @examples
#' \dontrun{
#'   first_yr <- 1931
#'   last_yr <- 2022
#'   N <- length(first_yr:last_yr)
#'   rdev <- arima.sim(list(order = c(1, 0, 0), ar = 0.5), n = N)
#'   get_rho(first_yr, last_yr, rdev)
#' }
#' 
get_rho <- function(first_yr = 1931, last_yr = 2022, rdev) {
  "[<-" <- ADoverload("[<-")
  # Model years: 1931-2022; Rec years: 1932-2023; n_years = n_recs: 92
  i1 <- 8 - first_yr # int i1 = 1965 but don't add 1 because years are offset (see above)
  i2 <- last_yr - first_yr - 5
  t1 <- rdev[i1:(i2 - 1)]
  t2 <- rdev[(i1 + 1):i2]
  t1m <- mean(t1)
  t2m <- mean(t2)
  phi <- sum((t1 - t1m) * (t2 - t2m)) / (sqrt(sum((t1 - t1m)^2)) * sqrt(sum((t2 - t2m)^2)))
  # phi <- cor(t1, t2) # same as above
  return(phi)
}

#' Recruitment prior
#'
#' Calculates the prior for the recruitment deviations.
#'
#' @param rdev a \code{vector} of recruitment deviations.
#' @param sigma_r recruitment standard deviation.
#' @return negative log-prior (scalar).
#' @importFrom RTMB dnorm dautoreg
#' @export
#'
get_recruitment_prior <- function(rdev, sigma_r) {
  "[<-" <- ADoverload("[<-")
  n_year <- length(rdev)
  # tau_ac2 <- get_rho(first_yr, last_yr, rdev)
  # tau_ac2 <- 0.05 # phi
  # r1 <- rdev[1:(n_year - 3)]
  # r2 <- rdev[(n_year - 2):n_year]
  lp1 <- -sum(dnorm(x = rdev, mean = 0, sd = sigma_r, log = TRUE))
  lp2 <- 0
  # lp <- n_year * log(sigma) + 0.5 * sum(r1^2) / sigma^2 + 0.5 * sum(r2^2) / (sigma^2 * (1 - phi^2))
  # lp1 <- -sum(dnorm(x = r1, mean = 0, sd = sigma, log = TRUE))
  # lp2 <- -dautoreg(x = r2, phi = phi, log = TRUE, scale = sigma)
  # lp2 <- -dautoreg(x = r2, phi = phi, log = TRUE, scale = sigma / sqrt(1 - phi^2))
  lp <- lp1 + lp2
  # REPORT(tau_ac2)
  return(lp)
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
