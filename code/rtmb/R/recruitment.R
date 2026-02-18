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
#' first_yr <- 1931
#' last_yr <- 2022
#' N <- length(first_yr:last_yr)
#' rdev <- arima.sim(list(order = c(1, 0, 0), ar = 0.5), n = N)
#' get_rho(first_yr, last_yr, rdev)
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
  tau_ac2 <- 0.05 # phi
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
#' @return Recruitment value (numeric).
#' @export
#' 
get_recruitment <- function(sbio, rdev, B0, alpha, beta, sigma_r = 0.6) {
  "[<-" <- ADoverload("[<-")
  rec <- (alpha * sbio) / (beta + sbio) * exp(rdev - 0.5 * sigma_r^2)
  return(rec)
}
