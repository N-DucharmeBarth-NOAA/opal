#' Estimate temporal autocorrelation in recruitment deviations
#'
#' Calculates the AR1 autocorrelation coefficient (phi) for recruitment deviations.
#'
#' @param first_yr First model year.
#' @param last_yr Last model year.
#' @param rdev_y Vector of recruitment deviations.
#' @return Estimated autocorrelation.
#' @export
#' @examples
#' first_yr <- 1931
#' last_yr <- 2022
#' N <- length(first_yr:last_yr)
#' rdev_y <- arima.sim(list(order = c(1, 0, 0), ar = 0.5), n = N)
#' get_rho(first_yr, last_yr, rdev_y)
#' 
get_rho <- function(first_yr = 1931, last_yr = 2022, rdev_y) {
  "[<-" <- ADoverload("[<-")
  # Model years: 1931-2022; Rec years: 1932-2023; n_years = n_recs: 92
  i1 <- 8 - first_yr # int i1 = 1965 but don't add 1 because years are offset (see above)
  i2 <- last_yr - first_yr - 5
  t1 <- rdev_y[i1:(i2 - 1)]
  t2 <- rdev_y[(i1 + 1):i2]
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
#' @param sigma recruitment standard deviation.
#' @param phi temporal autocorrelation squared.
#' @return negative log-prior (scalar).
#' @importFrom RTMB dnorm dautoreg
#' @export
#'
get_recruitment_prior <- function(rdev, sigma_r, phi) {
  "[<-" <- ADoverload("[<-")
  n_year <- length(rdev)
  # tau_ac2 <- get_rho(first_yr, last_yr, rdev)
  tau_ac2 <- 0.05
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
#' @param sr_dep Depensation parameter (default 1e-10).
#' @return Recruitment value (numeric).
#' @export
#' 
get_recruitment <- function(sbio, rdev, B0, alpha, beta, sigma_r = 0.6, sr_dep = 0) {
  "[<-" <- ADoverload("[<-")
  rec <- (alpha * sbio) / (beta + sbio) * (1 - exp(log(0.5) * sbio / (sr_dep * B0))) * exp(rdev - 0.5 * sigma_r^2)
  # ADREPORT(sigma_r)
  return(rec)
}

#' Plot recruitment deviates
#' 
#' Plot recruitment deviates by year.
#' 
#' @param data a \code{list} containing the data that was passed to \code{MakeADFun}.
#' @param object a \code{list} specifying the AD object created using \code{MakeADFun}.
#' @param proj a \code{list} containing the data that was passed to \code{MakeADFun}.
#' @param posterior an \code{rstan} objected created using the \code{tmbstan} function.
#' @param probs a numeric vector of probabilities with values in \code{[0,1]} for plotting quantiles of the posterior distribution.
#' @return a \code{ggplot2} object.
#' @import ggplot2
#' @importFrom stats median quantile
#' @importFrom reshape2 melt
#' @importFrom utils txtProgressBar
#' @export
#' 
plot_rec_devs <- function(data, object, proj = NULL, posterior = NULL, probs = c(0.025, 0.975)) {
  years <- data$first_yr:data$last_yr
  # num <- object$report()$par_rdev_y
  num <- object$env$last.par.best[names(object$par) %in% "par_rdev_y"]
  df <- data.frame(year = years, value = num)
  p <- ggplot(data = df, aes(x = .data$year, y = .data$value, color = "Reconstruction")) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_line() +
    labs(x = "Year", y = "Recruitment deviate", color = NULL)
  if (!is.null(proj)) {
    sim <- melt(proj$rdev_y) %>% filter(!.data$year %in% years)
    p <- p + geom_line(data = sim, aes(group = .data$iteration, color = "Projection"), alpha = 0.75)
  }
  return(p)
}

#' Plot recruitment
#' 
#' Plot recruitment by year.
#' 
#' @param data a \code{list} containing the data that was passed to \code{MakeADFun}.
#' @param object a \code{list} specifying the AD object created using \code{MakeADFun}.
#' @param posterior an \code{rstan} objected created using the \code{tmbstan} function.
#' @param probs a numeric vector of probabilities with values in \code{[0,1]} for plotting quantiles of the posterior distribution.
#' @return a \code{ggplot2} object.
#' @import ggplot2
#' @importFrom stats median quantile
#' @importFrom reshape2 melt
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @export
#' 
plot_recruitment <- function(data, object, posterior = NULL, probs = c(0.025, 0.975)) {
  
  rep <- object$report(object$env$last.par.best)
  years <- data$first_yr:(data$last_yr + 1)
  num <- rep$number_ysa
  N_rec1 <- data.frame(year = years, value = num[,1,1])
  R0 <- rep$R0
  
  p <- ggplot(data = N_rec1, aes(x = .data$year, y = .data$value / 1e6)) +
    geom_hline(yintercept = R0 / 1e6, linetype = "dashed") +
    geom_line(data = N_rec1, linetype = "dashed") +
    labs(x = "Year", y = "Recruitment (millions)") +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05)))
  
  if (!is.null(posterior)) {
    x <- extract_samples(fit = posterior)
    niter <- nrow(x)
    r1 <- obj$report(par = as.numeric(x[1,]))$recruitment_y
    rec <- matrix(NA, nrow = niter, ncol = length(r1))
    pb <- txtProgressBar(min = 1, max = niter, style = 3)
    for (i in 1:niter) {
      rec[i,] <- obj$report(par = as.numeric(x[i,]))$recruitment_y
      setTxtProgressBar(pb, i)
    }
    dimnames(rec) <- list(iteration = 1:niter, year = years)
    rec_mcmc <- melt(rec)
    p <- p + 
      stat_summary(data = rec_mcmc, geom = "ribbon", alpha = 0.5, 
                   # aes(fill = factor(mc_grid)),
                   fun.min = function(x) quantile(x, probs = probs[1]),
                   fun.max = function(x) quantile(x, probs = probs[2])) +
      stat_summary(data = rec_mcmc, geom = "line", fun = median)
  }
  
  return(p)
}
