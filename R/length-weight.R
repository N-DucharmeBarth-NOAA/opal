#' Obtain the length at age array
#' 
#' Obtain the mean length (cm) at age for each year and season. This function is 
#' used within \code{get_data} so is generally not needed directly. It simply 
#' structures the input \code{data.frame} into an \code{array}.
#' 
#' The input \code{data.frame} must contain the mean length (cm) at age arranged 
#' into the columns `Year`, `Season`, and ages 0 to 30:
#' 
#' | Year | Season |  0  |    1 | ... |    29 |    30 |
#' |------|--------|-----|------|-----|-------|-------|
#' | 1931 |      1 | 45  | 57.1 | ... |  74.1 |  88.8 |
#' | 1932 |      1 | 45  | 57.1 | ... |  74.1 |  88.8 |
#' | ...  |    ... | ... |  ... | ... |   ... |   ... |
#' | 2021 |      2 | 50  | 68.9 | ... | 184.0 | 184.0 |
#' | 2022 |      2 | 50  | 68.9 | ... | 184.0 | 184.0 |
#' 
#' @param length_mean a \code{data.frame} containing the mean length at age (a) for each year (y) and season (s).
#' @return an \code{array}.
#' @import dplyr
#' @export
#' @examples
#' length_mu_ysa <- get_length_at_age(length_mean = sbt::length_mean)
#' 
get_length_at_age <- function(length_mean) {
  
  years <- sort(unique(length_mean$Year))
  seasons <- sort(unique(length_mean$Season))
  ages <- names(length_mean)[-c(1:2)]

  length_mu_ysa <- array(0, 
                         dim = c(length(years), length(seasons), length(ages)), 
                         dimnames = list(year = years, season = seasons, age = ages))
  
  length_mu_ysa[,1,] <- length_mean %>%
    filter(.data$Season == 1) %>%
    arrange(.data$Year) %>%
    select(-Year, -Season) %>%
    as.matrix()
  
  length_mu_ysa[,2,] <- length_mean %>%
    filter(.data$Season == 2) %>%
    arrange(.data$Year) %>%
    select(-Year, -Season) %>%
    as.matrix()
  
  return(length_mu_ysa)
}

#' Obtain the mean weight at age
#' 
#' Obtain the mean weight (kg) at age for each fishery each year. This function 
#' is used within \code{get_data} so is generally not needed directly.
#' 
#' @param length_mu_ysa an \code{array} containing the mean length at age (a) for each year (y) and season (s). This \code{array} can be generated using the function \code{get_length_at_age}.
#' @param length_sd_a a \code{vector} containing the standard deviation of mean length at age.
#' @return an \code{array}.
#' @importFrom RTMB pnorm
#' @export
#' @examples
#' length_mu_ysa <- get_length_at_age(length_mean = sbt::length_mean)
#' weight_fya <- get_weight_at_age(length_mu_ysa = length_mu_ysa, 
#'                                 length_sd_a = sbt::length_sd$SD)
#' 
get_weight_at_age <- function(length_mu_ysa, length_sd_a) {

  n_year <- dim(length_mu_ysa)[1]
  n_age <- dim(length_mu_ysa)[3]
  
  if (length(length_sd_a) != n_age) stop("Dimension mismatch")
  
  # length-weight parameters
  prma <- prmb <- rep(NA, 4)
  
  prma[1] = 1.15 * 2.942e-06
  prma[2] = 3.13088e-05
  prma[3] = 1.5577e-05
  prma[4] = 1.15 * 1.8241e-07
  
  prmb[1] = 3.3438
  prmb[2] = 2.9058
  prmb[3] = 3.0214
  prmb[4] = 3.9056
  
  step <- 0.5
  
  n_fishery <- 6
  dn <- dimnames(length_mu_ysa)
  wt_age_fya <- array(NA, dim = c(n_fishery, n_year, n_age), dimnames = list(fishery = 1:n_fishery, year = dn$year, age = dn$age))
  
  # itype is the same as iff (fishery) for fisheries 1, 2 and 3, itype = 4 is the Aussie fishery (iff = 6)
  for (itype in 1:4) {
    if (itype <= 2) is <- 2 # season 2
    if (itype >= 3) is <- 1
    for (y in 1:n_year) {
      for (a in 1:n_age) {
        cumhld <- 0
        wt <- 0
        mu_len <- length_mu_ysa[y, is, a]
        len <- mu_len - 4 * length_sd_a[a]
        iprm <- itype
        if (itype == 3 && mu_len > 130) iprm <- 4
        if (itype == 4) iprm <- 2
        while (len < (mu_len + 4 * length_sd_a[a])) {
          cum <- pnorm((len + step - mu_len) / length_sd_a[a])
          dens <- cum - cumhld
          wt <- wt + dens * prma[iprm] * (len + 0.5 * step)^prmb[iprm]
          cumhld <- cum
          len <- len + step
        }
        iff <- itype
        if (itype == 4) iff <- 6 # Aussie fishery
        wt_age_fya[iff, y, a] <- wt
      }
    }
  }
  wt_age_fya[4,,] <- wt_age_fya[3,,] # LL4 = LL3 weights
  wt_age_fya[5,,] <- wt_age_fya[3,,] # Indonessian = LL3 weights
  # spwt_age_ya = wt_age_fya(4) # spawning biomass weight
  # spwt_age_ya /= 1000.;
  
  return(wt_age_fya / 1000)
}

#' Obtain dl
#' 
#' Obtain dl for use later in calculating phi. This function is used within 
#' the [get_data] function so it is generally not needed directly.
#' 
#' Obtain dl for use later in calculating phi(a, y). From the 
#' PRELIMINARY_CALCS_SECTION of \code{sbtmod.tpl}. More fine-scale than 
#' lenage_dist_syal. Needs a bit more detail to get it right.
#' Integrates over length-at-age distribution to get phi(age,year).
#' 
#' @param length_mu_ysa an \code{array} containing the mean length at age (a) for each year (y) and season (s).
#' @param length_sd_a a \code{vector} containing the standard deviation of the mean length at age.
#' @return an \code{array}.
#' @export
#' @examples
#' length_mu_ysa <- get_length_at_age(length_mean = sbt::length_mean)
#' dl_yal <- get_dl(length_mu_ysa = length_mu_ysa, length_sd_a = sbt::length_sd$SD)
#' 
get_dl <- function(length_mu_ysa, length_sd_a) {
  
  n_bins <- 15
  n_year <- dim(length_mu_ysa)[1]
  n_age <- dim(length_mu_ysa)[3]
  
  if (length(length_sd_a) != n_age) stop("Dimension mismatch")
  if (dim(length_mu_ysa)[2] != 2) stop("Input length_mu_ysa season dimension is not 2")

  dl_yal <- array(NA, dim = c(n_year, n_age, n_bins))
  
  for (iy in 1:n_year) { # for(int y=first_yr;y<=last_yr;y++) {
    for (ia in 1:n_age) { # for(int a=0;a<=last_age;a++) {
      lref <- length_mu_ysa[iy, 2, ia] # lref = mean_len_age(2,y,a);
      sdref <- length_sd_a[ia] # sdref = std_len(a);
      llq <- max(0, lref - 1.98 * sdref) # llq = lref-1.98*sdref;
      luq <- lref + 1.98 * sdref # luq = lref+1.98*sdref;
      ldel <- (luq - llq) / (n_bins - 1) # ldel = (luq-llq)/double(nbins2-1);
      dsum <- 0
      for (il in 1:n_bins) { # for(int il=1;il<=nbins2;il++) {
        ltmp <- llq + (il - 1) * ldel # ltmp = llq+double(il-1)*ldel;
        dl_yal[iy, ia, il] <- 1 / sqrt(2 * pi * sdref^2) * exp(-(ltmp - lref) * (ltmp - lref) / (2 * sdref^2)) # dl(y,a,il) = 1./sqrt(2.*M_PI*sdref*sdref) * exp(-(ltmp-lref)*(ltmp-lref)/(2.*sdref*sdref));
        dsum <- dsum + dl_yal[iy, ia, il] # dsum += dl(y,a,il);
      }
      # dl_yal[iy, ia,] <- dl_yal[iy, ia,] / sum(dl_yal[iy, ia,]) # for(int il=1;il<=nbins2;il++) dl(y,a,il) /= dsum;
      for (il in 1:n_bins) dl_yal[iy, ia, il] <- dl_yal[iy, ia, il] / dsum # for(int il=1;il<=nbins2;il++) dl(y,a,il) /= dsum;
    }
  }
  
  dn <- dimnames(length_mu_ysa)
  dimnames(dl_yal) <- list(year = dn$year, age = dn$age, bin = 1:15)
  
  return(dl_yal)
}

#' Plot the distribution of length at age
#' 
#' Plot the distribution of length (cm) at age by year and season.
#' 
#' @param data a \code{list} containing the data that was passed to \code{MakeADFun}.
#' @param years the years to show on the plot.
#' @return a \code{ggplot2} object.
#' @import ggplot2
#' @importFrom ggridges geom_density_ridges
#' @importFrom stats rnorm
#' @importFrom reshape2 melt
#' @export
#' 
plot_dist_length_at_age <- function(data, years = c(1931, 2022)) {
  
  ages <- data$min_age:data$max_age
  yrs <- data$first_yr:(data$last_yr + 1)
  
  df <- melt(data$length_mu_ysa) %>%
    mutate(sd = data$length_sd_a[Var3], Year = yrs[Var1], Season = Var2, Age = ordered(ages[Var3])) %>% 
    filter(Year %in% years) %>%
    slice(rep(1:n(), each = 1000)) %>%
    rowwise() %>%
    mutate(sim = rnorm(n = 1, mean = value, sd = sd), Season = paste0("Season: ", Season))
  
  p <- ggplot(data = df, aes(x = sim, y = Age, color = Age, fill = Age)) +
    geom_density_ridges(alpha = 0.75, rel_min_height = 1e-6) +
    facet_grid(Season ~ Year) +
    labs(x = "Length (cm)", y = "Age") +
    theme(legend.position = "none")
  
  return(p)
}

#' Plot length at age
#' 
#' Plot the mean length (cm) at age by year and season. Plus and minus one 
#' standard deviation is also plotted.
#' 
#' @param data a \code{list} containing the data that was passed to \code{MakeADFun}.
#' @param years the years to show on the plot.
#' @return a \code{ggplot2} object.
#' @import ggplot2
#' @importFrom reshape2 melt
#' @export
#' 
plot_length_at_age <- function(data, years = c(1931, 2022)) {
  
  df_sd <- data.frame(age = as.integer(names(data$length_sd_a)), SD = data$length_sd_a)
  
  df_mu <- melt(data$length_mu_ysa) %>%
    mutate(season = paste("Season:", season)) %>%
    filter(year %in% years) %>%
    left_join(df_sd, by = join_by("age"))
  
  ggplot(data = df_mu, aes(x = age, y = value, color = factor(year), fill = factor(year))) +
    geom_ribbon(aes(ymin = value - SD, ymax = value + SD), alpha = 0.25, color = NA) +
    geom_line() +
    facet_wrap(season ~ .) +
    labs(x = "Age", y = "Length (cm)", color = "Year", fill = "Year") +
    scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05)))
}

#' Plot weight at age
#' 
#' Plot the mean weight (kg) at age by fishery for a subset of years.
#' 
#' @param data a \code{list} containing the data that was passed to \code{MakeADFun}.
#' @param years the years to show on the plot.
#' @return a \code{ggplot2} object.
#' @import ggplot2
#' @importFrom reshape2 melt
#' @export
#' 
plot_weight_at_age <- function(data, years = c(1931, 2000, 2010, 2022)) {
  
  fisheries <- c("LL1", "LL2", "LL3, LL4, Indonesian", "Australian surface")
  
  df <- melt(data$weight_fya[c(1:3, 6),,]) %>%
    mutate(fishery = ifelse(.data$fishery == 6, 4, .data$fishery)) %>%
    mutate(fishery = factor(fisheries[.data$fishery], levels = fisheries)) %>%
    filter(.data$year %in% years)
  
  ggplot(data = df, aes(x = .data$age, y = .data$value * 1000, color = factor(.data$year))) +
    geom_line() +
    facet_wrap(fishery ~ .) +
    labs(x = "Age", y = "Weight (kg)", color = "Year") +
    scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05)))
}
