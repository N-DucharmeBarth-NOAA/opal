#' Compute selectivity-at-age by year (legacy SBT version)
#'
#' Constructs the selectivity-at-age array by fishery and year, incorporating initial values and changes.
#' This is the legacy version used for SBT model parameter loading in get_parameters().
#'
#' @param n_age Total number of ages.
#' @param max_age Maximum model age.
#' @param first_yr Model start year.
#' @param first_yr_catch Vector of first catch years per fishery.
#' @param sel_min_age_f,sel_max_age_f,sel_end_f Logical vector indicating whether to extend the final selectivity across remaining ages.
#' @param sel_change_year_fy Matrix fishery, year indicating change years.
#' @param par_log_sel Vector of changes in selectivity (log-space).
#' @return 3D array fishery, year, age of selectivity values.
#' @export
#' 
get_selectivity_v1 <- function(n_age, max_age, first_yr, first_yr_catch, 
                            sel_min_age_f, sel_max_age_f, sel_end_f, sel_change_year_fy,
                            par_log_sel) {
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  n_sel <- nrow(sel_change_year_fy)
  n_year <- ncol(sel_change_year_fy)
  ymin <- first_yr_catch - first_yr + 1
  sel_fya <- array(0, dim = c(n_sel, n_year, n_age))
  for (f in seq_len(n_sel)) {
    amin <- sel_min_age_f[f] + 1
    amax <- sel_max_age_f[f] + 1
    ipar <- 1
    for (y in ymin:n_year) {
      if (sel_change_year_fy[f, y] != 0) {
        sel_tmp <- exp(par_log_sel[[f]][ipar,])
        ipar <- ipar + 1
        sel_fya[f, y, amin:amax] <- sel_tmp / mean(sel_tmp)
        if (as.logical(sel_end_f[f]) && amax < max_age) {
          for (a in (amax + 1):n_age) {
            sel_fya[f, y, a] <- sel_fya[f, y, amax]
          }
        }
      } else {
        sel_fya[f, y, ] <- sel_fya[f, y - 1, ]
      }
    }
  }
  return(sel_fya)
}

#' Plot selectivity
#' 
#' Plot selectivity by fishery, year, and age.
#' 
#' @param data a \code{list} containing the data that was passed to \code{MakeADFun}.
#' @param object a \code{list} specifying the AD object created using \code{MakeADFun}.
#' @param posterior an \code{rstan} objected created using the \code{tmbstan} function.
#' @param probs a numeric vector of probabilities with values in \code{[0, 1]} for plotting quantiles of the posterior distribution.
#' @param years the years to show on the plot.
#' @param fisheries the fisheries to show on the plot.
#' @param ... options passed on to \code{geom_density_ridges}.
#' @return a \code{ggplot2} object.
#' @import ggplot2
#' @import dplyr
#' @importFrom rlang .data
#' @importFrom stats median quantile
#' @importFrom reshape2 melt
#' @importFrom ggridges geom_density_ridges
#' @importFrom scales pretty_breaks
#' @export
#' 
plot_selectivity <- function(data, object, posterior = NULL, probs = c(0.025, 0.975), 
                             years = 2013:2022, 
                             fisheries = c("LL1", "LL2", "LL3", "LL4", "Indonesian", "Australian", "CPUE"), ...) {
  
  ages <- data$min_age:data$max_age
  yrs <- data$first_yr:data$last_yr
  fsh <- c("LL1", "LL2", "LL3", "LL4", "Indonesian", "Australian", "CPUE")
  rem <- c(data$removal_switch_f, 0)
  fyr <- c(data$first_yr_catch_f, 1969)
  # data$catch_obs_ysf[,,3]
  # data$catch_obs_ysf[,,4]
  # cut out years with no catch, drop LF crosses in these years too?
  
  # Range of estimated ages for selectivity
  
  range_ages <- data.frame(fishery = fsh, min = data$sel_min_age_f, max = data$sel_max_age_f, removal = rem) %>% 
    filter(fishery %in% fisheries) %>%
    filter(removal == 0) %>%
    mutate(fishery = factor(fishery, levels = fsh))
  
  # Selectivity change years
  
  dimnames(data$sel_change_year_fy) <- NULL
  df_change <- melt(data$sel_change_year_fy) %>% 
    mutate(fishery = fsh[Var1], year = yrs[Var2]) %>%
    select(-Var1, -Var2) %>%
    pivot_wider(names_from = fishery) %>%
    mutate_at(fsh, function(x){1 + cumsum(x > 0)}) %>%
    pivot_longer(LL1:Australian, names_to = "fishery", values_to = "change") %>%
    mutate(change = ifelse(change %% 2, 1, 2))
  
  # if (is.null(posterior)) {
  # df_sel <- get_array(object$report()$sel_fya) %>%
  df1 <- object$report()$sel_fya
  # df2 <- proj
  
  df_sel <- df1 %>%
    melt() %>%
    mutate(fishery = fsh[Var1], removal = rem[Var1], 
           first_yr_catch = fyr[Var1], year = yrs[Var2], 
           age = ages[Var3]) %>%
    left_join(df_change, by = join_by("year", "fishery")) %>% 
    filter(fishery %in% fisheries, 
           year %in% years, 
           year >= first_yr_catch,
           removal == 0)
  
  # df_hrate <- get_array(object$report()$hrate_fya) %>%
  # df_hrate <- object$report()$hrate_fya %>%
  #   melt() %>%
  #   mutate(fishery = fsh[Var1], removal = rem[Var1], 
  #          first_yr_catch = fyr[Var1], year = yrs[Var2], 
  #          age = ages[Var3]) %>%
  #   filter(fishery %in% fisheries, year %in% years, 
  #          year >= first_yr_catch, removal == 1) %>%
  #   group_by(fishery, year) %>%
  #   mutate(value = value / sum(value, na.rm = TRUE)) %>%
  #   mutate(change = ifelse(year %% 2, 1, 2))
  # } else {
  #   # df0 <- get_posterior(object = object, posterior = posterior, pars = "sel_fya", iters = 1)
  #   post <- extract(object = posterior, pars = "lp__", permuted = FALSE, include = FALSE)
  #   chains <- dim(post)[2]
  #   iters <- nrow(post)
  #   # iters <- 5
  #   r1 <- object$report(par = post[1, 1, ])
  #   sel_jifya <- array(data = NA, dim = c(chains, iters, dim(get_array(r1$sel_fya))))
  #   hrate_jifya <- array(data = NA, dim = c(chains, iters, dim(get_array(r1$hrate_fya))))
  #   for (j in 1:chains) {
  #     for (i in 1:iters) {
  #       r1 <- object$report(par = post[i, j, ])
  #       sel_jifya[j,i,,,] <- get_array(r1$sel_fya)
  #       hrate_jifya[j,i,,,] <- get_array(r1$hrate_fya)
  #     }
  #   }
  #   
  #   df_sel <- sel_jifya %>%
  #     melt() %>%
  #     rename(chain = Var1, iter = Var2) %>%
  #     mutate(fishery = fsh[Var3], removal = rem[Var3], year = yrs[Var4], age = ages[Var5]) %>%
  #     group_by(fishery, removal, year, age) %>%
  #     mutate(value = median(value)) %>%
  #     left_join(df_change, by = join_by("year", "fishery")) %>% 
  #     filter(fishery %in% fisheries, year %in% years, year >= data$first_yr_catch) %>%
  #     filter(removal == 0) %>%
  #     mutate(fishery = factor(fishery, levels = fsh))
  #   
  #   df_hrate <- hrate_jifya %>%
  #     melt() %>%
  #     rename(chain = Var1, iter = Var2) %>%
  #     mutate(fishery = fsh[Var3], removal = rem[Var3], year = yrs[Var4], age = ages[Var5]) %>%
  #     group_by(fishery, removal, year, age) %>%
  #     mutate(value = median(value)) %>%
  #     filter(fishery %in% fisheries, year %in% years, year >= data$first_yr_catch) %>%
  #     filter(removal == 1) %>%
  #     mutate(fishery = factor(fishery, levels = fsh)) %>%
  #     group_by(fishery, year) %>%
  #     mutate(value = value / sum(value, na.rm = TRUE)) %>%
  #     mutate(change = ifelse(year %% 2, 1, 2))
  # }
  
  # df <- bind_rows(df_sel, df_hrate) %>%
  df <- df_sel %>%
    mutate(fishery = factor(fishery, levels = fsh))
  
  # Age and length composition specifications
  
  specs_lf <- data.frame(year = data$lf_year + data$first_yr, fishery = fsh[data$lf_fishery])
  specs_af <- data.frame(year = data$af_year + data$first_yr, fishery = fsh[data$af_fishery])
  specs <- bind_rows(specs_lf, specs_af) %>%
    inner_join(df, by = join_by(year, fishery)) %>%
    select(year, fishery) %>%
    distinct() %>%
    mutate(value = NA) %>%
    mutate(fishery = factor(fishery, levels = fsh))
  
  ggplot(data = df, aes(x = age, y = year, height = value, group = year)) +
    geom_vline(data = range_ages, aes(xintercept = min), linetype = "dashed") +
    geom_vline(data = range_ages, aes(xintercept = max), linetype = "dashed") +
    geom_density_ridges(aes(fill = factor(change), colour = factor(change)), stat = "identity", alpha = 0.5, rel_min_height = 0) +
    geom_point(data = specs, aes(x = 0.25, y = year), shape = 4) +
    facet_wrap(fishery ~ .) +
    labs(x = "Age", y = "Year") +
    theme(legend.position = "none") +
    scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
    scale_y_reverse(breaks = pretty_breaks())
}

#' Logistic selectivity as a function of length
#'
#' @param len Numeric vector of length-bin midpoints.
#' @param par Numeric vector of length 6 containing selectivity parameters.
#'   Only the first two entries are used: par[1] (a) and par[2] (b).
#'   `a` is the inflection point on the real line, transformed via
#'   `mean(len) + a * sd(len)`, so `a = 0` gives inflection at `mean(len)`.
#'   `b` is log-scale 95% width, transformed via `exp(b) * sd(len)`,
#'   so `b = 0` gives a 95% width equal to `sd(len)`.
#' @return Numeric vector of selectivity values in (0, 1].
#' @export
#'
sel_logistic <- function(len, par) {
  mu  <- mean(len)
  sd  <- sd(len)
  a   <- par[1]
  b   <- par[2]
  a50 <- mu + a * sd
  w95 <- exp(b) * sd
  sel <- 1 / (1 + exp(-log(19) * (len - a50) / w95))
  return(sel)
}

#' Double-normal selectivity as a function of length (SS3 pattern 24, full form)
#'
#' All parameters are on the real line, making this parameterization suitable
#' for unconstrained optimization and AD. The length-bin vector `x` is used
#' to define natural scales for position and width parameters.
#'
#' Only the full form is implemented (e and f both active). The simplified
#' forms where e <= -999 or f <= -999 are not supported.
#'
#' @param x Numeric vector of length-bin midpoints.
#' @param par Numeric vector of length 6 containing selectivity parameters:
#'   par[1] (a): Peak location (real line). Transformed via
#'     `mean(x) + a * sd(x)`, so `a = 0` places the peak at `mean(x)`.
#'   par[2] (b): Plateau width (real line). Controls the distance from peak to the
#'     start of the descending limb via logistic transform of the available
#'     range: `peak + bin_width + (0.99 * max(x) - peak - bin_width) / (1 + exp(-b))`.
#'   par[3] (c): Ascending width (real line, log-space). Actual width = `exp(c) * sd(x)`.
#'     `c = 0` gives an ascending width equal to `sd(x)`.
#'   par[4] (d): Descending width (real line, log-space). Actual width = `exp(d) * sd(x)`.
#'     `d = 0` gives a descending width equal to `sd(x)`.
#'   par[5] (e): Initial selectivity (real line, logit-space). Transformed via
#'     `1 / (1 + exp(-e))`, so `e = 0` gives initial selectivity of 0.5,
#'     large negative values give ~0, large positive values give ~1.
#'   par[6] (f): Final selectivity (real line, logit-space). Same transform as `e`.
#' @return Numeric vector of selectivity values in [0, 1].
#' @export
#'
sel_double_normal <- function(x, par) {
  mu <- mean(x)
  sd <- sd(x)
  # Extract parameters from vector
  a <- par[1]
  b <- par[2]
  c <- par[3]
  d <- par[4]
  e <- par[5]
  f <- par[6]
  # --- Transform parameters from real line to natural scale ---
  peak       <- mu + a * sd
  upselex    <- exp(c) * sd
  downselex  <- exp(d) * sd
  point1     <- 1 / (1 + exp(-e))
  point2     <- 1 / (1 + exp(-f))

  # --- Derived quantities ---
  startbin  <- 1
  j1        <- startbin - 1
  j2        <- length(x)
  bin_width <- x[2] - x[1]

  peak2  <- peak + bin_width + (0.99 * x[j2] - peak - bin_width) / (1 + exp(-b))

  t1min  <- exp(-(x[startbin] - peak)^2 / upselex)
  t2min  <- exp(-(x[j2] - peak2)^2 / downselex)

  # --- Compute selectivity ---
  t1 <- x - peak
  t2 <- x - peak2

  join1 <- 1 / (1 + exp(-(20 / (1 + abs(t1))) * t1))
  join2 <- 1 / (1 + exp(-(20 / (1 + abs(t2))) * t2))

  asc <- point1 + (1 - point1) * (exp(-t1^2 / upselex) - t1min) / (1 - t1min)
  dsc <- 1 + (point2 - 1) * (exp(-t2^2 / downselex) - 1) / (t2min - 1)

  sel <- rep(NA_real_, length(x))
  idx <- (j1 + 1):j2
  sel[idx] <- asc[idx] * (1 - join1[idx]) +
    join1[idx] * (1 - join2[idx] + dsc[idx] * join2[idx])

  if (startbin > 1) {
    sel[1:startbin] <- (x[1:startbin] / x[startbin])^2 * sel[startbin]
  }
  if (j2 < length(x)) {
    sel[(j2 + 1):length(x)] <- sel[j2]
  }

  return(sel)
}

#' Probability of length at age matrix (age-length key)
#'
#' Computes the probability distribution of fish lengths across different
#' age classes. AD-compatible: uses RTMB::pnorm and ADoverload so that
#' gradients propagate if growth parameters (mu_a, sd_a) are estimated.
#'
#' @param len_lower Numeric vector of lower bounds of length bins (length L). Data only.
#' @param len_upper Numeric vector of upper bounds of length bins (length L). Data only.
#' @param mu_a Numeric vector of mean length at age (length A). May be AD.
#' @param sd_a Numeric vector of SD of length at age (length A). May be AD.
#' @return Matrix of dimensions L x A where columns sum to 1.
#' @importFrom RTMB ADoverload pnorm
#' @export
#'
get_pla <- function(len_lower, len_upper, mu_a, sd_a) {
  "[<-" <- ADoverload("[<-")
  n_len <- length(len_lower)
  n_age <- length(mu_a)
  pla <- matrix(0, nrow = n_len, ncol = n_age)
  for (a in seq_len(n_age)) {
    for (z in seq_len(n_len)) {
      pla[z, a] <- pnorm((len_upper[z] - mu_a[a]) / sd_a[a]) -
                   pnorm((len_lower[z] - mu_a[a]) / sd_a[a])
    }
    col_sum <- sum(pla[, a])
    # Normalize — use epsilon instead of if-branch for AD safety
    pla[, a] <- pla[, a] / (col_sum + 1e-12)
  }
  return(pla)
}

#' Compute selectivity-at-age from length-based selectivity curves
#'
#' Defines selectivity as a parametric function of length (logistic or
#' double-normal), then converts to selectivity-at-age by matrix-multiplying
#' with the probability-of-length-at-age (PLA/ALK).
#'
#' The PLA is computed internally from \code{mu_a} and \code{sd_a} so that
#' AD gradients propagate correctly if growth parameters are estimated in
#' the future.
#'
#' Selectivity is currently time-invariant within each fishery (constant
#' across years).
#'
#' @param data A list containing model data. Required elements:
#'   n_fishery, n_year, n_age, sel_type_f, sel_lengths, sel_len_lower,
#'   sel_len_upper.
#' @param par_sel Numeric matrix of dimensions [n_fishery, 6]. Each row is
#'   a real-line parameter vector. For logistic (sel_type_f == 1), only
#'   columns 1:2 are used. For double-normal (sel_type_f == 2), all 6 are used.
#' @param mu_a Numeric vector (length n_age) of mean length at age.
#'   May be AD if growth parameters are estimated.
#' @param sd_a Numeric vector (length n_age) of SD of length at age.
#'   May be AD if growth parameters are estimated.
#' @return 3D array sel_fya of dimensions [n_fishery, n_year, n_age].
#' @importFrom RTMB ADoverload
#' @export
#'
get_selectivity <- function(data, par_sel, mu_a, sd_a) {
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")

  n_fishery <- data$n_fishery
  n_year <- data$n_year
  n_age <- data$n_age

  # Compute PLA inside the function (AD-safe)
  pla <- get_pla(data$sel_len_lower, data$sel_len_upper, mu_a, sd_a)

  sel_fya <- array(0, dim = c(n_fishery, n_year, n_age))

  for (f in seq_len(n_fishery)) {
    par_f <- par_sel[f, ]

    # Branch on data value (not AD) — safe for AD
    if (data$sel_type_f[f] == 1L) {
      sel_at_length <- sel_logistic(data$sel_lengths, par_f)
    } else {
      sel_at_length <- sel_double_normal(data$sel_lengths, par_f)
    }

    # Convert to selectivity-at-age via PLA
    sel_at_age <- as.vector(t(pla) %*% sel_at_length)

    # Normalize to max = 1
    sel_at_age <- sel_at_age / max(sel_at_age)

    # Time-invariant: replicate across years
    for (y in seq_len(n_year)) {
      sel_fya[f, y, ] <- sel_at_age
    }
  }

  return(sel_fya)
}

#' Convert SS3 selectivity parameters to RTMB real-line parameterization
#'
#' Takes SS3 natural-scale selectivity parameters and converts them to the
#' centered real-line parameterization used by the RTMB \code{sel_logistic}
#' and \code{sel_double_normal} functions.
#'
#' @param ss3_pars A data.frame or matrix with one row per fishery. For
#'   double-normal (pattern 24): columns are peak, top_logit, ascend_se,
#'   descend_se, start_logit, end_logit. For logistic (pattern 1): columns
#'   are inflection, width.
#' @param sel_type_f Integer vector (length n_fishery). 1 = logistic, 2 = double-normal.
#' @param sel_lengths Numeric vector of selectivity length-bin midpoints
#'   (same vector that will be passed to sel_logistic/sel_double_normal).
#' @return Numeric matrix [n_fishery, 6] of RTMB real-line parameters.
#' @export
#'
convert_ss3_selex_to_rtmb <- function(ss3_pars, sel_type_f, sel_lengths) {
  mu_len <- mean(sel_lengths)
  sd_len <- sd(sel_lengths)
  n_fishery <- length(sel_type_f)
  par_sel <- matrix(0, nrow = n_fishery, ncol = 6)

  for (f in seq_len(n_fishery)) {
    if (sel_type_f[f] == 1L) {
      # Logistic (SS3 pattern 1): inflection (cm), width (cm)
      ss3_inflection <- ss3_pars[f, 1]
      ss3_width <- ss3_pars[f, 2]
      par_sel[f, 1] <- (ss3_inflection - mu_len) / sd_len
      par_sel[f, 2] <- log(ss3_width / sd_len)
      # par_sel[f, 3:6] remain 0 (unused)
    } else {
      # Double-normal (SS3 pattern 24): peak, top_logit, ascend_se, descend_se, start_logit, end_logit
      ss3_peak        <- ss3_pars[f, 1]
      ss3_top_logit   <- ss3_pars[f, 2]
      ss3_ascend_se   <- ss3_pars[f, 3]
      ss3_descend_se  <- ss3_pars[f, 4]
      ss3_start_logit <- ss3_pars[f, 5]
      ss3_end_logit   <- ss3_pars[f, 6]

      par_sel[f, 1] <- (ss3_peak - mu_len) / sd_len        # a: peak location
      par_sel[f, 2] <- ss3_top_logit                         # b: plateau (same space)
      par_sel[f, 3] <- ss3_ascend_se - log(sd_len)           # c: ascending width
      par_sel[f, 4] <- ss3_descend_se - log(sd_len)          # d: descending width

      # Handle start_logit = -999 (SS3 convention for "fix at 0")
      if (ss3_start_logit <= -999) {
        par_sel[f, 5] <- -9
      } else {
        par_sel[f, 5] <- ss3_start_logit                     # e: initial selectivity
      }

      # Handle end_logit = -999
      if (ss3_end_logit <= -999) {
        par_sel[f, 6] <- -9
      } else {
        par_sel[f, 6] <- ss3_end_logit                       # f: final selectivity
      }
    }
  }
  return(par_sel)
}

#' Convert RTMB real-line selectivity parameters back to SS3 natural scale
#'
#' Inverse of \code{convert_ss3_selex_to_rtmb}. Useful for verifying
#' round-trip conversion and reporting parameter values in natural units.
#'
#' @param par_sel Numeric matrix [n_fishery, 6] of RTMB real-line parameters.
#' @param sel_type_f Integer vector. 1 = logistic, 2 = double-normal.
#' @param sel_lengths Numeric vector of selectivity length-bin midpoints.
#' @return Data.frame with SS3-scale parameter values per fishery.
#' @export
#'
convert_rtmb_selex_to_ss3 <- function(par_sel, sel_type_f, sel_lengths) {
  mu_len <- mean(sel_lengths)
  sd_len <- sd(sel_lengths)
  n_fishery <- nrow(par_sel)
  result <- data.frame(
    fleet = 1:n_fishery,
    type = sel_type_f,
    peak_or_inflection = NA_real_,
    top_logit_or_width = NA_real_,
    ascend_se = NA_real_,
    descend_se = NA_real_,
    start_logit = NA_real_,
    end_logit = NA_real_
  )
  for (f in seq_len(n_fishery)) {
    if (sel_type_f[f] == 1L) {
      result$peak_or_inflection[f] <- mu_len + par_sel[f, 1] * sd_len
      result$top_logit_or_width[f] <- exp(par_sel[f, 2]) * sd_len
    } else {
      result$peak_or_inflection[f] <- mu_len + par_sel[f, 1] * sd_len
      result$top_logit_or_width[f] <- par_sel[f, 2]
      result$ascend_se[f] <- par_sel[f, 3] + log(sd_len)
      result$descend_se[f] <- par_sel[f, 4] + log(sd_len)
      result$start_logit[f] <- par_sel[f, 5]
      result$end_logit[f] <- par_sel[f, 6]
    }
  }
  return(result)
}
