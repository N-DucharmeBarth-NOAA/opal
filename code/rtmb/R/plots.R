#' Plot catch
#' 
#' Plot catch (in thousands of tonnes) by year, season, and fishery.
#' 
#' @param data a \code{list} containing the data that was passed to \code{MakeADFun}.
#' @param obj a \code{list} specifying the AD object created using \code{MakeADFun}.
#' @param posterior an \code{rstan} objected created using the \code{tmbstan} function.
#' @param proj an \code{rstan} objected created using the \code{tmbstan} function.
#' @param probs a numeric vector of probabilities with values in \code{[0,1]} for plotting quantiles of the posterior distribution.
#' @param plot_resid plot the residual or not.
#' @return a \code{ggplot2} object.
#' @import ggplot2
#' @import dplyr
#' @importFrom stats median quantile
#' @importFrom reshape2 melt
#' @importFrom scales breaks_pretty
#' @export
#' 
plot_catch <- function(data, obj, posterior = NULL, proj = NULL, probs = c(0.05, 0.95), plot_resid = FALSE) {
  
  yrs1 <- data$years
  # yrs2 <- data$first_yr_catch:data$last_yr
  # fsh <- c("LL1", "LL2", "LL3", "LL4", "Indonesia", "Australia")
  fsh <- paste0("Fishery: ", 1:data$n_fishery)
  
  df_obs <- reshape2::melt(data$catch_obs_ysf, value.name = "obs") %>%
    filter(.data$obs > 0) %>%
    mutate(season = paste("Season:", 1), fishery = fsh[.data$fishery], Type = "Observed")

  # if (!is.null(proj)) {
  #   df_proj <- melt(proj, value.name = "obs") %>%
  #     filter(.data$obs > 0) %>%
  #     mutate(Season = paste("Season:", .data$Season), Type = "Projected")
  #   df_obs <- bind_rows(df_obs, df_proj)
  # }
  
  df_pred <- obj$report()$catch_pred_ysf %>%
    reshape2::melt(value.name = "pred") %>%
    mutate(year = yrs1[.data$Var1], season = paste("Season:", .data$Var2), fishery = fsh[.data$Var3]) %>%
    right_join(df_obs, by = join_by("year", "season", "fishery")) %>%
    mutate(Fishery = factor(.data$fishery, levels = fsh)) %>%
    mutate(resid = .data$obs - .data$pred)
  
  print(paste0("The maximum catch difference was: ", max(df_pred$resid)))
  
  if (plot_resid) {
    p <- ggplot(data = df_pred, aes(x = .data$year, y = .data$resid)) +
      geom_point(color = "red") +
      labs(x = "Year", y = "Catch residual (tonnes)")
  } else {
    p <- ggplot(data = df_pred, aes(x = .data$year, y = .data$obs)) +
      geom_point(aes(color = .data$Type)) +
      # geom_point(aes(y = .data$pred), shape = 16) +
      geom_line(aes(y = .data$pred), group = 1) +
      labs(x = "Year", y = "Catch (tonnes)", color = NULL) +
      scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05)))
  }
  
  p <- p + 
    facet_wrap(fishery ~ season, scales = "free_y")
    # scale_x_continuous(breaks = breaks_pretty())
  
  return(p)
}

#' Plot age frequency
#' 
#' Plot fit to age frequency (AF) observations.
#' 
#' @param data a \code{list} containing the data that was passed to \code{MakeADFun}.
#' @param object a \code{list} specifying the AD object created using \code{MakeADFun}.
#' @param posterior an \code{rstan} objected created using the \code{tmbstan} function.
#' @param probs a numeric vector of probabilities with values in \code{[0,1]} for plotting quantiles of the posterior distribution.
#' @param years the years to plot.
#' @param fishery the fishery to plot (Indonesian or Australian).
#' @param iters The number of iterations to be extracted.
#' @param ... additional parameters passed on to \code{facet_wrap}.
#' @return a \code{ggplot2} object.
#' @import ggplot2
#' @importFrom scales pretty_breaks
#' @importFrom stats rmultinom
#' @export
#' 
plot_af <- function(data, object, posterior = NULL, probs = c(0.025, 0.975), years = NULL, fishery = "Indonesian", iters = NULL, ...) {
  
  specs <- data.frame(Year = data$af_year + data$first_yr, 
                      Fishery = c("Indonesian", "Australian")[data$af_fishery - 4], 
                      N = data$af_n, min = data$af_min_age, max = data$af_max_age) %>%
    mutate(id = 1:n())
  
  obs <- cbind(specs, data$af_obs) %>%
    data.frame() %>%
    pivot_longer(cols = starts_with("X"), names_to = "Age", values_to = "obs") %>%
    mutate(Age = parse_number(.data$Age))
  
  pred <- cbind(specs, object$report()$af_pred) %>%
    data.frame() %>%
    pivot_longer(cols = starts_with("X"), names_to = "Age", values_to = "pred") %>%
    mutate(Age = parse_number(.data$Age) - 1)
  
  df <- full_join(obs, pred, by = join_by("Year", "Fishery", "N", "min", "max", "Age", "id")) %>% 
    filter(.data$Fishery == fishery, .data$Age >= min, .data$Age <= max)
  
  if (!is.null(years)) df <- df %>% filter(.data$Year %in% years)
  
  dfN <- df %>% 
    select(.data$Year, .data$Fishery, .data$N) %>%
    distinct() %>%
    mutate(N = paste0("N=", .data$N))
  
  p <- ggplot(data = df, aes(x = .data$Age, y = .data$obs)) +
    geom_label(data = dfN, aes(x = -Inf, y = Inf, label = N), hjust = 0, vjust = 1, label.r = unit(0, "lines")) +
    geom_point(colour = "red") +
    geom_line(aes(y = pred), linetype = "dashed") +
    labs(x = "Age", y = "Proportion") +
    # facet_wrap(Year ~ ., ...) +
    facet_wrap(Year ~ .) +
    scale_x_continuous(breaks = pretty_breaks()) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05)))
  
  if (!is.null(posterior)) {
    df0 <- get_posterior(object = object, posterior = posterior, pars = "af_pred", iters = 10) %>%
      mutate(Age = rep(data$age_a, each = 85)[id]) %>%
      mutate(id = rep(1:85, 31)[id]) %>%
      left_join(specs, by = join_by("id")) %>%
      select(-id, -output) %>%
      rename(pred = .data$value) %>%
      filter(.data$Fishery == fishery)
    
    df1 <- df0 %>% pivot_wider(names_from = Age, values_from = pred)
    prob <- as.matrix(df1 %>% select(`0`:`30`))
    
    if (data$af_switch == 1) {
      # Dirichlet
      df_ppred <- t(mapply(rmultinom, n = 1, size = df1$N, prob = split(x = prob, f = c(row(prob)))))
    } else if (data$af_switch == 2) {
      # DM
      df_ppred <- t(mapply(rmultinom, n = 1, size = df1$N, prob = split(x = prob, f = c(row(prob)))))
    } else {
      df_ppred <- t(mapply(rmultinom, n = 1, size = df1$N, prob = split(x = prob, f = c(row(prob)))))
    }
    df_ppred <- df_ppred / rowSums(df_ppred)
    
    dfpp <- cbind(df1 %>% select(chain, iter, Year, Fishery, N, min, max), df_ppred) %>%
      pivot_longer(cols = !chain:max, names_to = "Age", values_to = "ppred") %>%
      mutate(Age = as.numeric(Age) - 1)
    
    df_mcmc <- full_join(df0, dfpp, by = join_by("chain", "iter", "Age", "Year", "Fishery", "N", "min", "max")) %>%
      filter(.data$Age >= min, .data$Age <= max)
    
    p <- p + 
      stat_summary(data = df_mcmc, aes(y = .data$ppred), geom = "ribbon", alpha = 0.5,
                   fun.min = function(x) quantile(x, probs = probs[1]),
                   fun.max = function(x) quantile(x, probs = probs[2])) +
      stat_summary(data = df_mcmc, aes(y = .data$pred), geom = "ribbon", alpha = 0.5,
                   fun.min = function(x) quantile(x, probs = probs[1]),
                   fun.max = function(x) quantile(x, probs = probs[2])) +
      stat_summary(data = df_mcmc, aes(y = .data$pred), geom = "line", fun = median)
  }
  
  return(p)
}

#' Plot LF
#' 
#' Plot LF fit.
#' 
#' @param data a \code{list} containing the data that was passed to \code{MakeADFun}.
#' @param object a \code{list} specifying the AD object created using \code{MakeADFun}.
#' @param posterior an \code{rstan} objected created using the \code{tmbstan} function.
#' @param probs a numeric vector of probabilities with values in \code{[0,1]} for plotting quantiles of the posterior distribution.
#' @param years the years to plot.
#' @param fishery the fishery to plot (LL1, LL2, LL3, or LL4).
#' @param ... additional parameters passed on to \code{facet_wrap}.
#' @return a \code{ggplot2} object.
#' @import ggplot2
#' @importFrom rlang .data
#' @importFrom scales pretty_breaks
#' @importFrom stats rmultinom
#' @export
#' 
plot_lf <- function(data, object, posterior = NULL, probs = c(0.025, 0.975), years = NULL, fishery = "LL1", ...) {
  
  specs <- data.frame(Year = data$lf_year + data$first_yr, 
                      Fishery = c("LL1", "LL2", "LL3", "LL4")[data$lf_fishery], 
                      N = data$lf_n, min = data$lf_minbin) %>%
    mutate(id = 1:n())
  
  obs <- cbind(specs, data$lf_obs) %>%
    data.frame() %>%
    pivot_longer(cols = starts_with("X"), names_to = "Length", values_to = "obs") %>%
    mutate(Length = parse_number(Length))
  
  pred <- cbind(specs, object$report()$lf_pred) %>%
    data.frame() %>%
    pivot_longer(cols = starts_with("X"), names_to = "Length", values_to = "pred") %>%
    mutate(Length = parse_number(Length))
  
  df <- full_join(obs, pred, by = join_by("Year", "Fishery", "N", "min", "Length", "id")) %>% 
    mutate(Length = seq(87.5, 184, 4)[Length]) %>%
    filter(Fishery == fishery, Length >= min)
  
  if (!is.null(years)) df <- df %>% filter(Year %in% years)
  
  dfN <- df %>% 
    select(Year, Fishery, N) %>%
    distinct() %>%
    mutate(N = paste0("N=", N))
  
  p <- ggplot(data = df, aes(x = .data$Length, y = .data$obs)) +
    geom_label(data = dfN, aes(x = -Inf, y = Inf, label = .data$N), hjust = 0, vjust = 1, label.r = unit(0, "lines")) +
    geom_point(colour = "red") +
    geom_line(aes(y = .data$pred), linetype = "dashed") +
    labs(x = "Length (cm)", y = "Proportion") +
    facet_wrap(Year ~ ., ...) +
    scale_x_continuous(breaks = pretty_breaks()) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05)))
  
  if (!is.null(posterior)) {
    df0 <- get_posterior(object = object, posterior = posterior, pars = "lf_pred") %>%
      mutate(Length = rep(1:25, each = 85)[id]) %>%
      mutate(id = rep(1:85, 31)[id]) %>%
      left_join(specs, by = join_by("id")) %>%
      select(-id, -output) %>%
      rename(pred = value) %>%
      filter(Fishery == fishery)
    
    df1 <- df0 %>% pivot_wider(names_from = Age, values_from = pred)
    prob <- as.matrix(df1 %>% select(`0`:`30`))
    
    df_ppred <- t(mapply(rmultinom, n = 1, size = df1$N, prob = split(x = prob, f = c(row(prob)))))
    df_ppred <- df_ppred / rowSums(df_ppred)
    
    dfpp <- cbind(df1 %>% select(chain, iter, Year, Fishery, N, min, max), df_ppred) %>%
      pivot_longer(cols = !chain:max, names_to = "Age", values_to = "ppred") %>%
      mutate(Age = as.numeric(Age) - 1)
    
    df_mcmc <- full_join(df0, dfpp, by = join_by("chain", "iter", "Age", "Year", "Fishery", "N", "min", "max")) %>%
      filter(Age >= min, Age <= max)
    
    p <- p +
      stat_summary(data = df_mcmc, geom = "ribbon", alpha = 0.5, aes(y = .data$ppred),
                   fun.min = function(x) quantile(x, probs = probs[1]),
                   fun.max = function(x) quantile(x, probs = probs[2])) +
      stat_summary(data = df_mcmc, geom = "ribbon", alpha = 0.5, aes(y = .data$pred),
                   fun.min = function(x) quantile(x, probs = probs[1]),
                   fun.max = function(x) quantile(x, probs = probs[2])) +
      stat_summary(data = df_mcmc, aes(y = .data$pred), geom = "line", fun = median)
  }
  
  return(p)
}

#' Plot LF
#' 
#' Plot LF fit.
#' 
#' @param data a \code{list} containing the data that was passed to \code{MakeADFun}.
#' @param object a \code{list} specifying the AD object created using \code{MakeADFun}.
#' @param posterior an \code{rstan} objected created using the \code{tmbstan} function.
#' @param probs a numeric vector of probabilities with values in \code{[0,1]} for plotting quantiles of the posterior distribution.
#' @param years the years to plot.
#' @return a \code{ggplot2} object.
#' @import ggplot2
#' @importFrom rlang .data
#' @importFrom scales pretty_breaks
#' @importFrom stats rmultinom
#' @export
#' 
plot_cpue_lf <- function(data, object, posterior = NULL, probs = c(0.025, 0.975), years = NULL) {
  
  specs <- data.frame(Year = data$cpue_year + data$first_yr - 1, N = data$cpue_n) %>%
    mutate(id = 1:n())
  
  obs <- cbind(specs, data$cpue_lf) %>%
    data.frame() %>%
    pivot_longer(cols = starts_with("X"), names_to = "Length", values_to = "obs") %>%
    mutate(Length = parse_number(.data$Length))
  
  pred <- cbind(specs, object$report()$cpue_lf_pred) %>%
    data.frame() %>%
    pivot_longer(cols = starts_with("X"), names_to = "Length", values_to = "pred") %>%
    mutate(Length = parse_number(.data$Length))
  
  df <- full_join(obs, pred, by = join_by("Year", "N", "Length", "id")) %>% 
    mutate(Length = seq(87.5, 184, 4)[.data$Length])
  
  if (!is.null(years)) df <- df %>% filter(.data$Year %in% years)
  
  dfN <- df %>% 
    select(Year, N) %>%
    distinct() %>%
    mutate(N = paste0("N=", .data$N))
  
  p <- ggplot(data = df, aes(x = .data$Length, y = .data$obs)) +
    geom_label(data = dfN, aes(x = -Inf, y = Inf, label = N), hjust = 0, vjust = 1, label.r = unit(0, "lines")) +
    geom_point(colour = "red") +
    geom_line(aes(y = .data$pred), linetype = "dashed") +
    labs(x = "Length (cm)", y = "Proportion") +
    facet_wrap(Year ~ .) +
    scale_x_continuous(breaks = pretty_breaks()) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05)))
  
  # if (!is.null(posterior)) {
  #   df0 <- get_posterior(object = object, posterior = posterior, pars = "lf_pred") %>%
  #     mutate(Length = rep(1:25, each = 85)[id]) %>%
  #     mutate(id = rep(1:85, 31)[id]) %>%
  #     left_join(specs, by = join_by("id")) %>%
  #     select(-id, -output) %>%
  #     rename(pred = value) %>%
  #     filter(Fishery == fishery)
  #   
  #   df1 <- df0 %>% pivot_wider(names_from = Age, values_from = pred)
  #   prob <- as.matrix(df1 %>% select(`0`:`30`))
  #   
  #   df_ppred <- t(mapply(rmultinom, n = 1, size = df1$N, prob = split(x = prob, f = c(row(prob)))))
  #   df_ppred <- df_ppred / rowSums(df_ppred)
  #   
  #   dfpp <- cbind(df1 %>% select(chain, iter, Year, Fishery, N, min, max), df_ppred) %>%
  #     pivot_longer(cols = !chain:max, names_to = "Age", values_to = "ppred") %>%
  #     mutate(Age = as.numeric(Age) - 1)
  #   
  #   df_mcmc <- full_join(df0, dfpp, by = join_by("chain", "iter", "Age", "Year", "Fishery", "N", "min", "max")) %>%
  #     filter(Age >= min, Age <= max)
  #   
  #   p <- p +
  #     stat_summary(data = df_mcmc, geom = "ribbon", alpha = 0.5,
  #                  aes(y = ppred),
  #                  fun.min = function(x) quantile(x, probs = probs[1]),
  #                  fun.max = function(x) quantile(x, probs = probs[2])) +
  #     stat_summary(data = df_mcmc, geom = "ribbon", alpha = 0.5,
  #                  aes(y = pred),
  #                  fun.min = function(x) quantile(x, probs = probs[1]),
  #                  fun.max = function(x) quantile(x, probs = probs[2])) +
  #     stat_summary(data = df_mcmc, aes(y = pred), geom = "line", fun = median)
  # }
  return(p)
}

#' Plot CPUE
#' 
#' Plot CPUE fit.
#' 
#' @param data a \code{list} containing the data that was passed to \code{MakeADFun}.
#' @param object a \code{list} specifying the AD object created using \code{MakeADFun}.
#' @param posterior an \code{rstan} objected created using the \code{tmbstan} function.
#' @param probs a numeric vector of probabilities with values in \code{[0,1]} for plotting quantiles of the posterior distribution.
#' @param nsim number of simulations to plot.
#' @return a \code{ggplot2} object.
#' @import ggplot2
#' @importFrom rlang .data
#' @importFrom scales pretty_breaks
#' @importFrom stats rlnorm
#' @export
#' 
plot_cpue <- function(data, object, posterior = NULL, probs = c(0.025, 0.975), nsim = 10) {
  
  rep <- object$report(object$env$last.par.best)
  
  df_mle <- data$cpue_data %>%
    mutate(pred = rep$cpue_pred, sigma = rep$cpue_sigma)
  
  # data.frame(year = yrs, obs = data$cpue_obs, pred = object$report()$cpue_pred, 
  #                      sigma = sqrt(data$cpue_sd^2 + object$report()$cpue_sigma^2))

  # mat_sim <- matrix(NA, nrow = nrow(df_mle), ncol = nsim)
  # for (i in seq_len(nsim)) mat_sim[,i] <- exp(obj$simulate()$cpue_log_obs)
  # df_sim <- bind_cols(df_mle, as.data.frame(mat_sim)) %>%
  #   pivot_longer(cols = -names(df_mle), names_to = "iter", values_to = "sim")
  
  p <- ggplot(data = df_mle, aes(x = .data$year, y = .data$value)) +
    # geom_line(data = df_sim, aes(y = .data$sim, group = .data$iter, color = "Simulated"), alpha = 0.5) +
    geom_point(aes(color = "Observed")) +
    geom_linerange(aes(ymin = exp(log(.data$value) - .data$sigma), ymax = exp(log(.data$value) + .data$sigma), color = "Observed")) +
    geom_line(aes(y = .data$pred, color = "Predicted")) +
    labs(x = "Year", y = "CPUE", color = NULL) +
    scale_x_continuous(breaks = pretty_breaks()) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
    facet_wrap(season ~ fishery)
  
  # if (!is.null(posterior)) {
  #   sigma <- get_posterior(object = object, posterior = posterior, pars = "cpue_sigma")
  #   df_mcmc <- get_posterior(object = object, posterior = posterior, pars = "cpue_pred") %>%
  #     mutate(year = yrs[id]) %>%
  #     rename(pred = value) %>% 
  #     left_join(sigma, by = join_by("chain", "iter")) %>%
  #     rowwise() %>%
  #     mutate(ppred = rlnorm(n = 1, meanlog = log(.data$pred), sdlog = value))
  #   # df_ppred <- mapply(rlnorm, n = 1, meanlog = log(df_pred$pred), sdlog = df_pred$value)
  #   # df_mcmc <- df_pred %>% mutate(ppred = df_ppred)
  #   
  #   p <- p + 
  #     stat_summary(data = df_mcmc, geom = "ribbon", alpha = 0.5, aes(y = .data$ppred),
  #                  fun.min = function(x) quantile(x, probs = probs[1]),
  #                  fun.max = function(x) quantile(x, probs = probs[2])) +
  #     stat_summary(data = df_mcmc, geom = "ribbon", alpha = 0.5, aes(y = .data$pred),
  #                  fun.min = function(x) quantile(x, probs = probs[1]),
  #                  fun.max = function(x) quantile(x, probs = probs[2])) +
  #     stat_summary(data = df_mcmc, aes(y = .data$pred), geom = "line", fun = median)
  # }
  
  return(p)
}

#' Plot spawning biomass
#' 
#' Plot the spawning biomass (tonnes) or relative spawning biomass by year for a 
#' single model run, a grid of model runs, or an MCMC.
#' 
#' @param data_list a \code{list} containing the data that was passed to \code{MakeADFun}.
#' @param object_list a \code{list} specifying the AD object created using the \code{MakeADFun} function.
#' @param posterior an \code{rstan} objected created using the \code{tmbstan} function.
#' @param probs a numeric vector of probabilities with values in \code{[0,1]} for plotting quantiles of the posterior distribution. Defaults to the 90% credible interval.
#' @param relative if the plot should be relative spawning biomass.
#' @param labels a \code{vector} of labels for the model runs.
#' @return a \code{ggplot2} object.
#' @import ggplot2
#' @import dplyr
#' @importFrom methods is
#' @importFrom stats median quantile
#' @importFrom scales pretty_breaks
#' @export
#' 
plot_biomass_spawning <- function(data_list, object_list, posterior = NULL, probs = c(0.025, 0.975), relative = TRUE, labels = NULL) {
  
  n_model <- length(data_list)
  if (is.null(labels)) labels <- 1:n_model
  
  dfj <- vector("list", n_model)
  for (j in 1:n_model) {
    yrs <- data_list[[j]]$first_yr:(data_list[[j]]$last_yr + 1)
    rep <- object_list[[j]]$report()
    dfj[[j]] <- data.frame(Model = labels[j], year = yrs, B0 = rep$B0, value = rep$spawning_biomass_y)
  }
  df1 <- bind_rows(dfj) %>% mutate(Model = factor(Model, levels = labels))
  
  # dfi <- list()
  # for (i in 1:n_grid) {
  #   rep <- obj[[i]]$report()
  #   dfi[[i]] <- data.frame(grid = i, year = yrs, B0 = rep$B0, value = rep$spawning_biomass_y)
  # }
  # df1 <- bind_rows(dfi)
  
  if (relative) {
    df1 <- df1 %>% mutate(value = .data$value / .data$B0)
    ylab <- "Relative spawning biomass"
  } else {
    df1 <- df1 %>% mutate(value = .data$value / 1e6)
    ylab <- "Spawning biomass (millions of tonnes)"
  }
  
  p <- ggplot(data = df1, aes(x = .data$year, y = .data$value, color = .data$Model)) +
    # geom_hline(yintercept = B01 / 1e6, linetype = "dashed") +
    # geom_line(aes(colour = factor(h), group = grid), linetype = "dashed") +
    # geom_line(aes(group = .data$grid), linetype = "dashed") +
    geom_line() +
    scale_x_continuous(breaks = pretty_breaks()) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
    labs(x = "Year", y = ylab)
  if (n_model == 1) p <- p + theme(legend.position = "none")
  
  # MCMC
  # if (!is.null(posterior)) {
  #   if (is(posterior, "adfit")) {
  #     B02 <- get_posterior(object = object, posterior = posterior, pars = "par_B0")
  #     df_mcmc <- get_posterior(object = object, posterior = posterior, pars = "spawning_biomass_y") %>%
  #       mutate(year = yrs[id]) %>% 
  #       left_join(B02, by = join_by("chain", "iter")) %>%
  #       mutate(mc_grid = 1)
  #   # } else if (is(posterior, "stanfit")) {
  #   #   B02 <- get_posterior(object = object, posterior = posterior, pars = "par_B0")
  #   #   df_mcmc <- get_posterior(object = object, posterior = posterior, pars = "spawning_biomass_y") %>%
  #   #     mutate(year = yrs[id]) %>% 
  #   #     left_join(B02, by = join_by("chain", "iter")) %>%
  #   #     mutate(mc_grid = 1)
  #   } else if (is.list(posterior)) {
  #     B0_mcmc <- NULL
  #     bio_mcmc <- NULL
  #     for (j in 1:length(posterior)) {
  #       B0_j <- get_posterior(object = object, posterior = posterior[[j]], pars = "par_B0") %>%
  #         mutate(mc_grid = j)
  #       bio_j <- get_posterior(object = object, posterior = posterior[[j]], pars = "spawning_biomass_y") %>%
  #         mutate(mc_grid = j)
  #       B0_mcmc <- bind_rows(B0_mcmc, B0_j)
  #       bio_mcmc <- bind_rows(bio_mcmc, bio_j)
  #     }
  #     df_mcmc <- bio_mcmc %>%
  #       mutate(year = yrs[id]) %>% 
  #       left_join(B0_mcmc, by = join_by("chain", "iter", "mc_grid"))
  #   }
  #   
  #   if (relative) {
  #     df_mcmc <- df_mcmc %>% mutate(value = .data$value.x / .data$value.y)
  #   } else {
  #     df_mcmc <- df_mcmc %>% mutate(value = .data$value.x / 1e6)
  #   }
  #   
  #   p <- p + 
  #     stat_summary(data = df_mcmc, geom = "ribbon", alpha = 0.5, aes(fill = factor(.data$mc_grid)),
  #                  fun.min = function(x) quantile(x, probs = probs[1]),
  #                  fun.max = function(x) quantile(x, probs = probs[2])) +
  #     stat_summary(data = df_mcmc, geom = "line", fun = median, aes(colour = factor(.data$mc_grid)))
  # }
  
  return(p)
}

#' Plot initial numbers
#' 
#' Plot initial numbers by age.
#' 
#' @param data a \code{list} containing the data that was passed to \code{MakeADFun}.
#' @param object a \code{list} specifying the AD object created using \code{MakeADFun}.
#' @param posterior an \code{rstan} objected created using the \code{tmbstan} function.
#' @param probs a numeric vector of probabilities with values in \code{[0,1]} for plotting quantiles of the posterior distribution.
#' @return a \code{ggplot2} object.
#' @import ggplot2
#' @importFrom stats median quantile
#' @export
#' 
plot_initial_numbers <- function(data, object, posterior = NULL, probs = c(0.025, 0.975)) {
  ages <- data$min_age:data$max_age
  num <- object$report()$number_ysa
  N_age1 <- data.frame(age = ages, value = num[1,1,])
  p <- ggplot(data = N_age1, aes(x = .data$age, y = .data$value / 1e6)) +
    geom_line(data = N_age1, linetype = "dashed") +
    labs(x = "Age", y = "Initial numbers (millions)") +
    scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05)))
  # MCMC
  # if (!is.null(posterior)) {
  #   M_age2 <- get_posterior(object = object, posterior = posterior, parameter = "M_age")
  #   
  #   p <- p + 
  #     stat_summary(data = M_age2, geom = "ribbon", alpha = 0.5, 
  #                  fun.min = function(x) quantile(x, probs = probs[1]), 
  #                  fun.max = function(x) quantile(x, probs = probs[2])) + 
  #     stat_summary(data = M_age2, geom = "line", fun = median)
  # }
  return(p)
}

#' Plot harvest rate
#' 
#' Plot harvest rate by year, season, and age.
#' 
#' @param data a \code{list} containing the data that was passed to \code{MakeADFun}.
#' @param object a \code{list} specifying the AD object created using \code{MakeADFun}.
#' @param posterior an \code{rstan} objected created using the \code{tmbstan} function.
#' @param probs a numeric vector of probabilities with values in \code{[0,1]} for plotting quantiles of the posterior distribution.
#' @param years the years to show on the plot.
#' @param ... options passed on to \code{geom_density_ridges}.
#' @return a \code{ggplot2} object.
#' @import ggplot2
#' @import dplyr
#' @importFrom stats median quantile
#' @importFrom reshape2 melt
#' @importFrom ggridges geom_density_ridges
#' @importFrom scales pretty_breaks
#' @export
#' 
plot_hrate <- function(data, object, posterior = NULL, probs = c(0.025, 0.975), years = 2013:2022, ...) {
  
  yrs <- data$first_yr:data$last_yr
  ages <- data$min_age:data$max_age

  df <- object$report()$hrate_ysa %>%
    melt() %>%
    mutate(year = yrs[Var1], season = Var2, age = ages[Var3]) %>%
    filter(year %in% years, year >= data$first_yr_catch)

  p <- ggplot(data = df, aes(x = age, y = year, height = value, group = year)) +
    geom_density_ridges(stat = "identity", alpha = 0.75, rel_min_height = 0, color = NA, ...) + 
    facet_wrap(season ~ .) +
    labs(x = "Age", y = "Year") +
    scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
    scale_y_reverse(breaks = pretty_breaks())
  
  return(p)
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
