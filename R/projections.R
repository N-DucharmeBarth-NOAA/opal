#' Project dynamics
#' 
#' @param data a \code{list} of parameter values.
#' @param object a \code{list} of parameter values.
#' @param mcmc a \code{list} of parameter values.
#' @param n_proj the number of projection years.
#' @param n_iter a \code{list} of inputs.
#' @param rdev_y a \code{list} of inputs.
#' @param sel_fya a \code{list} of inputs.
#' @param catch_ysf a \code{list} of inputs.
#' @return the negative log-likelihood (NLL) value.
#' @importFrom SparseNUTS extract_samples
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export
#' 
project_dynamics <- function(data, object, mcmc = NULL, n_proj = 5, n_iter = 1, rdev_y, sel_fya, catch_ysf) {
  
  if (!is.null(mcmc)) {
    post <- extract_samples(fit = mcmc)
    n_post <- nrow(post)
    if (n_iter > n_post) stop("Cannot simulate this many iterations from the 'mcmc' input.")
  }
  
  dr <- dim(rdev_y)
  if (length(dr) != 2) stop("Incorrect dimensions for the array 'rdev_y'.")
  if (n_iter > dr[1]) stop("Cannot simulate this many iterations from the array 'rdev_y'.")
  
  ds <- dim(sel_fya)
  if (length(ds) != 4) stop("Incorrect dimensions for the array 'sel_fya'.")
  if (n_iter > ds[1]) stop("Cannot simulate this many iterations from the array 'sel_fya'.")
  
  proj_years <- (data$last_yr + 1):(data$last_yr + n_proj)
  
  # dn <- dimnames(data$catch_obs_ysf)
  # proj_catch_ysf <- array(0, dim = list(n_proj, data$n_season, data$n_fishery), dimnames = list(Year = proj_years, Season = dn$Season, Fishery = dn$Fishery))
  # for (y in 1:n_proj) proj_catch_ysf[y,,] <- data$catch_obs_ysf[71,,]
  
  dn <- dimnames(data$weight_fya)
  proj_weight_fya <- array(0, dim = list(data$n_fishery, n_proj, data$n_age), dimnames = list(fishery = dn$fishery, year = proj_years, age = dn$age))
  for (y in 1:n_proj) proj_weight_fya[,y,] <- data$weight_fya[,data$n_year,]
  
  dn <- dimnames(data$af_sliced_ysfa)
  proj_af_sliced_ysfa <- array(0, dim = list(n_proj, data$n_season, data$n_fishery, data$n_age), dimnames = list(year = proj_years, season = dn$Season, fishery = dn$Fishery, age = dn$Age))
  for (y in 1:n_proj) proj_af_sliced_ysfa[y,,,] <- data$af_sliced_ysfa[data$n_year,,,]
  
  for (f in seq_len(data$n_fishery)) {
    if (data$removal_switch_f[f] > 0) {
      # x1 <- data$catch_obs_ysf[,,f]
      # x2 <- apply(data$af_sliced_ysfa[,,f,], 1:2, sum)
      x1 <- proj_catch_ysf[,,f]
      x2 <- apply(proj_af_sliced_ysfa[,,f,], 1:2, sum)
      df1 <- as.data.frame(x1) %>%
        rownames_to_column("Year") %>%
        pivot_longer(cols = c(`1`, `2`), names_to = "Season", values_to = "catch")
      df2 <- as.data.frame(x2) %>%
        rownames_to_column("Year") %>%
        pivot_longer(cols = c(`1`, `2`), names_to = "Season", values_to = "AF")
      df <- left_join(df1, df2, by = c("Year", "Season")) %>%
        mutate(Year = as.numeric(Year)) %>%
        mutate(violation = (catch > 0 & AF == 0) | (AF > 0 & catch == 0)) %>%
        filter(violation)
      if (nrow(df) > 0) stop("You have specified catch for a fishery using direct removals but no associated sliced AFs.")
      # data$catch_obs_ysf[39:44,1,f]
      # rowSums(data$af_sliced_ysfa[60:65,1,f,])
    }
  }
  
  dyn <- vector("list", n_iter)
  if (n_iter > 1) pb <- txtProgressBar(min = 1, max = n_iter, style = 3)
  rep1 <- object$report(object$env$last.par.best)
  
  for (i in seq_len(n_iter)) {
    if (!is.null(mcmc)) {
      rep <- object$report(as.numeric(post[i,]))
    } else {
      rep <- rep1
    }
    
    # phi_ya <- get_phi(log_psi = rep$par_log_psi, data$length_m50, data$length_m95, data$length_mu_ysa, data$length_sd_a, data$dl_yal)
    phi_ya <- rep$phi_ya
    dn <- list(year = data$first_yr:(data$last_yr + 1), age = data$age_a)
    dimnames(phi_ya) <- dn
    proj_phi_ya <- array(0, dim = list(n_proj + 1, data$n_age), dimnames = list(year = c(proj_years, max(proj_years) + 1), age = dn$age))
    for (y in 1:(n_proj + 1)) proj_phi_ya[y,] <- phi_ya[data$n_year + 1,]
    
    dyn[[i]] <- do_dynamics(
      first_yr = 1, first_yr_catch = 1, 
      B0 = rep$B0, R0 = rep$R0, alpha = rep$alpha, beta = rep$beta, h = rep$par_h, sigma_r = rep$sigma_r,
      rdev_y = rdev_y[i,], # projected recruitment deviates
      M_a = rep$M_a, 
      phi_ya = proj_phi_ya, # projected phi
      init_number_a = rep$number_ysa[data$n_year + 1, 1,], # different N
      removal_switch_f = data$removal_switch_f, 
      catch_obs_ysf = catch_ysf, # projected catch
      sel_fya = sel_fya[i,,,], # projected selectivity
      weight_fya = proj_weight_fya, # projected weight
      af_sliced_ysfa = data$af_sliced_ysfa)
    
    if (n_iter > 1) setTxtProgressBar(pb, i)
  }

  return(dyn)
}

#' Project recruitment deviates
#' 
#' @param data a \code{list} of parameter values.
#' @param obj a \code{list} of parameter values.
#' @param mcmc a \code{list} of parameter values.
#' @param first_yr the first year.
#' @param last_yr the last year.
#' @param n_proj the number of projection years.
#' @param n_iter a \code{list} of inputs.
#' @param arima default = TRUE, FALSE = "lognormal"
#' @return the negative log-likelihood (NLL) value.
#' @importFrom SparseNUTS extract_samples
#' @importFrom forecast auto.arima
#' @importFrom stats simulate sd
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export
#' 
project_selectivity <- function(data, obj, mcmc = NULL, 
                                first_yr = 2000, last_yr = NULL, n_proj = 5, n_iter = NULL, arima = TRUE) {
  
  if (is.null(last_yr)) last_yr <- data$last_yr
  proj_years <- (data$last_yr + 1):(data$last_yr + n_proj)
  samp_years <- (first_yr - data$first_yr + 1):(last_yr - data$first_yr + 1)
  sel_fya <- obj$report()$sel_fya[,samp_years,]
  
  n_fishery <- dim(sel_fya)[1]
  n_age <- dim(sel_fya)[3]
  
  # if (!is.null(mcmc)) {
  #   post <- extract_samples(fit = mcmc)
  #   rdevs1 <- as.matrix(post[grepl("par_rdev_y", names(post))])
  #   if (!is.null(n_iter)) {
  #     ii <- sample(x = 1:nrow(post), size = n_iter)
  #     rdevs1 <- rdevs1[ii,]
  #   } else {
  #     n_iter <- nrow(post)
  #   }
  # } else {
  # if (is.null(n_iter)) n_iter <- 1
  # x <- obj$env$last.par.best[names(obj$par) %in% "par_rdev_y"]
  # rdevs1 <- matrix(x, nrow = n_iter, ncol = length(x), byrow = TRUE)
  # }
  # rdevs <- rdevs1[,samp_years]
  
  sim_ifya <- array(0, dim = list(n_iter, n_fishery, n_proj, n_age), 
                    dimnames = list(iteration = 1:n_iter, fishery = 1:n_fishery, year = proj_years, age = data$age_a))
  
  removal_switch_f <- c(data$removal_switch_f, 0)
  
  for (f in seq_len(n_fishery)) {
    for (a in seq_len(n_age)) {
      for (i in seq_len(n_iter)) {
        y <- log(sel_fya[f,,a])
        if (all(is.finite(y)) & removal_switch_f[f] == 0) {
          if (sd(y) > 0) {
            # fit <- ar(y)
            # sim_ifya[i, f,, a] <- exp(arima.sim(n = n_proj, list(ar = fit$ar), sd = sqrt(fit$var.pred)))
            sim_ifya[i, f,, a] <- exp(rnorm(n = n_proj, mean(y), sd(y)))
          } else {
            sim_ifya[i, f,, a] <- exp(mean(y))
          }
        }
      }
    }
  }
  
  return(sim_ifya)
}

#' Project recruitment deviates
#' 
#' @param data a \code{list} of parameter values.
#' @param obj a \code{list} of parameter values.
#' @param mcmc a \code{list} of parameter values.
#' @param first_yr a \code{list} of inputs.
#' @param last_yr a \code{list} of inputs.
#' @param n_proj a \code{list} of inputs.
#' @param n_iter a \code{list} of inputs.
#' @param max.p Maximum value of p, or the maximum value of p (the AR order) to consider.
#' @param max.d Maximum value of d, or the maximum value of q (the MA order) to consider.
#' @param max.q Maximum value of q
#' @param arima default = TRUE, FALSE = "lognormal"
#' @return a \code{list} of projected recruitment deviates and ARIMA specifications.
#' @importFrom SparseNUTS extract_samples
#' @importFrom forecast auto.arima
#' @importFrom stats simulate ar arima.sim
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export
#' 
project_rec_devs <- function(data, obj, mcmc = NULL, first_yr = 1931, last_yr = NULL, n_proj = 5, n_iter = NULL,
                             max.p = 5, max.d = 5, max.q = 5, arima = TRUE) {
  
  if (!is.null(mcmc)) {
    post <- extract_samples(fit = mcmc)
    rdevs1 <- as.matrix(post[grepl("par_rdev_y", names(post))])
    if (!is.null(n_iter)) {
      ii <- sample(x = 1:nrow(post), size = n_iter)
      rdevs1 <- rdevs1[ii,]
    } else {
      n_iter <- nrow(post)
    }
  } else {
    if (is.null(n_iter)) n_iter <- 1
    x <- obj$env$last.par.best[names(obj$par) %in% "par_rdev_y"]
    rdevs1 <- matrix(x, nrow = n_iter, ncol = length(x), byrow = TRUE)
    # rdevs1 <- array(obj$par[names(obj$par) %in% "par_rdev_y"], dim = c(n_iter, n_year))
    # rdevs1 <- t(as.matrix(obj$par[names(obj$par) %in% "par_rdev_y"]))
    # rdevs1 <- t(as.matrix(obj$report()$par_rdev_y))
  }
  dimnames(rdevs1) <- list(iteration = 1:n_iter, year = data$first_yr:data$last_yr)
  
  if (is.null(last_yr)) last_yr <- data$last_yr
  samp_years <- (first_yr - data$first_yr + 1):(last_yr - data$first_yr + 1)
  rdevs <- rdevs1[,samp_years]
  
  proj_years <- (data$last_yr + 1):(data$last_yr + n_proj)
  sim <- matrix(NA, nrow = n_iter, ncol = n_proj, dimnames = list(iteration = 1:n_iter, year = proj_years))
  arima_pars <- matrix(NA, nrow = n_iter, ncol = 3)
  
  if (arima) {
    if (n_iter > 1) pb <- txtProgressBar(min = 1, max = n_iter, style = 3)
    for (i in seq_len(n_iter)) {
      fit <- auto.arima(y = rdevs[i,], max.p = max.p, max.d = max.d, max.q = max.q, approximation = FALSE, stepwise = FALSE, ic = "bic", trace = FALSE)
      # fc <- forecast(object = fit, h = h)
      sim[i,] <- simulate(object = fit, nsim = n_proj, future = TRUE, bootstrap = TRUE)
      arima_pars[i,] <- fit$arma[1:3]
      if (n_iter > 1) setTxtProgressBar(pb, i)
    }
  } else {
    for (i in seq_len(n_iter)) {
      fit <- ar(x = rdevs[i,], order.max = max.p)
      sim[i,] <- arima.sim(n = n_proj, list(ar = fit$ar), sd = sqrt(fit$var.pred))
      arima_pars[i, 1] <- fit$order
      # devs <- rdevs[i,]
      # rho <- get_rho(first_yr = first_yr, last_yr = last_yr, rdev = devs)
      # prev_dev <- devs[length(devs)]
      # for (y in seq_len(n_proj)) {
      #   sim[i, y] <- rho * prev_dev + sqrt(1 - rho^2) * rnorm(n = 1, mean = mean(devs), sd = sd(devs))
      #   prev_dev <- sim[i, y]
      # }
    }
  }
  
  # rdev_y <- sim
  # dimnames(rdev_y) <- list(iteration = 1:n_iter, year = (data$last_yr + 1):(data$last_yr + n_proj))
  return(list(rdev_y = sim, arima_pars = arima_pars))
}

#' #' Extend year-varying value into projection years
#' #' 
#' #' @param value a \code{array} where the first dimension is years y
#' #' @param par parameter name used to determine dimensions and indexing (e.g. ytrsfl)
#' #' @param n_proj number of years to add to array
#' #' @return a \code{array} including all iterations of the derived value
#' #' @export
#' #' 
#' extend_years <- function(value, par, n_proj) {
#'   d <- dim(value)
#'   ## vector
#'   if (is.null(d)) {
#'     val_new <- value[length(value)]
#'     val_out <- rep(val_new, n_proj)
#'   } else {
#'     ## array
#'     d2 <- strsplit(par, "_y")[[1]][2]
#'     d2 <- paste0("y", d2)
#'     index <- strsplit(d2, "")[[1]]
#'     if (length(index) < length(d)) {
#'       index <- c("i", index)
#'     }
#'     
#'     ## remove final year
#'     last_year <- d[which(index == "y")]
#'     if (index[1] == "i") {
#'       if (nchar(d2) == 6) val_new <- value[,last_year,,,,,, drop = FALSE]
#'       if (nchar(d2) == 5) val_new <- value[,last_year,,,,, drop = FALSE]
#'       if (nchar(d2) == 4) val_new <- value[,last_year,,,, drop = FALSE]
#'       if (nchar(d2) == 3) val_new <- value[,last_year,,, drop = FALSE]
#'       if (nchar(d2) == 2) val_new <- value[,last_year,, drop = FALSE]
#'     } else {
#'       if (nchar(d2) == 6) val_new <- value[last_year,,,,,, drop = FALSE]
#'       if (nchar(d2) == 5) val_new <- value[last_year,,,,, drop = FALSE]
#'       if (nchar(d2) == 4) val_new <- value[last_year,,,, drop = FALSE]
#'       if (nchar(d2) == 3) val_new <- value[last_year,,, drop = FALSE]
#'       if (nchar(d2) == 2) val_new <- value[last_year,, drop = FALSE]
#'     }
#'     
#'     ## projection dimensions and values
#'     dn <- dim(val_new)
#'     dn[which(index == "y")] <- n_proj
#'     
#'     val_out <- array(NA, dim = dn)
#'     for (y in 1:n_proj) {
#'       if (index[1] == "i") {
#'         if (nchar(d2) == 6) val_out[,y,,,,,] <- val_new
#'         if (nchar(d2) == 5) val_out[,y,,,,] <- val_new
#'         if (nchar(d2) == 4) val_out[,y,,,] <- val_new
#'         if (nchar(d2) == 3) val_out[,y,,] <- val_new
#'         if (nchar(d2) == 2) val_out[,y,] <- val_new
#'       } else {
#'         if (nchar(d2) == 6) val_out[y,,,,,] <- val_new
#'         if (nchar(d2) == 5) val_out[y,,,,] <- val_new
#'         if (nchar(d2) == 4) val_out[y,,,] <- val_new
#'         if (nchar(d2) == 3) val_out[y,,] <- val_new
#'         if (nchar(d2) == 2) val_out[y,] <- val_new
#'       }
#'     }
#'   }
#'   
#'  return(val_out)
#' }
#' 
#' #' Prepare data and values lists for projections
#' #' 
#' #' @param data \code{list} data values
#' #' @param parameters \code{list} input parameters
#' #' @param obj from initial fit
#' #' @param adfit NULL if not included
#' #' @param n_proj number of years to add to array
#' #' @param iter sample of posterior distribution (single value only)
#' #' @return a \code{list} with years extended for projection
#' #' @importFrom SparseNUTS extract_samples
#' #' @importFrom parallel detectCores makeCluster stopCluster
#' #' @importFrom doParallel registerDoParallel
#' #' @importFrom abind abind
#' #' @export
#' #' 
#' prepare_proj <- function(data, parameters, obj, adfit, n_proj, iter = NULL) {
#'   data_out <- list()
#'   for (i in 1:length(data)) {
#'     par <- names(data[i])
#'     if (par == "handling_mortality_y") {
#'       data_out[[i]] <- extend_years(data[[i]], par = par, n_proj = n_proj)
#'     } else {
#'       if (grepl("_y", par) && !(grepl("_year", par))) {
#'         data_out[[i]] <- extend_years(data[[i]], par = par, n_proj = n_proj)
#'       } else {
#'         data_out[[i]] <- data[[i]]
#'       }
#'     }
#'   }
#'   names(data_out) <- names(data)
#'   
#'   # MAP
#'   if (all(is.null(adfit))) {
#'     if (any(grepl('logit_U0', names(obj$par)))) {
#'       map_list <- list(
#'         logit_h = obj$report()$logit_h,
#'         F0 = as.numeric(exp(obj$par["log_F0"])),
#'         U0 = as.numeric(exp(obj$par['logit_U0'])),
#'         R0 = as.numeric(exp(obj$par["log_R0"])),
#'         recruitment_size_sl = obj$report()$recruitment_size_sl,
#'         numbers_rsl = obj$report()$numbers_ytrsl[data$n_year, 2,,,],
#'         Rdev_yr = obj$report()$Rdev_yr,
#'         Rsigma = obj$report()$Rsigma,
#'         growth_ytrsll = obj$report()$growth_ytrsll,
#'         M_ytrsl = obj$report()$M_ytrsl,
#'         F_ytrf = obj$report()$F_ytrf,
#'         U_ytrf = obj$report()$U_ytrf,
#'         selectivity_ytrsfl = obj$report()$selectivity_ytrsfl
#'       )      
#'     } else {
#'       map_list <- list(
#'         logit_h = obj$report()$logit_h,
#'         F0 = as.numeric(exp(obj$par["log_F0"])),
#'         U0 = plogis(parameters$logit_U0),
#'         R0 = as.numeric(exp(obj$par["log_R0"])),
#'         recruitment_size_sl = obj$report()$recruitment_size_sl,
#'         numbers_rsl = obj$report()$numbers_ytrsl[data$n_year, 2,,,],
#'         Rdev_yr = obj$report()$Rdev_yr,
#'         Rsigma = obj$report()$Rsigma,
#'         growth_ytrsll = obj$report()$growth_ytrsll,
#'         M_ytrsl = obj$report()$M_ytrsl,
#'         F_ytrf = obj$report()$F_ytrf,
#'         U_ytrf = obj$report()$U_ytrf,
#'         selectivity_ytrsfl = obj$report()$selectivity_ytrsfl
#'       )
#'     }
#' 
#'     map_out <- list()
#'     for (i in 1:length(map_list)) {
#'       par <- names(map_list)[i]
#'       if (grepl("_ytr", par) | grepl("_yr", par)) {
#'         map_out[[i]] <- extend_years(map_list[[i]], par = par, n_proj = n_proj)
#'       } else {
#'         map_out[[i]] <- map_list[[i]]
#'       }
#'     }
#'     names(map_out) <- names(map_list)
#'     mcmc_out <- NULL
#'   } else { 
#'     # MCMC
#'     map_out <- NULL
#'     post <- extract_samples(adfit)
#' 
#'     if (any(grepl('logit_U0', names(obj$par)))) {
#'       mcmc_list <- list(
#'         U0 = as.numeric(plogis(post[,grep("logit_U0", colnames(post))])),
#'         R0 = as.numeric(exp(post[,grep("log_R0", colnames(post))])),
#'         recruitment_size_sl = obj$report()$recruitment_size_sl
#'       )      
#'     } else {
#'       mcmc_list <- list(
#'         U0 = plogis(parameters$logit_U0),
#'         R0 = as.numeric(exp(post[,grep("log_R0", colnames(post))])),
#'         recruitment_size_sl = obj$report()$recruitment_size_sl
#'       ) 
#'     }
#' 
#'     
#'     post2 <- get_posterior(object = obj, 
#'                             posterior = mcmc,
#'                             pars = c("logit_h", "Rsigma", "numbers_ytrsl", "growth_ytrsll", "M_ytrsl", "U_ytrf", "selectivity_ytrsfl"),
#'                             option = 2,
#'                             type = "list")
#'     
#'     n_chains <- length(post2)
#'     n_par <- length(post2[[1]])
#'     
#'     ## to have length par
#'     outlist <- list()
#'     for (j in 1:n_chains) {
#'       parlist <- post2[[j]]
#'       for (k in 1:n_par) {
#'         iterlist <- parlist[[k]]
#'         out1 <- abind(iterlist, along = 0)
#'         if (j == 1) {
#'           outlist[[k]] <- out1
#'         } else {
#'           outlist[[k]] <- abind(outlist[[k]], out1, along = 1)
#'         }
#'       }
#'     }
#'     names(outlist) <- names(post2[[1]])
#'     
#'     num_ytrsl <- outlist[["numbers_ytrsl"]]
#'     outlist[["numbers_ytrsl"]] <- array(num_ytrsl[,data$n_year,2,,,], dim = c(dim(num_ytrsl)[1], data$n_region, data$n_sex, n_length))
#'     names(outlist)[which(names(outlist) == "numbers_ytrsl")] <- "numbers_rsl"
#'     
#'     mcmc_list <- c(mcmc_list, outlist)
#'     mcmc_out <- list()
#'     for (i in 1:length(mcmc_list)) {
#'       par <- names(mcmc_list)[i]
#'       if (grepl("_y", par)) {
#'         mcmc_out[[i]] <- extend_years(value = mcmc_list[[i]], par = par, n_proj = n_proj)
#'       } else {
#'         mcmc_out[[i]] <- mcmc_list[[i]]
#'       }
#'     }
#'     names(mcmc_out) <- names(mcmc_list)
#'   }
#'   proj_values <- c(data_out, map_out, mcmc_out)
#'   return(proj_values)
#' }
