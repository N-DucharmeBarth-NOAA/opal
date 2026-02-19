#' Obtain samples from the posterior distribution of reported quantities
#' 
#' When samples from the posterior distribution are obtained using \code{tmbstan}, 
#' only the model parameters can be accessed directly from the \code{stanfit} 
#' object (i.e., derived quantities such as `par_B0` or `M_a` cannot be accessed 
#' directly). Accessing derived quantities must be done using the this function.
#' 
#' @param object A \code{list} specifying the AD object created using \code{MakeADFun}.
#' @param posterior An \code{rstan} objected created using the \code{tmbstan} function.
#' @param pars The parameter(s) to be extracted.
#' @param iters The number of iterations to be extracted.
#' @param option The parallel option to use.
#' @param type default data frame "df" or "array"
#' @return A \code{data.frame} or \code{array}
#' @import foreach
#' @importFrom doParallel registerDoParallel
#' @importFrom SparseNUTS extract_samples
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom utils txtProgressBar
#' @importFrom reshape2 melt
#' @export
#' 
get_posterior <- function(object, posterior, pars = "par_B0", iters = NULL, option = 2, type = "df") {
  post <- extract_samples(fit = posterior, as.list = TRUE)
  chains <- length(post)
  if (is.null(iters)) iters <- dim(post[[1]])[1]
  n_pars <- length(pars)
  n_idx <- integer(n_pars)
  
  r1 <- object$report(par = as.numeric(post[[1]][1,]))
  for (k in 1:n_pars) {
    n_idx[k] <- length(r1[[pars[k]]])
  }
  
  if (option == 1) { # Option 1 (non parallel)
    clist <- vector("list", chains)
    for (j in 1:chains) {
      clist[[j]] <- vector("list", n_pars)
      names(clist[[j]]) <- pars
      for (k in 1:n_pars) {
        clist[[j]][[k]] <- matrix(NA, nrow =  iters, ncol = n_idx[k], dimnames = list(iteration = 1:iters, index = 1:n_idx[k]))
      }
    }
    pb <- txtProgressBar(min = 1, max = chains * iters, style = 3)
    ij <- 0
    for (j in 1:chains) {
      for (i in 1:iters) {
        r1 <- object$report(par = as.numeric(post[[j]][i,])) # Generate a report for each iteration
        for (k in 1:n_pars) {
          clist[[j]][[k]][i,] <- r1[[pars[k]]]
        }
        ij <- ij + 1
        setTxtProgressBar(pb = pb, value = ij)
      }
    }
  }
  
  if (option == 2) { # Option 2 (parallel chains)
    chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
    if (nzchar(chk) && chk == "TRUE") {
      n_cores <- 2L # use 2 cores in CRAN/Travis/AppVeyor
    } else {
      n_cores <- min(detectCores() - 1, chains)
    }
    my_cluster <- makeCluster(n_cores, type = "PSOCK") # create the cluster
    registerDoParallel(cl = my_cluster) # register it to be used by %dopar%
    clist <- foreach(j = 1:chains, .packages = c("decamod")) %dopar% {
      plist <- vector("list", n_pars)
      names(plist) <- pars
      for (k in 1:n_pars) {
        plist[[k]] <- matrix(NA, nrow =  iters, ncol = n_idx[k], dimnames = list(iteration = 1:iters, index = 1:n_idx[k]))
      }
      for (i in 1:iters) {
        r1 <- object$report(par = as.numeric(post[[j]][i,])) # Generate a report for each iteration
        for (k in 1:n_pars) {
          plist[[k]][i,] <- r1[[pars[k]]]
        }
      }
      plist
    }
    stopCluster(cl = my_cluster)
  }
  
  if (type == "df") {
    df <- melt(clist) %>%
      select(parameter = L2, index, chain = L1, iteration, value) %>%
      tibble()
    return(df)
  } else if (type == "array") {
    plist <- vector("list", n_pars)
    for (k in 1:n_pars) {
      plist[[k]] <- vector("list", chains)
      for (j in 1:chains) {
        plist[[k]][[j]] <- vector("list", iters)
        for (i in 1:iters) {
          plist[[k]][[j]][[i]] <- clist[[j]][[k]][[i]]
        }
      }
    }
    larray <- alist <- vector("list", n_pars)
    names(larray) <- pars
    for (k in 1:n_pars) {
      alist[[k]] <- vector("list", chains)
      d <- dim(clist[[1]][[k]][[1]])
      if (is.null(d)) d <- length(clist[[1]][[k]][[1]])
      ld <- length(d)
      alist[[k]] <- array(unlist(plist[[k]]), dim = c(d, iters, chains))
      l1 <- length(dim(alist[[k]]))
      larray[[k]] <- aperm(alist[[k]], c(l1, l1 - 1, 1:ld))
      dimnames(larray[[k]]) <- list(chain = 1:chains, iteration = 1:iters)
    }
    return(larray)
  } else {
    return(clist)
  }
  # head(alist[[1]][,1,1])
  # head(larray[[1]][1,1,])
  # head(clist[[1]][[1]][[1]])
  # head(df$value)
  # tail(alist[[2]][,500,4])
  # tail(larray[[2]][4,500,])
  # tail(clist[[4]][[2]][[500]])
  # tail(df$value)
}

#' Obtain samples from the posterior distribution of reported quantities
#' 
#' When samples from the posterior distribution are obtained using \code{tmbstan}, 
#' only the model parameters can be accessed directly from the \code{stanfit} 
#' object (i.e., derived quantities such as `par_B0` or `M_a` cannot be accessed 
#' directly). Accessing derived quantities must be done using the this function.
#' 
#' @param object A \code{list} specifying the AD object created using \code{MakeADFun}.
#' @param posterior An \code{rstan} objected created using the \code{tmbstan} function.
#' @param pars The parameter(s) to be extracted.
#' @param iters The number of iterations to be extracted.
#' @param option The parallel option to use.
#' @return A \code{data.frame}.
#' @importFrom rstan extract
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @export
#' 
get_posterior3 <- function(object, posterior, pars = "par_B0", iters = NULL, option = 1) {
  
  # Extract the posterior samples for model parameters (note that this does
  # not include any derived quantities).
  post <- extract(object = posterior, pars = "lp__", permuted = FALSE, include = FALSE)
  # post <- as.matrix(posterior)
  # post <- post[, -ncol(post)] # remove the lp__ column
  
  chains <- dim(post)[2]
  if (is.null(iters)) iters <- nrow(post)
  r1 <- object$report(par = post[1, 1, ])
  n_vars <- length(r1[[pars]])
  
  # Option 1 - non parallel
  if (option == 1) {
    mout <- matrix(NA, nrow = chains * iters, ncol = n_vars + 2)
    for (j in 1:chains) {
      for (i in 1:iters) {
        r1 <- object$report(par = post[i, j, ]) # Generate a report for each iteration
        # df <- data.frame(chain = j, iteration = i, output = pars, value = r1[[pars]]) %>%
        #   mutate(id = 1:n())
        mout[(j - 1) * iters + i,] <- c(j, i, r1[[pars]])
        # output <- bind_rows(output, df)
      }
    }
    output <- as.data.frame(mout) %>% 
      pivot_longer(cols = !V1:V2, names_to = "id") %>%
      mutate(id = parse_number(id) - 2)
    names(output) <- c("chain", "iter", "id", "value")
  }
  
  if (option > 1) {
    chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
    
    if (nzchar(chk) && chk == "TRUE") {
      n_cores <- 2L # use 2 cores in CRAN/Travis/AppVeyor
    } else {
      n_cores <- detectCores() - 1
    }
    
    sbt_cluster <- makeCluster(n_cores, type = "PSOCK") # create the cluster
    registerDoParallel(cl = sbt_cluster) # register it to be used by %dopar%
    
    # Option 2 - parallel chains
    if (option == 2) {
      output <- foreach(j = 1:chains, .packages = c("sbt", "dplyr"), .combine = rbind) %dopar% {
        df <- NULL
        for (i in 1:iters) {
          r1 <- object$report(par = post[i, j, ]) # Generate a report for each iteration
          df_chain <- data.frame(chain = j, iteration = i, output = pars, value = r1[[pars]]) %>%
            mutate(id = 1:n())
          df <- bind_rows(df, df_chain)
        }
        df
      }
    }
    
    # Option 3 - parallel chains and iterations
    if (option == 3) {
      output <- foreach(j = 1:chains, .packages = c("sbt", "dplyr"), .combine = rbind) %:%
        foreach(i = 1:iters, .combine = rbind) %dopar% {
          r1 <- object$report(par = post[i, j, ]) # Generate a report for each iteration
          df <- data.frame(chain = j, iteration = i, output = pars, value = r1[[pars]]) %>%
            mutate(id = 1:n())
          df
        }
    }    
    
    stopCluster(cl = sbt_cluster)
  }
  
  df <- output %>% 
    mutate(output = pars) %>%
    relocate(value, .after = last_col())
  
  return(df)
}


#' Obtain samples from the posterior distribution of reported quantities
#' 
#' When samples from the posterior distribution are obtained using \code{tmbstan}, 
#' only the model parameters can be accessed directly from the \code{stanfit} 
#' object (i.e., derived quantities such as `par_B0` or `M_a` cannot be accessed 
#' directly). Accessing derived quantities must be done using the this function.
#' 
#' @param object A \code{list} specifying the AD object created using \code{MakeADFun}.
#' @param posterior An \code{rstan} objected created using the \code{tmbstan} function.
#' @return A \code{data.frame}.
#' @importFrom rstan extract
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @export
#' 
get_posterior2 <- function(object, posterior) {

  reps <- c("lp_lf", "lp_af", "lp_cpue", "lp_tags", "lp_aerial",
            "cpue_pred", "cpue_resid", "aerial_pred", "aerial_resid", "troll_pred", "troll_resid",
            "tag_pred", "tag_resid", "lf_pred", "af_pred",
            "M_a", "spawning_biomass_y")
  
  # Extract the posterior samples for model parameters (note that this does
  # not include any derived quantities).
  post <- extract(posterior, pars = "lp__", permuted = FALSE, include = FALSE)
  
  chains <- dim(post)[2]
  iters <- nrow(post)
  n_vars <- length(reps)
  
  output <- vector("list", n_vars)
  
  for (j in 1:chains) {
    for (i in 1:iters) {
      r1 <- object$report(par = post[i, j, ]) # Generate a report for each iteration
      for (k in 1:n_vars) {
        output[[k]] <- bind_rows(output[[k]], r1[[reps[k]]])
      }
    }
  }
  
  # output <- as.data.frame(mout) %>% 
  #   pivot_longer(cols = !V1:V2, names_to = "id") %>%
  #   mutate(id = parse_number(id) - 2)
  # names(output) <- c("chain", "iter", "id", "value")
  
  df <- output %>% 
    mutate(output = pars) %>%
    relocate(value, .after = last_col())
  
  return(df)
}
