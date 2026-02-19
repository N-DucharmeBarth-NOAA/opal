#' Set up a grid
#' 
#' @param parameters a \code{list} containing the parameter inputs.
#' @param m0 the M0 values to be included in the grid.
#' @param m10 the M10 values to be included in the grid.
#' @param h the h values to be included in the grid.
#' @param psi the psi values to be included in the grid.
#' @return a \code{list} of parameter inputs ready to be passed to \code{MakeADFun}.
#' @export
#' @examples
#' #parameters <- get_parameters(data = data)
#' #grid_parameters <- get_grid(parameters = parameters)
#' 
get_grid <- function(parameters,
                     m0 = c(0.5, 0.4, 0.45), 
                     m10 = c(0.065, 0.085, 0.105), 
                     h = c(0.55, 0.63, 0.72, 0.8), 
                     psi = c(1.75, 1.5, 2)) {
  
  grid <- expand.grid(m0 = m0, m10 = m10, h = h, psi = psi)
  N <- nrow(grid)
  par_list <- rep(list(parameters), times = N)
  for (i in 1:N) {
    par_list[[i]]$par_log_psi <- log(grid$psi[i])
    par_list[[i]]$par_log_m0 <- log(grid$m0[i])
    par_list[[i]]$par_log_m10 <- log(grid$m10[i])
    par_list[[i]]$par_log_h <- log(grid$h[i])
  }
  names(par_list) <- apply(X = grid, MARGIN = 1, FUN = function(x) paste(paste0(names(grid), "=", x), collapse = ", "))
  return(par_list)
}

#' Run a grid
#' 
#' Fit a model for each grid cell defined in \code{grid} where each grid cell has a
#' different combination of fixed parameter values.
#' 
#' @param data a \code{list} containing the data created using the \code{get_data} function.
#' @param grid_parameters a \code{list} of parameter inputs for each grid cell. This can be created using \code{get_grid}.
#' @param bounds the lower and upper parameter bounds created using the \code{get_bounds} function.
#' @param map a \code{list} defining how to optionally collect and fix parameters.
#' @param random a character \code{vector} defining the random effect parameters.
#' @param control a \code{list} defining the controls passed to \code{nlminb}.
#' @return a \code{list} of grid cells.
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom stats nlminb
#' @import foreach
#' @export
#' 
run_grid <- function(data, grid_parameters, bounds, map = list(), random = c(), control = list(eval.max = 10000, iter.max = 10000)) {
  
  n_grid <- length(grid_parameters)
  # n_cores <- detectCores() - 1
  # sbt_cluster <- makeCluster(n_cores, type = "PSOCK") # create the cluster
  # registerDoParallel(cl = sbt_cluster) # register it to be used by %dopar%
  # # print(sbt_cluster) # check cluster definition (optional)
  # # getDoParRegistered() # check if it is registered (optional)
  # # getDoParWorkers() # how many workers are available? (optional)
  grid_obj <- vector("list", length = n_grid)
  
  # grid_obj <- foreach(i = 1:n_grid, .packages = "sbt") %dopar% {
  for (i in 1:n_grid) {
    obj <- MakeADFun(func = cmb(sbt_model, data), parameters = grid_parameters[[i]], map = map)
    opt <- nlminb(start = obj$par, objective = obj$fn, gradient = obj$gr, hessian = obj$he, lower = bounds$lower, upper = bounds$upper, control = control)
    opt <- nlminb(start = opt$par, objective = obj$fn, gradient = obj$gr, hessian = obj$he, lower = bounds$lower, upper = bounds$upper, control = control)
    opt <- nlminb(start = opt$par, objective = obj$fn, gradient = obj$gr, hessian = obj$he, lower = bounds$lower, upper = bounds$upper, control = control)
    # ce <- check_estimability(obj = obj)
    # grd[1:i,]
    # unique(names(obj$par)) # List of parameters that are "on"
    # obj$fn()
    # obj$report()$nll
    # obj$report()$ll_prior
    # obj$report()$ll_data
    # obj$report()$ll_penalty
    # obj$report()$lp_rec
    # obj$report()$lp_sel
    # obj$report()$lp_sel_old
    # obj$report()$lp_sel_par
    # obj$report()$lp_af
    # obj$report()$lp_lf
    # obj$report()$lp_cpue
    # obj$report()$lp_troll
    # obj$report()$lp_aerial
    # opt <- tryCatch({
    #   if (i == 1) { start <- obj$env$last.par.best } else { start <- opt$par }
    #   opt <- nlminb(start = start, objective = obj$fn, gradient = obj$gr, hessian = obj$he, lower = bounds$lower, upper = bounds$upper, control = control)
    #   opt
    # }, error = function(e) {
    #   cat("An error occurred:", conditionMessage(e), "\n")
    #   NA
    # })
    # if (!is.list(opt)) {
    #   opt <- tryCatch({
    #     start <- obj$env$last.par.best
    #     opt <- nlminb(start = start, objective = obj$fn, gradient = obj$gr, hessian = obj$he, lower = bounds$lower, upper = bounds$upper, control = control)
    #     opt
    #   }, error = function(e) {
    #     cat("An error occurred:", conditionMessage(e), "\n")
    #     NA
    #   })
    # }
    # if (!is.list(opt)) {
    #   opt <- tryCatch({
    #     start <- obj$env$last.par.best
    #     start[names(start) == "par_log_B0"] <- start[names(start) == "par_log_B0"] * 2
    #     opt <- nlminb(start = start, objective = obj$fn, gradient = obj$gr, hessian = obj$he, lower = bounds$lower, upper = bounds$upper, control = control)
    #     opt
    #   }, error = function(e) {
    #     cat("An error occurred:", conditionMessage(e), "\n")
    #     NA
    #   })
    # }
    obj$par <- opt$par
    print(paste0("Done: ", i, " of ", n_grid, ". Convergence: ", c("true", "false")[opt$convergence + 1]))
    grid_obj[[i]] <- obj
  }
  names(grid_obj) <- names(grid_parameters)
  # stopCluster(cl = sbt_cluster)
  return(grid_obj)
}

#' Run a grid again
#' 
#' Fit a model for each grid cell defined in \code{grid} where each grid cell has a
#' different combination of fixed parameter values.
#' 
#' @param grid a \code{list} of parameter inputs for each grid cell. This can be created using \code{get_grid}.
#' @param bounds the lower and upper parameter bounds created using the \code{get_bounds} function.
#' @param control a \code{list} defining the controls passed to \code{nlminb}.
#' @return a \code{list} of grid cells.
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom stats nlminb
#' @import foreach
#' @export
#' 
rerun_grid <- function(grid, bounds, control = list(eval.max = 10000, iter.max = 10000)) {
  n_grid <- length(grid)
  for (i in 1:n_grid) {
    obj <- grid[[i]]
    opt <- nlminb(start = obj$par, objective = obj$fn, gradient = obj$gr, hessian = obj$he, lower = bounds$lower, upper = bounds$upper, control = control)
    obj$par <- opt$par
    print(paste0("Done: ", i, " of ", n_grid, ". Convergence: ", c("true", "false")[opt$convergence + 1]))
    grid[[i]] <- obj
  }
  return(grid)
}

#' Sample from the grid
#' 
#' @param grid a \code{list} of model runs.
#' @param n_samples the number of samples from the grid.
#' @return a \code{vector} of integers.
#' @export
#' 
sample_grid <- function(grid, n_samples = 2000) {
  n_grid <- length(grid)
  gdf <- matrix(NA, nrow = n_grid, ncol = 5)
  colnames(gdf) <- c("nll", "m0", "m10", "h", "psi")
  for (i in 1:n_grid) {
    rep <- grid[[i]]$report()
    # gdf[i, 1] <- rep$nll
    gdf[i, 1] <- grid[[i]]$fn()
    gdf[i, 2] <- rep$par_m0
    gdf[i, 3] <- rep$par_m10
    gdf[i, 4] <- rep$par_h
    gdf[i, 5] <- exp(rep$par_log_psi)
  }
  
  gdf <- data.frame(gdf) %>% 
    mutate(ll = -(.data$nll - min(.data$nll)))
  
  df <- gdf %>%
    group_by(.data$h, .data$psi) %>%
    summarise(sum_exp_ll = sum(exp(.data$ll))) %>%
    full_join(gdf) %>%
    mutate(prob = exp(.data$ll) / .data$sum_exp_ll)

  grid_ints <- sample(x = 1:n_grid, size = n_samples, replace = TRUE, prob = df$prob)
  
  grid_freq <- data.frame(gdf) %>% 
    mutate(grid_ints = as.factor(1:n())) %>% 
    left_join(data.frame(table(grid_ints)), by = join_by(grid_ints)) %>%
    mutate(Freq = ifelse(is.na(.data$Freq), 0, .data$Freq)) %>%
    rename(Cell = grid_ints) %>%
    relocate(.data$Cell) %>%
    tibble()
  
  return(list(grid_cells = grid_ints, grid_freq = grid_freq))
}

#' Convert grid to an snutsfit object
#' 
#' @param data a \code{list} of model runs.
#' @param parameters a \code{list} of model runs.
#' @param grid a \code{list} of model runs.
#' @param grid_parameters a \code{list} of model runs.
#' @param grid_cells a \code{list} of model runs.
#' @return a \code{vector} of integers.
#' @export
#' 
grid_to_snutsfit <- function(data, parameters, grid, grid_parameters, grid_cells) {
  map <- get_map(parameters = parameters)
  map$par_log_psi <- NULL
  map$par_log_m0 <- NULL
  map$par_log_m10 <- NULL
  map$par_log_h <- NULL
  obj <- MakeADFun(func = cmb(sbt_model, data), parameters = parameters, map = map)
  parnames <- names(obj$par)
  n_samples <- length(grid_cells$grid_cells)
  samples <- array(NA, dim = c(n_samples, 1, length(obj$par)), dimnames = list(NULL, NULL, parnames))
  for (i in 1:n_samples) {
    j <- grid_cells$grid_cells[i]
    # obj$env$parList() alternatively could use this below?
    samples[i, 1,] <- c(grid[[j]]$par[1], # par_log_B0
                        grid_parameters[[j]]$par_log_psi, 
                        grid_parameters[[j]]$par_log_m0, 
                        grid[[j]]$par[2], # par_log_m4
                        grid_parameters[[j]]$par_log_m10,
                        grid[[j]]$par[3], # par_log_m30
                        grid_parameters[[j]]$par_log_h,
                        grid[[j]]$par[-c(1:3)])
  }
  # getS3method("print", "snutsfit")
  # Create snutsfit object
  x <- list(samples = samples, sampler_params = list(), 
            monitor = NULL, model = "RTMB", metric = "", par_names = parnames, 
            warmup = 0, thin = 1, time.total = 1, algorithm = "grid")
  class(x) <- c("snutsfit", "list")
  return(x)
}

#' Check estimability of the grid
#' 
#' @param grid a \code{list} of model runs.
#' @return a \code{list} of \code{list}s.
#' @export
#' 
check_grid <- function(grid) {
  n_grid <- length(grid)
  n_par <- length(grid[[1]]$par)
  grid_check <- vector("list", length = n_grid)
  # grid_df <- array(NA, dim = c(n_grid, n_par), dimnames = list(cell = 1:n_grid, par = names(grid[[1]]$par)))
  grid_vec <- numeric(n_grid)
  for (i in 1:n_grid) {
    grid_check[[i]] <- tryCatch({ check_estimability(obj = grid[[i]]) }, error = function(e) { NA })
    # names(grid_check[[i]])
    if (all(is.na(grid_check[[i]]))) {
      grid_vec[i] <- "Not estimable"
    } else if (length(grid_check[[i]]$WhichBad) > 0) {
      grid_vec[i] <- "Some parameters are not estimable"
    } else {
      grid_vec[i] <- "All parameters are estimable"
    }
    # grid_check[[i]] <- check_estimability(obj = grid[[i]])
    # if (is.list(grid_check[[i]])) {
    #   if (!is.null(grid_check[[i]]$BadParams)) {
    #     df <- grid_check[[i]]$BadParams %>% select(Param_check) %>% t()
    #     grid_df[i,] <- df
    #   } else {
    #     grid_df[i,] <- "OK"
    #   }
    # }
  }
  return(list(grid_check = grid_check, grid_vec = grid_vec))
}
