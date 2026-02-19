#' Get default initial parameter values
#' 
#' Get \code{list} of default initial parameter values.
#' 
#' @param data a \code{list} containing the data inputs.
#' @return a \code{list} of initial parameter values.
#' @importFrom RTMB qlogis
#' @export
#' 
get_parameters <- function(data) {
  
  # Get the old selectivity by age
  
  par_sels_init_i <- c(
    -4.56821478363, -3.85771953452, -3.0819114554, -2.23809099937, -1.32688067984, 
    -0.372217226998, 0.541778639307, 1.2393598106, 1.49477420269, 1.17902583073, 
    0.344374802831, -0.845471904109, -2.22931052927, -3.69982201386, -5.20612347659, 
    -6.73414421013, -0.171924122901, 0.695000877635, 0.73288296987, 0.290667868452, 
    -0.331941143169, -0.875224352921, -1.18617671285, -1.21682358711, -8.03695875624, 
    -6.50585978192, -5.79765887463, -5.81545101544, -5.47874258254, -3.00362884413, 
    0.981593912528, 1.30600103493, 2.1560337093, -0.210627387524, -2.62954236554, 
    -3.79425886058, -4.16423731152, -4.40807503186, -4.96056373235, -6.18207757277, 
    -0.652549552877, 0.500669643795, 1.18101769223, 1.38791404606, 1.1444079597, 
    0.496277433637, -0.518012422029, -1.87032784302, -3.54469077155, -5.53168273895, 
    -7.82495837692, -10.4203450156, -13.3153424186, -16.5087006191, -19.9999987938, 
    -5.87695944696, -4.7417328781, -3.72137571571, -2.81660537418, -2.02886718755, 
    -1.35916673512, -0.806158577016, -0.36443095021, -0.0238352119902, 0.229517627172, 
    0.410810757223, 0.533750980138, 0.608668100411, 0.642040500784, 0.63728361577, 
    0.595790502631, 0.51728055247, 0.399892682471, 0.241289437278, 0.0401771450023, 
    -7.02605494021, -2.16948335905, 1.11192895795, 1.27795468149, 0.0396597464904, 
    -1.69235463727, -3.61517154569, -5.54363491703)
  old_change_year_fy <- ifelse(t(as.matrix(data_csv1$sel_change_sd[,-1])) > 0, 1, 0)
  old_change_year_fy <- rbind(old_change_year_fy, old_change_year_fy[1,])
  dimnames(old_change_year_fy) <- list(fishery = 1:7, year = 1931:2022)
  old_min_age_f <- c(2, 2, 2, 8, 6, 0, 2)
  old_max_age_f <- c(17, 9, 17, 22, 25, 7, 17)
  old_end_f <- c(1, 0, 1, 1, 1, 0, 1)
  xx <- get_selectivity_v1(n_age = 31, max_age = 30, first_yr = 1931, first_yr_catch = 1952, 
                           old_min_age_f, old_max_age_f, old_end_f, old_change_year_fy,
                           par_sels_init_i, data_par1$par_sels_change_i)
  xx2 <- array(0.001, dim = c(7, 92, 31))
  xx2[1:6,,] <- xx
  xx2[7,,] <- xx[1,,] # replicate LL1 for CPUE
  xx <- xx2
  
  # Load the old selectivity by age into the new selectivity arrays
  
  par_sel <- list()
  for (f in 1:7) {
    ny <- sum(data$sel_change_year_fy[f,])
    na <- length(data$sel_min_age_f[f]:data$sel_max_age_f[f])
    yrs <- names(data$sel_change_year_fy[f,])[as.numeric(data$sel_change_year_fy[f,]) > 0]
    par_sel[[f]] <- array(0, dim = c(ny, na), dimnames = list(year = yrs, age = data$sel_min_age_f[f]:data$sel_max_age_f[f]))
  }
  for (f in 1:7) {
    fy <- data$first_yr_catch_f[f] - data$first_yr + 1
    old_change_year_fy[f, fy] <- 1
    # years
    yrs1 <- as.integer(names(old_change_year_fy[f,])[old_change_year_fy[f,] > 0])
    yrs2 <- as.integer(rownames(par_sel[[f]]))
    iy1 <- yrs1[match(yrs2, yrs1)] - data$first_yr + 1
    # ages
    a1 <- c(old_min_age_f[f]:old_max_age_f[f])
    a2 <- as.integer(colnames(par_sel[[f]]))
    ia1 <- a1[match(a2, a1)] + 1
    par_sel[[f]][,] <- log(xx[f, iy1, ia1])
  }
  par_sel[[4]] <- t(as.matrix(par_sel[[4]])) # to force as matrix
  
  Reps <- numeric(data$n_year)
  Reps[1:92] <- c(
    -0.026773339774, -0.0265371705947, -0.0264576934765, -0.0264809589676, -0.0253172755517, 
    -0.0225807621224, -0.0152193096452, 0.0228267570305, 0.149129101872, 0.244440123736, 
    -0.131170166687, -0.337394652412, 0.0139181283331, -0.0953179819162, -0.554012324153, 
    -0.717815138593, -0.890277962517, -1.24584226503, -0.905928553415, 0.0481193096045, 
    0.509883919219, 0.369178139928, 0.0200393915701, 0.0161347196642, 0.128405281831, 
    0.495241152693, 0.518621733004, 0.189229321779, 0.657620685529, -0.00309171033436, 
    -0.053090507029, 0.15123246077, -0.377796628172, -0.641485789278, 0.290862506439, 
    -0.0196666152771, 0.10662981116, 0.390241374295, 0.402286424201, -0.106742792853, 
    -0.317531521408, 0.0256285331607, 0.247039815418, 0.222881709463, 0.137140884252, 
    0.154704450454, 0.0105926667428, 0.119344869401, 0.230625290669, 0.138727457219, 
    0.0443299967093, -0.133352374547, -0.0602148772241, -0.446518702758, -0.211646422605, 
    -0.0906406842507, 0.0611385294873, 0.101199787292, 0.237065562766, 0.314697198794, 
    0.151524132403, 0.0598717817895, -0.00221903012586, 0.104054625262, 0.391133122996, 
    0.289952950947, 0.293194952246, 0.0755435407127, -0.292975637822, -0.313422770309, 
    -0.259930464052, -0.101546146365, 0.203357953914, 0.561661597369, 0.800124179485, 
    0.523613418819, 0.372302477003, 0.434354048441, 0.765522020841, 0.637879440475, 
    0.557827709277, 0.591485774246, 0.593723317885, -0.00508298501671, -0.164733210452, 
    0.0722954829899, 0.209102911964, 0.00278805353826, -0.399509917682, -0.235390996894, -0.128141681135, 0)
  
  parameters <- list(
    par_log_B0 = 16.19836, 
    par_log_m0 = log(0.4), 
    par_log_m4 = log(0.1670507),
    par_log_m10 = log(0.065), 
    par_log_m30 = log(0.45741),
    par_log_h = log(0.72),
    par_log_sigma_r = log(0.6), 
    par_log_cpue_q = -0.02033773, 
    par_cpue_creep = 0.005,
    par_log_cpue_sigma = log(0.2), 
    par_log_cpue_omega = log(1),
    
    par_log_af_alpha = c(1.5, 1),
    par_log_lf_alpha = c(1.7, 2.5, 1.5, 1, 1.8),
    
    par_sel_rho_y = c(0.7, 0.7, 0.5, 0.7, 0.5, 0.5, 0.7),
    par_sel_rho_a = c(0.9, 0.9, 0.5, 0.9, 0.9, 0.5, 0.9),
    par_log_sel_sigma = log(c(0.31, 0.25, 0.75, 0.25, 0.38, 1.13, 0.31)),
    par_log_sel_1 = par_sel[[1]], 
    par_log_sel_2 = par_sel[[2]],
    par_log_sel_3 = par_sel[[3]], 
    par_log_sel_4 = par_sel[[4]],
    par_log_sel_5 = par_sel[[5]], 
    par_log_sel_6 = par_sel[[6]],
    par_log_sel_7 = par_sel[[7]], 
    par_rdev_y = Reps
  )
  
  return(parameters)
}

#' Get default parameter mapping
#' 
#' Get a default parameter mapping. Parameter mapping is used by \code{MakeADFun} 
#' to turn parameters on/off.
#' 
#' @param parameters a \code{list} containing the initial parameter values to be passed to \code{MakeADFun}.
#' @return a named \code{list} of parameter mapping.
#' @export
#' 
get_map <- function(parameters) {
  map <- list()
  map[["par_log_psi"]] <- factor(NA)
  map[["par_log_m0"]] <- factor(NA)
  map[["par_log_m10"]] <- factor(NA)
  map[["par_log_h"]] <- factor(NA)
  map[["par_log_sigma_r"]] <- factor(NA)
  map[["par_log_cpue_sigma"]] <- factor(NA)
  map[["par_log_cpue_omega"]] <- factor(NA)
  map[["par_cpue_creep"]] <- factor(NA)
  map[["par_log_af_alpha"]] <- factor(rep(NA, 2))
  map[["par_log_lf_alpha"]] <- factor(rep(NA, 5))
  map[["par_sel_rho_y"]] <- factor(rep(NA, length(parameters$par_sel_rho_y)))
  map[["par_sel_rho_a"]] <- factor(rep(NA, length(parameters$par_sel_rho_a)))
  map[["par_log_sel_sigma"]] <- factor(rep(NA, length(parameters$par_log_sel_sigma)))
  map[["par_log_sel_4"]] <- factor(matrix(NA, nrow = nrow(parameters$par_log_sel_4), ncol = ncol(parameters$par_log_sel_4)))
  # map[["par_rec_dev_y"]] <- factor(rep(NA, length(parameters$par_rdev_y)))
  return(map)
}

#' Get default parameter bounds
#' 
#' Get \code{data.frame} of default parameter bounds.
#' 
#' @param obj a \code{list} specifying the AD object created using the \code{MakeADFun} function.
#' @param parameters a \code{list} specifying the AD object created using the \code{MakeADFun} function.
#' @return a \code{data.frame} of parameter bounds.
#' @importFrom RTMB qlogis
#' @export
#' 
get_bounds <- function(obj, parameters) {
  
  Lwr <- rep(-Inf, length(obj$par))
  Upr <- rep(Inf, length(obj$par))
  
  Lwr[grep("par_log_psi", names(obj$par))] <- log(0.5)
  Upr[grep("par_log_psi", names(obj$par))] <- log(3)
  # These were the old M bounds
  # Lwr[grep("par_log_m0", names(obj$par))] <- log(0.2)
  # Upr[grep("par_log_m0", names(obj$par))] <- log(0.55)
  # Lwr[grep("par_log_m4", names(obj$par))] <- parameters$par_log_m10
  # Upr[grep("par_log_m4", names(obj$par))] <- log(0.333 * exp(parameters$par_log_m10) + 0.667 * exp(parameters$par_log_m0))
  # Lwr[grep("par_log_m10", names(obj$par))] <- log(0.029)
  # Upr[grep("par_log_m10", names(obj$par))] <- log(0.21)
  # Lwr[grep("par_log_m30", names(obj$par))] <- log(0.2)
  # Upr[grep("par_log_m30", names(obj$par))] <- log(0.7)
  Lwr[grep("par_log_m0", names(obj$par))] <- log(1e-6)
  Upr[grep("par_log_m0", names(obj$par))] <- log(1)
  Lwr[grep("par_log_m4", names(obj$par))] <- log(1e-6)
  Upr[grep("par_log_m4", names(obj$par))] <- log(1)
  Lwr[grep("par_log_m10", names(obj$par))] <- log(1e-6)
  Upr[grep("par_log_m10", names(obj$par))] <- log(1)
  Lwr[grep("par_log_m30", names(obj$par))] <- log(1e-6)
  Upr[grep("par_log_m30", names(obj$par))] <- log(1)
  
  # Lwr[grep("par_log_cpue_tau", names(obj$par))] <- log(0.20)
  # Upr[grep("par_log_cpue_tau", names(obj$par))] <- log(0.20)
  Lwr[grep("par_log_sigma_r", names(obj$par))] <- log(0.1)
  Upr[grep("par_log_sigma_r", names(obj$par))] <- log(2.0)
  Lwr[grep("par_log_h", names(obj$par))] <- log(0.21)
  Upr[grep("par_log_h", names(obj$par))] <- log(1.0)
  Lwr[grep("par_rdev_y", names(obj$par))] <- rep(-5, length(parameters$par_rdev_y))
  Upr[grep("par_rdev_y", names(obj$par))] <- rep(5, length(parameters$par_rdev_y))
  
  check_bounds(opt = obj, lower = Lwr, upper = Upr)
  
  df <- data.frame(parameter = names(obj$par), init = obj$par, lower = Lwr, upper = Upr)
  
  return(df)
}

#' Check if parameters are up against the bounds
#' 
#' A \code{data.frame} containing the parameters that are against their lower or 
#' upper bound.
#' 
#' @param opt an optimized TMB object.
#' @param lower a vector of lower bounds.
#' @param upper a vector of upper bounds.
#' @return a \code{data.frame}.
#' @export
#' 
check_bounds <- function(opt, lower, upper) {
  df <- data.frame(par = names(opt$par), lb = lower, value = opt$par, ub = upper) %>%
    mutate(index = 1:n())
  rownames(df) <- NULL
  ilb <- which(df$value <= df$lower)
  iub <- which(df$value >= df$upper)
  return(df[c(ilb, iub), ])
}
