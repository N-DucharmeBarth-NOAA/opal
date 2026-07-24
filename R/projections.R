#' Project dynamics
#'
#' Forward-projects population dynamics for \code{n_proj} years using either
#' MCMC posterior draws (when \code{mcmc} is supplied) or multivariate-normal
#' (MVN) draws derived from the Hessian-based variance-covariance matrix at the
#' MLE (when \code{mcmc = NULL}).
#'
#' @param data A \code{list} of model data (as passed to \code{opal_model}).
#' @param object The RTMB AD object returned by \code{RTMB::MakeADFun}, after
#'   optimisation.
#' @param mcmc Optional. MCMC fit object returned by \code{SparseNUTS}.  When
#'   supplied, posterior draws are used for the projection.  When \code{NULL}
#'   (default), MVN draws are generated from the Hessian-derived
#'   variance-covariance matrix.
#' @param n_proj Integer. Number of projection years.
#' @param n_iter Integer. Number of iterations (posterior draws or MVN samples).
#' @param rdev_y Numeric matrix \code{[n_iter, n_proj]}.  Projected recruitment
#'   deviates (e.g., from \code{\link{project_rec_devs}}).
#' @param sel_fya Numeric array \code{[n_iter, n_fishery, n_proj, n_age]}.
#'   Projected selectivity (e.g., from \code{\link{project_selectivity}}).
#' @param catch_ysf Numeric array \code{[n_proj, n_season, n_fishery]}.
#'   Projected observed catch by year, season, and fishery.
#' @param return_hist Logical (default \code{FALSE}).  When \code{TRUE} the
#'   function returns a named list with elements \code{dyn} (the projection
#'   results) and \code{hist_sbio}, a numeric matrix
#'   \code{[n_iter, n_year + 1]} containing the full spawning-biomass
#'   trajectory from each parameter draw over the historical period.  The
#'   last column of \code{hist_sbio} corresponds to the same terminal state
#'   as the \code{proj_year = 0} bridge point in the projection, so the two
#'   can be plotted seamlessly without any join discontinuity.
#' @return When \code{return_hist = FALSE} (default): a \code{list} of length
#'   \code{n_iter}, each element being the named list returned by
#'   \code{\link{do_dynamics}} for that iteration.  When
#'   \code{return_hist = TRUE}: a named list with elements \code{dyn} and
#'   \code{hist_sbio}.
#' @importFrom SparseNUTS extract_samples
#' @importFrom stats optimHess plogis rnorm
#' @importFrom RTMB qlogis
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export
#'
project_dynamics <- function(data, object, mcmc = NULL, n_proj = 5, n_iter = 1,
                              rdev_y, sel_fya, catch_ysf, return_hist = FALSE) {

  # Input validation ----
  dr <- dim(rdev_y)
  if (length(dr) != 2) stop("'rdev_y' must be a 2-D matrix [n_iter, n_proj].")
  if (n_iter > dr[1]) stop("'n_iter' exceeds the number of rows in 'rdev_y'.")

  ds <- dim(sel_fya)
  if (length(ds) != 4) stop("'sel_fya' must be a 4-D array [n_iter, n_fishery, n_proj, n_age].")
  if (n_iter > ds[1]) stop("'n_iter' exceeds the first dimension of 'sel_fya'.")

  dc <- dim(catch_ysf)
  if (length(dc) != 3 || dc[1] != n_proj || dc[2] != data$n_season || dc[3] != data$n_fishery)
    stop("'catch_ysf' must be a 3-D array [n_proj, n_season, n_fishery].")

  # Parameter draws ----
  if (!is.null(mcmc)) {
    # MCMC path: use posterior draws
    post <- extract_samples(fit = mcmc)
    if (n_iter > nrow(post))
      stop("'n_iter' exceeds the number of available MCMC draws.")
  } else {
    # MVN path: simulate from Hessian-derived variance-covariance matrix.
    # Parameters already named par_log_* are on log scale; the Hessian from
    # obj$he() is therefore already in log-space for those.  For any positive
    # parameter NOT already on log scale (e.g. a creep rate, a rho), we apply
    # a log-transform and adjust the Hessian via the delta method so that MVN
    # draws remain in a valid range and back-transform cleanly.
    mu_vec   <- object$env$last.par.best
    par_names <- names(mu_vec)

    # Classify each free parameter's working scale for MVN simulation:
    #   par_log_* or log_* → already on log scale   (identity, jac = 1)
    #   rdev_y, init_rdev_a → normal deviates, unbounded (identity, jac = 1)
    #   par_sel            → selectivity working-scale params, may be negative
    #                        (identity, jac = 1)
    #   *rho*              → any autocorrelation param in (0,1) (logit, jac = theta*(1-theta))
    #   anything else positive                           (log,   jac = theta)
    # Deviation params (rdev_y, init_rdev_a) and selectivity working-scale params
    # (par_sel) are unbounded and can be negative, so they must NOT be
    # log-transformed even when a particular MLE value happens to be positive.
    needs_logit <- grepl("rho", par_names) & !grepl("logit", par_names)
    needs_log   <- mu_vec > 0 &
      !grepl("log_|logit", par_names) &
      !grepl("rdev",       par_names) &
      !grepl("par_sel",    par_names) &
      !needs_logit

    # Working (unconstrained) MLE vector
    mu_unc <- mu_vec
    mu_unc[needs_log]   <- log(mu_vec[needs_log])
    mu_unc[needs_logit] <- qlogis(mu_vec[needs_logit])

    # Hessian in the native (optimizer) parameter space
    H_nat <- if (!is.null(object$he)) {
      object$he()
    } else {
      stats::optimHess(par = mu_vec, fn = object$fn, gr = object$gr)
    }

    # Delta-method Jacobian: d theta_k / d phi_k for each transform.
    # H_unc[k,j] = H_nat[k,j] * jac[k] * jac[j]
    # (second-derivative correction vanishes at the MLE where gradient = 0)
    jac <- rep(1.0, length(mu_vec))
    jac[needs_log]   <- mu_vec[needs_log]                             # d exp(phi)/d phi = theta
    jac[needs_logit] <- mu_vec[needs_logit] * (1 - mu_vec[needs_logit]) # d plogis(phi)/d phi = theta*(1-theta)
    H_unc  <- H_nat * outer(jac, jac)

    Sigma_unc <- tryCatch(
      solve(H_unc),
      error = function(e) stop(
        "Hessian is singular after scale transformation; cannot invert for MVN simulation. ",
        "Run check_estimability() to diagnose non-estimable parameters."
      )
    )

    # Cholesky MVN draws in working space, then back-transform to native space
    z         <- matrix(rnorm(n_iter * length(mu_unc)), nrow = n_iter)
    draws_unc <- sweep(z %*% chol(Sigma_unc), 2, mu_unc, "+")
    mvn_draws <- draws_unc
    mvn_draws[, needs_log]   <- exp(draws_unc[, needs_log])
    mvn_draws[, needs_logit] <- plogis(draws_unc[, needs_logit])

    # Replace MVN-drawn rdev_y with their conditional expectations given the
    # scalar-parameter draws (log_B0, log_cpue_q).
    #
    # Rationale: SE(rdev_y) ≈ σ_r (constrained mainly by the prior).  Drawing
    # all 268 simultaneously inflates E[SBY] ~76 % above the MLE via Jensen's
    # inequality.  Pinning at the *unconditional* MLE ignores the negative
    # log_B0 ↔ rdev_y correlation in the Hessian (higher B0 requires lower
    # rdev_y to fit the same declining trajectory), so the marginal B0
    # uncertainty inflates the bridge-point ribbon above the sdreport
    # delta-method ribbon.
    #
    # Setting rdev_y to E[rdev_y | scalar_draw] = μ + Σ_{rdev,sc} Σ_{sc,sc}⁻¹
    # (draw − μ) preserves this negative correlation, eliminates Jensen bias
    # (individual rdev_y shifts are O(0.17) vs σ_r = 0.6), and aligns the
    # bridge-point ribbon with the sdreport historical uncertainty.
    # Recruitment-type deviations (rdev_y and, when estimated, init_rdev_a) are
    # only weakly identified (SE ≈ σ_r, constrained mainly by the prior), so the
    # conditional-expectation treatment below is applied to all of them.
    rdev_idx   <- grepl("rdev", par_names)
    scalar_idx <- which(!rdev_idx)                        # log_B0, log_cpue_q, par_sel, ...
    B_cond     <- Sigma_unc[rdev_idx, scalar_idx, drop = FALSE] %*%
                  solve(Sigma_unc[scalar_idx, scalar_idx, drop = FALSE])  # [n_rdev x 2]
    delta_sc   <- sweep(mvn_draws[, scalar_idx, drop = FALSE], 2,
                        mu_vec[scalar_idx], "-")                           # [n_iter x 2]
    mvn_draws[, rdev_idx] <-
      matrix(mu_vec[rdev_idx], nrow = n_iter, ncol = sum(rdev_idx), byrow = TRUE) +
      t(B_cond %*% t(delta_sc))                                            # [n_iter x n_rdev]
  }

  # Minimal data list for the projection call to do_dynamics() ----
  proj_data <- list(
    first_yr       = 1L,
    first_yr_catch = 1L,
    n_year         = n_proj,
    n_season       = data$n_season,
    n_fishery      = data$n_fishery,
    n_age          = data$n_age,
    catch_obs_ysf  = catch_ysf,
    catch_units_f  = data$catch_units_f
  )

  # Quantities that depend only on fixed (mapped) parameters are identical
  # for every draw, so compute them once from the MLE report.
  rep_mle         <- object$report()
  M_a_mle         <- rep_mle$M_a
  spa_mle         <- rep_mle$spawning_potential_a
  proj_weight_fya <- array(0, dim = c(data$n_fishery, n_proj, data$n_age))
  for (y in seq_len(n_proj)) proj_weight_fya[, y, ] <- rep_mle$weight_fya_mod[, data$n_year, ]

  # Projection loop ----
  dyn          <- vector("list", n_iter)
  hist_sbio    <- if (return_hist) matrix(NA_real_, n_iter, data$n_year + 1) else NULL
  if (n_iter > 1) pb <- txtProgressBar(min = 1, max = n_iter, style = 3)

  for (i in seq_len(n_iter)) {
    rep <- if (!is.null(mcmc)) {
      object$report(as.numeric(post[i, ]))
    } else {
      object$report(mvn_draws[i, ])
    }

    if (return_hist) hist_sbio[i, ] <- rep$spawning_biomass_y

    # Terminal state from this draw's full historical trajectory — reflects
    # parameter uncertainty in the initial condition.  Future process
    # uncertainty (rdev_y, sel_fya) is then added on top.
    dyn[[i]] <- do_dynamics(
      data                 = proj_data,
      parameters           = list(rdev_y = rdev_y[i, ]),
      B0                   = rep$B0,
      R0                   = rep$R0,
      alpha                = rep$alpha,
      beta                 = rep$beta,
      sigma_r              = rep$sigma_r,
      M_a                  = M_a_mle,
      spawning_potential_a = spa_mle,
      weight_fya           = proj_weight_fya,
      init_number_a        = rep$number_ysa[data$n_year + 1, 1, ],
      init_number0_a       = rep$number0_ysa[data$n_year + 1, 1, ],
      sel_fya              = sel_fya[i, , , ]
    )

    if (n_iter > 1) setTxtProgressBar(pb, i)
  }

  if (return_hist) return(list(dyn = dyn, hist_sbio = hist_sbio))
  return(dyn)
}

#' Project selectivity
#' 
#' @param data a \code{list} of parameter values.
#' @param obj a \code{list} of parameter values.
#' @param first_yr the first year sampled. Defaults to the first model year.
#' @param last_yr the last year.
#' @param n_proj the number of projection years.
#' @param n_iter the number of simulated trajectories.
#' @return An array of projected selectivity by iteration, fishery, year, and
#'   age.
#' @importFrom stats rnorm sd
#' @export
#' 
project_selectivity <- function(data, obj, first_yr = NULL, last_yr = NULL,
                                n_proj = 5, n_iter = 1) {
  if (is.null(first_yr)) first_yr <- data$first_yr
  if (is.null(last_yr)) last_yr <- data$last_yr
  proj_years <- (data$last_yr + 1):(data$last_yr + n_proj)
  samp_years <- (first_yr - data$first_yr + 1):(last_yr - data$first_yr + 1)
  sel_fya <- obj$report()$sel_fya[,samp_years,]
  
  n_fishery <- dim(sel_fya)[1]
  n_age <- dim(sel_fya)[3]
  
  sim_ifya <- array(0, dim = list(n_iter, n_fishery, n_proj, n_age), 
                    dimnames = list(iteration = 1:n_iter, fishery = 1:n_fishery, year = proj_years, age = data$age_a))
  
  removal_switch_f <- c(data$removal_switch_f, 0)
  
  for (f in seq_len(n_fishery)) {
    for (a in seq_len(n_age)) {
      for (i in seq_len(n_iter)) {
        y <- log(sel_fya[f,,a])
        if (all(is.finite(y)) & removal_switch_f[f] == 0) {
          if (sd(y) > 0) {
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
#' @param first_yr the first year sampled. Defaults to the first model year.
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
project_rec_devs <- function(data, obj, mcmc = NULL, first_yr = NULL, last_yr = NULL, n_proj = 5, n_iter = NULL,
                             max.p = 5, max.d = 5, max.q = 5, arima = TRUE) {
  if (is.null(first_yr)) first_yr <- data$first_yr
  if (!is.null(mcmc)) {
    post <- extract_samples(fit = mcmc)
    rdevs1 <- as.matrix(post[grepl("rdev_y", names(post))])
    if (!is.null(n_iter)) {
      ii <- sample(x = 1:nrow(post), size = n_iter)
      rdevs1 <- rdevs1[ii,]
    } else {
      n_iter <- nrow(post)
    }
  } else {
    if (is.null(n_iter)) n_iter <- 1
    x <- obj$env$last.par.best[names(obj$par) %in% "rdev_y"]
    rdevs1 <- matrix(x, nrow = n_iter, ncol = length(x), byrow = TRUE)
  }
  # rdevs1 may cover fewer years than the full model span when early rdev_y are
  # fixed in the map (e.g., spin-up period fixed at zero).  Derive the first
  # estimated year from the actual column count rather than assuming n_year.
  n_rdev <- ncol(rdevs1)
  first_rdev_yr <- data$last_yr - n_rdev + 1
  dimnames(rdevs1) <- list(iteration = 1:n_iter, year = first_rdev_yr:data$last_yr)

  if (is.null(last_yr)) last_yr <- data$last_yr
  samp_first_yr <- max(first_yr, first_rdev_yr)
  samp_years <- (samp_first_yr - first_rdev_yr + 1):(last_yr - first_rdev_yr + 1)
  rdevs <- rdevs1[, samp_years, drop = FALSE]
  
  proj_years <- (data$last_yr + 1):(data$last_yr + n_proj)
  sim <- matrix(NA, nrow = n_iter, ncol = n_proj, dimnames = list(iteration = 1:n_iter, year = proj_years))
  arima_pars <- matrix(NA, nrow = n_iter, ncol = 3)
  
  if (arima) {
    if (n_iter > 1) pb <- txtProgressBar(min = 1, max = n_iter, style = 3)
    for (i in seq_len(n_iter)) {
      fit <- auto.arima(y = rdevs[i,], max.p = max.p, max.d = max.d, max.q = max.q, approximation = FALSE, stepwise = FALSE, ic = "bic", trace = FALSE)
      sim[i,] <- simulate(object = fit, nsim = n_proj, future = TRUE, bootstrap = TRUE)
      arima_pars[i,] <- fit$arma[1:3]
      if (n_iter > 1) setTxtProgressBar(pb, i)
    }
  } else {
    for (i in seq_len(n_iter)) {
      fit <- ar(x = rdevs[i,], order.max = max.p)
      sim[i,] <- arima.sim(n = n_proj, list(ar = fit$ar), sd = sqrt(fit$var.pred))
      arima_pars[i, 1] <- fit$order
    }
  }
  
  return(list(rdev_y = sim, arima_pars = arima_pars))
}
