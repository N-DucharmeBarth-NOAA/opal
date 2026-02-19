#' Adaptive likelihood profiling
#'
#' Calculate 1D likelihood profiles with respect to single parameters or more
#' generally, with respect to arbitrary linear combinations of parameters (e.g. contrasts).
#'
#' @param obj Object from MakeADFun that has been optimized.
#' @param name Name or index of a parameter to profile.
#' @param lincomb Optional linear combination of parameters to profile. By default a unit vector corresponding to name.
#' @param h Initial adaptive stepsize on parameter axis.
#' @param ytol Adjusts the range of the likelihood values.
#' @param ystep Adjusts the resolution of the likelihood profile.
#' @param maxit Max number of iterations for adaptive algorithm.
#' @param parm.range Valid parameter range.
#' @param slice Do slicing rather than profiling?
#' @param adaptive Logical; Use adaptive step size?
#' @param trace Trace progress? (TRUE, or a numeric value of 1, gives basic tracing: numeric values > 1 give more information).
#' @return a \code{vector} penalty
#' @importFrom utils tail
#' @importFrom stats nlminb
#' @importFrom Matrix Diagonal
#' @export
#'
sbtprofile <- function(obj, name, lincomb, h = 1e-04, ytol = 2, ystep = 0.1,
                       maxit = ceiling(5 * ytol / ystep), parm.range = c(-Inf, Inf),
                       slice = FALSE, adaptive = TRUE, trace = TRUE) {
  
  restore.on.exit <- c("last.par.best", "random.start", "value.best", "last.par", "inner.control", "tracemgc")
  oldvars <- sapply(restore.on.exit, get, envir = obj$env, simplify = FALSE)
  restore.oldvars <- function() {
    for (var in names(oldvars)) assign(var, oldvars[[var]], envir = obj$env)
  }
  on.exit(restore.oldvars())
  par <- obj$env$last.par.best
  if (!is.null(obj$env$random)) par <- par[-obj$env$random]
  if (missing(lincomb)) {
    if (missing(name)) stop("No 'name' or 'lincomb' specified")
    stopifnot(length(name) == 1)
    if (is.numeric(name)) {
      lincomb <- as.numeric(1:length(par) == name)
      name <- names(par)[name]
    } else if (is.character(name)) {
      if (sum(names(par) == name) != 1) stop("'name' is not unique")
      lincomb <- as.numeric(names(par) == name)
    } else {
      stop("Invalid name argument")
    }
  } else {
    if (missing(name)) name <- "parameter"
  }
  stopifnot(length(lincomb) == length(par))
  X <- Diagonal(length(lincomb))
  i <- which(lincomb != 0)[1]
  X[i, ] <- lincomb
  invX <- solve(X)
  direction <- invX[, i]
  C <- invX[, -i, drop = FALSE]
  that <- sum(lincomb * par)
  if (slice) {
    f <- function(x) {
      par <- par + x * direction
      dlist <- list()
      dlist$objective <- obj$fn(par)
      dlist$ll <- c(
        sum(obj$report()$lp_prior),
        sum(obj$report()$lp_sel),
        sum(obj$report()$lp_rec),
        sum(obj$report()$lp_penalty),
        sum(obj$report()$lp_af),
        sum(obj$report()$lp_lf),
        sum(obj$report()$lp_cpue_lf),
        sum(obj$report()$lp_cpue),
        sum(obj$report()$lp_aerial),
        sum(obj$report()$lp_troll),
        sum(obj$report()$lp_tags),
        sum(obj$report()$lp_pop),
        sum(obj$report()$lp_hsp),
        sum(obj$report()$lp_gt)
      )
      return(dlist)
    }
  } else {
    f <- function(x) {
      par <- par + x * direction
      if (length(C) == 0) {
        return(obj$fn(par))
        # return(sum(-obj$report(par)$cpue_ll))
      }
      newfn <- function(par0) {
        par <- par + as.vector(C %*% par0)
        obj$fn(par)
      }
      newgr <- function(par0) {
        par <- par + as.vector(C %*% par0)
        as.vector(obj$gr(par) %*% C)
      }
      obj$env$value.best <- Inf
      obj$env$inner.control$trace <- FALSE
      obj$env$tracemgc <- FALSE
      control <- list(step.min = 0.001)
      ans <- nlminb(start, newfn, newgr, control = control)
      start <<- ans$par
      if (trace > 0) {
        if (trace > 1) cat("Profile displacement:", x * direction, "\n")
        cat("Profile value:", ans$objective, "\n")
      }
      dlist <- list()
      dlist$objective <- ans$objective
      dlist$ll <- c(
        sum(obj$report()$lp_prior),
        sum(obj$report()$lp_sel),
        sum(obj$report()$lp_rec),
        sum(obj$report()$lp_penalty),
        sum(obj$report()$lp_af),
        sum(obj$report()$lp_lf),
        sum(obj$report()$lp_cpue_lf),
        sum(obj$report()$lp_cpue),
        sum(obj$report()$lp_aerial),
        sum(obj$report()$lp_troll),
        sum(obj$report()$lp_tags),
        sum(obj$report()$lp_pop),
        sum(obj$report()$lp_hsp),
        sum(obj$report()$lp_gt)
      )
      return(dlist)
    }
  }
  f.original <- f
  f <- function(x) {
    y <- try(f.original(x), silent = TRUE)
    if (is(y, "try-error")) y <- NA
    return(y)
  }
  start <- NULL
  
  evalAlongLine <- function(h) {
    start <<- rep(0, length(par) - 1)
    x <- 0
    ylist <- f(x)
    y <- ylist$objective
    z <- ylist$ll
    if (slice) obj$env$random.start <- expression(last.par[random])
    for (it in 1:maxit) {
      yinit <- y[1]
      xcurrent <- tail(x, 1)
      ycurrent <- tail(y, 1)
      xnext <- xcurrent + h
      if (xnext + that < parm.range[1]) {
        if (trace > 1) {
          cat("below minimum value: break\n")
        }
        break
      }
      if (parm.range[2] < xnext + that) {
        if (trace > 1) {
          cat("above maximum value: break\n")
        }
        break
      }
      ylist <- f(xnext)
      ynext <- ylist$objective
      znext <- ylist$ll
      x <- c(x, xnext)
      y <- c(y, ynext)
      z <- rbind(z, znext)
      if (is.na(ynext)) {
        if (trace > 1) {
          cat("y is NA: break\n")
        }
        break
      }
      if ((ydiff <- abs(ynext - yinit)) > ytol) {
        if (trace > 1) {
          cat(sprintf("delta y=%f > %f: break\n", ydiff, ytol))
        }
        break
      }
      if (adaptive) {
        speedMax <- ystep
        speedMin <- ifelse(ynext >= yinit, ystep / 4, ystep / 8)
        if (abs(ynext - ycurrent) > speedMax) {
          h <- h / 2
          if (trace > 1) cat(sprintf("halve step size (to %f)\n", h))
        }
        if (abs(ynext - ycurrent) < speedMin) {
          h <- h * 2
          if (trace > 1) cat(sprintf("double step size (to %f)\n", h))
        }
      }
    }
    ans <- cbind(data.frame(x = x + that, y = y), z)
    names(ans) <- c(name, "value", "Prior", "Selectivity", "Recruitment", "Penalty", "AF", "LF", "CPUE LF", "CPUE", "Aerial", "Troll", "Tag", "POP", "HSP", "GT")
    ans
  }
  
  if (trace > 1) cat("profile up\n")
  ans1 <- evalAlongLine(h)
  restore.oldvars()
  if (trace > 1) cat("profile down\n")
  ans2 <- evalAlongLine(-h)
  ans <- rbind(ans1, ans2)
  ord <- order(ans[[1]])
  ans <- ans[ord, ]
  rownames(ans) <- NULL
  class(ans) <- c("tmbprofile", class(ans))
  return(ans)
}

#' Run the dynamics
#'
#' @param obj last year of model before projection
#' @param mcmc number of projection years
#' @param name number of projection years
#' @return a \code{list} of derived values.
#' @importFrom SparseNUTS extract_samples
#' @export
#'
get_mcmc_profile <- function(obj, mcmc, name) {
  # lls <- get_posterior2(object = obj, posterior = mcmc1, pars = c("tag_ll", "cpue_ll", "sexr_ll", "lf_ll", "priors"))
  lls <- NULL
  lls$tag_ll <- lls$tag_ll %>%
    group_by(chain, iter, output) %>%
    summarise(value = -sum(value))
  lls$cpue_ll <- lls$cpue_ll %>%
    group_by(chain, iter, output) %>%
    summarise(value = -sum(value))
  lls$sexr_ll <- lls$sexr_ll %>%
    group_by(chain, iter, output) %>%
    summarise(value = -sum(value))
  lls$lf_ll <- lls$lf_ll %>%
    group_by(chain, iter, output) %>%
    summarise(value = -sum(value))
  post1 <- extract_samples(fit = mcmc, inc_lp = TRUE)[, c(name, "lp__")] %>%
    data.frame() %>%
    rename(par = 1, total = lp__) %>%
    mutate(total = -1 * total) %>%
    bind_cols(
      tag = lls$tag_ll$value, cpue = lls$cpue_ll$value,
      sexr = lls$sexr_ll$value, lf = lls$lf_ll$value, prior = -lls$priors$value
    )
  return(post1)
}

#' Run the dynamics
#'
#' @param x MLE profile
#' @param y MCMC profile
#' @param xlab number of projection years
#' @param rescale number of projection years
#' @return a \code{list} of derived values.
#' @importFrom SparseNUTS extract_samples
#' @importFrom mgcv gam
#' @export
#'
plot_profile <- function(x, y = NULL, xlab = NULL, rescale = TRUE) {
  par_name <- names(x)[1]
  if (is.null(xlab)) xlab <- par_name
  
  # if (!is.null(y)) {
  #   post2 <- y %>% pivot_longer(cols = !par)
  #   newdata <- data.frame(par = seq(min(y$par), max(y$par), length.out = 1000))
  #   par_gam <- list()
  #   par_gam[[1]] <- gam(total ~ s(par, bs = "cs"), data = y)
  #   par_gam[[2]] <- gam(tag ~ s(par, bs = "cs"), data = y)
  #   par_gam[[3]] <- gam(sexr ~ s(par, bs = "cs"), data = y)
  #   par_gam[[4]] <- gam(lf ~ s(par, bs = "cs"), data = y)
  #   par_gam[[5]] <- gam(cpue ~ s(par, bs = "cs"), data = y)
  #   par_gam[[6]] <- gam(prior ~ s(par, bs = "cs"), data = y)
  #   par_smooth <- list()
  #   for (i in 1:length(par_gam)) par_smooth[[i]] <- predict(par_gam[[i]], newdata = newdata)
  #   mcmc_min <- data.frame(
  #     name = c("total", "tag", "sexr", "lf", "cpue", "prior"),
  #     value = as.numeric(lapply(par_smooth, min))
  #   )
  #   dev_exp <- data.frame(
  #     name = c("total", "tag", "sexr", "lf", "cpue", "prior"),
  #     value = paste0("Dev. exp.=", round(100 * as.numeric(lapply(par_gam, function(x) summary(x)$dev.expl)), 1), "%")
  #   )
  # }
  
  df <- data.frame(x) %>%
    rename(par = 1, Total = value) %>%
    pivot_longer(cols = !par)# %>%
    # filter(!name %in% c("catch", "jacobian", "penalty"))
  
  if (rescale) {
    ylab <- "Rescaled log-likelihood"
    mle_prof_min <- df %>% group_by(name) %>% summarise(diff = min(.data$value))
    df <- df %>% left_join(mle_prof_min, by = join_by(name)) %>% mutate(value = .data$value - .data$diff) %>% select(-diff)
  } else { 
    ylab <- "Log-likelihood"
  }
  
  mle_max <- df %>% 
    group_by(name) %>% 
    summarise(max = max(.data$value)) %>%
    arrange(-max) %>%
    head(6)
  
  df <- df %>% 
    # right_join(mle_max, by = join_by(name)) %>% 
    # filter(max > 2) %>%
    # select(-max) %>%
    mutate(par = exp(par))
  # df <- df %>% filter(max(value) > 2)
  
  # if (!is.null(y) & rescale) {
  # mle_prof_min <- df %>%
  #   group_by(name) %>%
  #   summarise(value = min(value))
  # ll_diff <- left_join(mcmc_min, mle_prof_min, by = join_by(name)) %>%
  #   mutate(diff = value.x - value.y) %>%
  #   select(name, diff)
  # df <- df %>%
  #   left_join(ll_diff, by = join_by(name)) %>%
  #   mutate(value = value + diff)
  # }
  
  # if (!is.null(y)) {
  #   p <- ggplot(data = post2, aes(x = par, y = value, color = name)) +
  #     geom_vline(xintercept = obj$env$last.par.best[par_name], linetype = "dashed") +
  #     geom_point(color = "black", alpha = 0.25) +
  #     geom_smooth(se = FALSE, color = "black", alpha = 0.25) +
  #     geom_label(data = dev_exp, aes(x = -Inf, y = Inf, label = value), hjust = 0, vjust = 1, color = "black", label.r = unit(0, "lines")) +
  #     geom_line(data = df, linewidth = 1.5)
  # } else {
  p <- ggplot(data = df, aes(x = .data$par, y = .data$value, color = .data$name)) +
    geom_vline(xintercept = exp(obj$env$last.par.best[par_name]), linetype = "dashed") +
    geom_hline(yintercept = 2, linetype = "dashed") +
    geom_line()
  # }
  p <- p + 
    # facet_wrap(.data$name ~ ., scales = "free_y") +
    facet_wrap(.data$name ~ .) +
    labs(x = xlab, y = ylab, color = "Component")
  
  return(p)
}
