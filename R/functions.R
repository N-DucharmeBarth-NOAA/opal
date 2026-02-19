#' Helper to make closure
#' 
#' @param f a \code{vector} of midpoints.
#' @param d natural mortality at the reference size.
#' @return a \code{vector}.
#' @export
#' 
cmb <- function(f, d) function(p) f(p, d)

#' Check for identifiability of fixed effects
#'
#' Calculates the matrix of second-derivatives of the marginal likelihood
#' with respect to fixed effects, to see if any linear combinations are not estimable (i.e. cannot be
#' uniquely estimated conditional upon model structure and available data, e.g., resulting
#' in a likelihood ridge and singular, non-invertable Hessian matrix)
#'
#' @param obj The compiled object
#' @param h optional argument containing pre-computed Hessian matrix
#' @return A tagged list of the hessian and the message
#' @importFrom stats optimHess
#' @export
#' 
check_estimability <- function(obj, h) {
  ParHat <- extract_fixed(obj) # Extract fixed effects
  # Check for problems
  Gr <- obj$gr(ParHat)
  if (any(Gr > 0.01)) stop("Some gradients are high, please improve optimization and only then use `Check_Identifiable`")
  # Finite-different hessian
  List <- NULL
  if (missing(h)) {
    List[["Hess"]] <- optimHess(par = ParHat, fn = obj$fn, gr = obj$gr)
  } else {
    List[["Hess"]] <- h
  }
  # Check eigen decomposition
  List[["Eigen"]] <- eigen(List[["Hess"]])
  List[["WhichBad"]] <- which(List[["Eigen"]]$values < sqrt(.Machine$double.eps))
  # Check result
  if (length(List[["WhichBad"]]) == 0) {
    message("All parameters are estimable")
  } else {
    # Check for parameters
    RowMax <- apply(List[["Eigen"]]$vectors[,List[["WhichBad"]], drop = FALSE], MARGIN = 1, FUN = function(vec){max(abs(vec))})
    List[["BadParams"]] <- data.frame("Param" = names(obj$par), "MLE" = ParHat, "Param_check" = ifelse(RowMax > 0.1, "Bad", "OK"))
    print(List[["BadParams"]])
  }
  return(invisible(List))
}

#' Extract fixed effects
#'
#' \code{extract_fixed} extracts the best previous value of fixed effects, in a way that works for both mixed and fixed effect models
#' 
#' @param obj The compiled object
#' @return A vector of fixed-effect estimates
#' 
extract_fixed <- function(obj) {
  if (length(obj$env$random) == 0) {
    Return <- obj$env$last.par.best
  } else {
    Return <- obj$env$last.par.best[-c(obj$env$random)]
  }
  return( Return )
}

#' Positive Constraint Penalty Function
#' 
#' Ensures a value remains above a small threshold using a smooth approximation and penalty.
#' See https://github.com/kaskr/adcomp/issues/7 for discussion.
#'
#' @param x Numeric value to constrain.
#' @param eps Minimum allowable value (default 0.001).
#' @return A list with:
#' \describe{
#'   \item{new}{Transformed value.}
#'   \item{penalty}{Penalty applied if \code{x < eps}.}
#' }
#' @importFrom RTMB ADoverload logspace_add
#' @export
#'
posfun <- function(x, eps = 0.001) {
  "[<-" <- ADoverload("[<-")
  pen <- eps * (1 / (1 - (x - eps) / eps + (x - eps)^2 / eps^2 - (x - eps)^3 / eps^3 + (x - eps)^4 / eps^4 - (x - eps)^5 / eps^5))
  out <- list()
  out$new <- eps * logspace_add(x / eps, 0)
  out$penalty <- pen
  return(out)
}
