#' Resolve a biology vector to age-basis
#'
#' If the vector length matches \code{n_age}, it is returned unchanged
#' (assumed to already be on an age basis).  If the vector length matches
#' \code{n_len} (and \code{n_len != n_age}), it is converted to age using the
#' probability-of-length-at-age matrix: \code{vec_a = t(pla) \%*\% vec_l}.
#' Otherwise an informative error is thrown.
#'
#' **Note:** When \code{n_age == n_len} the vector is always treated as
#' age-basis.  If you need to pass a length-basis vector in that situation you
#' must convert it to age externally before calling the model.
#'
#' @param vec Numeric vector (length \code{n_age} or \code{n_len}).
#'   May be an AD type inside \code{opal_model()}.
#' @param n_age Integer.  Number of age classes.
#' @param n_len Integer.  Number of length bins.
#' @param pla Matrix (\code{n_len} x \code{n_age}).  Probability of length at
#'   age (columns sum to 1), as returned by \code{\link{get_pla}}.
#' @param name Character.  Name of the vector used in warning/error messages.
#' @return Numeric vector of length \code{n_age} on an age basis.
#' @importFrom RTMB ADoverload
#' @export
#'
resolve_bio_vector <- function(vec, n_age, n_len, pla, name = "vector") {
  "c" <- ADoverload("c")
  if (length(vec) == n_age) {
    if (n_age == n_len) {
      warning(sprintf(
        paste0("%s has length %d which matches both n_age and n_len. ",
               "Assuming age-basis. If this is a length-basis vector, ",
               "convert it to age externally before passing to the model."),
        name, length(vec)
      ))
    }
    return(vec)
  } else if (length(vec) == n_len) {
    return(c(t(pla) %*% vec))
  } else {
    stop(sprintf(
      "%s has length %d, which matches neither n_age (%d) nor n_len (%d)",
      name, length(vec), n_age, n_len
    ))
  }
}

#' Summarise model parameters in a table
#'
#' Builds a data.frame combining initial values, estimated values, bounds,
#' gradient, gradient check, and bounds check for all parameters. Each row
#' is one scalar element of the unlisted parameter vector.
#'
#' @param obj An RTMB AD object from `MakeADFun`.
#' @param parameters Named list of initial parameter values passed to `MakeADFun`.
#' @param map Named list of map factors passed to `MakeADFun`.
#' @param lower Numeric vector of lower bounds (length `length(obj$par)`),
#'   or `NULL` for `-Inf` everywhere.
#' @param upper Numeric vector of upper bounds (length `length(obj$par)`),
#'   or `NULL` for `Inf` everywhere.
#' @param grad_tol Threshold for gradient check: absolute gradient below this
#'   is `"OK"`, otherwise `"BAD"`. Default `1e-4`.
#' @param bnds_tol Fraction of finite bound range defining the proximity margin
#'   for the bounds check. Default `0.025` (2.5%).
#' @param include Which rows to return: `"core"` (default) returns estimated
#'   parameters excluding `rdev_y`; `"all_est"` includes `rdev_y`; `"all"`
#'   includes fixed parameters too.
#' @param digits Number of significant digits for numeric columns. Default `3`.
#'   Use `NULL` to suppress rounding.
#' @param show_map Logical. Include `map` and `fixed` columns. Default `FALSE`.
#' @return A data.frame with columns `par`, `init`, `est`, `lwr`,
#'   `upr`, `grd`, `gr_chk`, `bd_chk`, and optionally `group`, `map`
#'   and `fixed` (when `show_map = TRUE`).
#' @export
get_par_table <- function(obj, parameters, map,
                          lower    = NULL,
                          upper    = NULL,
                          grad_tol = 1e-4,
                          bnds_tol = 0.025,
                          include  = c("core", "all_est", "all"),
                          digits   = 3,
                          show_map = FALSE) {

  include <- match.arg(include)

  # --- Unlist all parameter vectors -----------------------------------------

  init_vec  <- unlist(parameters)
  est_vec   <- unlist(obj$env$parList(obj$env$last.par.best))
  n_full    <- length(init_vec)
  par_names <- names(init_vec)

  n_est <- length(obj$par)
  if (is.null(lower)) lower <- rep(-Inf, n_est)
  if (is.null(upper)) upper <- rep(Inf,  n_est)

  gr_vec <- as.vector(obj$gr())

  # --- Build group membership (which top-level list element) ----------------

  group_vec <- character(n_full)
  idx <- 1L
  for (nm in names(parameters)) {
    len <- length(unlist(parameters[[nm]]))
    group_vec[idx:(idx + len - 1L)] <- nm
    idx <- idx + len
  }

  # --- Build map value and fixed indicator ----------------------------------

  map_val_vec <- rep(NA_character_, n_full)
  fixed_vec   <- rep(FALSE, n_full)

  for (nm in names(parameters)) {
    sel <- which(group_vec == nm)
    if (nm %in% names(map)) {
      mv <- as.character(unlist(map[[nm]]))   # factor -> character
      map_val_vec[sel] <- mv
      fixed_vec[sel]   <- is.na(mv)
    }
  }

  # --- Identify which full-parameter positions map to obj$par ---------------
  # obj$par contains: freely estimated params + first occurrence of each
  # unique non-NA map level. Fixed (map = NA) params are excluded.

  in_obj_par    <- rep(FALSE, n_full)
  seen_map_vals <- character(0)
  obj_par_idx   <- integer(0)

  for (i in seq_len(n_full)) {
    if (!fixed_vec[i]) {
      mv <- map_val_vec[i]
      if (is.na(mv)) {
        # No map entry -> freely estimated, always present in obj$par
        in_obj_par[i]  <- TRUE
        obj_par_idx    <- c(obj_par_idx, i)
      } else if (!(mv %in% seen_map_vals)) {
        # First occurrence of this shared map level
        in_obj_par[i]  <- TRUE
        obj_par_idx    <- c(obj_par_idx, i)
        seen_map_vals  <- c(seen_map_vals, mv)
      }
    }
  }

  # --- Fill bounds / gradient into full-length vectors ----------------------

  lower_full <- rep(NA_real_, n_full)
  upper_full <- rep(NA_real_, n_full)
  gr_full    <- rep(NA_real_, n_full)

  lower_full[obj_par_idx] <- lower
  upper_full[obj_par_idx] <- upper
  gr_full[obj_par_idx]    <- gr_vec

  # --- Gradient check -------------------------------------------------------

  grad_check_vec <- rep(NA_character_, n_full)
  has_grad <- !is.na(gr_full)
  grad_check_vec[has_grad] <- ifelse(
    abs(gr_full[has_grad]) < grad_tol, "OK", "BAD"
  )

  # --- Bounds check ---------------------------------------------------------

  bnds_check_vec <- rep(NA_character_, n_full)
  for (i in which(!fixed_vec)) {
    lo <- lower_full[i]
    hi <- upper_full[i]
    if (is.na(lo) || is.na(hi)) next
    rng <- hi - lo
    val <- est_vec[i]
    if (is.finite(rng) && rng > 0) {
      margin <- bnds_tol * rng
      if (!is.na(val)) {
        if (is.finite(lo) && val <= lo + margin) {
          bnds_check_vec[i] <- "LO"
        } else if (is.finite(hi) && val >= hi - margin) {
          bnds_check_vec[i] <- "HI"
        } else {
          bnds_check_vec[i] <- "OK"
        }
      }
    } else {
      # At least one infinite bound: check each side independently
      if (!is.na(val)) {
        lo_flag <- is.finite(lo) && val <= lo + abs(lo) * bnds_tol + bnds_tol
        hi_flag <- is.finite(hi) && val >= hi - abs(hi) * bnds_tol - bnds_tol
        if (lo_flag) {
          bnds_check_vec[i] <- "LO"
        } else if (hi_flag) {
          bnds_check_vec[i] <- "HI"
        } else {
          bnds_check_vec[i] <- "OK"
        }
      }
    }
  }

  # --- Assemble full table --------------------------------------------------

  out <- data.frame(
    par        = par_names,
    group      = group_vec,
    init       = init_vec,
    est        = est_vec,
    lwr        = lower_full,
    upr        = upper_full,
    map        = map_val_vec,
    fixed      = fixed_vec,
    gr   = gr_full,
    gr_chk = grad_check_vec,
    bd_chk = bnds_check_vec,
    row.names  = NULL,
    stringsAsFactors = FALSE
  )

  # --- Round numeric columns -----------------------------------------------

  if (!is.null(digits)) {
    num_cols <- c("init", "est", "lwr", "upr", "gr")
    for (col in num_cols) {
      out[[col]] <- signif(out[[col]], digits = digits)
    }
  }

  # --- Filter rows based on `include` ---------------------------------------

  if (include == "core") {
    out <- out[!out$fixed & out$group != "rdev_y", , drop = FALSE]
  } else if (include == "all_est") {
    out <- out[!out$fixed, , drop = FALSE]
  }
  # "all" -> return everything

  # --- Drop map / fixed columns if not requested ----------------------------

  if (!show_map) {
    out$group <- NULL
    out$map   <- NULL
    out$fixed <- NULL
  }

  rownames(out) <- NULL
  out
}

#' Return strongly-correlated parameter pairs
#'
#' Computes the parameter correlation matrix from the (inverse) Hessian and
#' returns a data.frame of pairs whose absolute correlation exceeds
#' \code{threshold}.  Pairs are sorted by descending absolute correlation.
#' If no pairs meet the criterion the function returns an informative string
#' instead of a data.frame.
#'
#' @param obj An RTMB AD object from \code{MakeADFun} after optimisation.
#' @param h Optional pre-computed Hessian matrix (square, same dimension as
#'   \code{length(obj$par)}).  When \code{NULL} (default) the Hessian is
#'   obtained via \code{obj$he()} if available, otherwise
#'   \code{\link[stats]{optimHess}}.
#' @param threshold Absolute correlation threshold.  Only pairs with
#'   \code{|cor| > threshold} are returned.  Default \code{0.9}.
#' @param include Which parameters to consider: \code{"core"} (default)
#'   excludes \code{rdev_y} parameters; \code{"all"} includes everything in
#'   \code{obj$par}.
#' @param digits Number of significant digits applied to the correlation
#'   column.  Default \code{3}.  Use \code{NULL} to suppress rounding.
#' @return A data.frame with columns \code{par1}, \code{par2},
#'   \code{correlation} (sorted by descending \code{|correlation|}), or a
#'   character string when no pairs exceed \code{threshold}.
#' @importFrom stats ave cov2cor optimHess
#' @export
#'
get_cor_pairs <- function(obj,
                          h         = NULL,
                          threshold = 0.9,
                          include   = c("core", "all"),
                          digits    = 3) {

  include <- match.arg(include)

  # --- Build indexed parameter names ----------------------------------------
  # obj$par gives repeated bare names for vector params (e.g. "rdev_y" x 40).
  # Append a sequential index to any name that appears more than once so that
  # pairs are identifiable (e.g. "rdev_y23"), mirroring names(unlist(parameters))
  # used in get_par_table().

  raw_names  <- names(obj$par)
  within_grp <- stats::ave(seq_along(raw_names), raw_names, FUN = seq_along)
  dup_names  <- names(which(table(raw_names) > 1L))
  par_names  <- ifelse(raw_names %in% dup_names,
                       paste0(raw_names, as.integer(within_grp)),
                       raw_names)

  # --- Parameter subset selection -------------------------------------------

  keep_idx <- seq_along(par_names)
  if (include == "core") {
    keep_idx <- which(!grepl("^rdev_y", par_names))
  }

  # --- Compute or use supplied Hessian --------------------------------------

  if (!is.null(h)) {
    H <- h
  } else if (!is.null(obj$he)) {
    H <- obj$he()
  } else {
    ParHat <- extract_fixed(obj)
    H      <- stats::optimHess(par = ParHat, fn = obj$fn, gr = obj$gr)
  }

  H_sub    <- H[keep_idx, keep_idx, drop = FALSE]
  par_sub  <- par_names[keep_idx]

  # --- Covariance -> correlation ---------------------------------------------

  cov_mat <- tryCatch(
    solve(H_sub),
    error = function(e) {
      stop(paste0(
        "Hessian is singular and cannot be inverted. ",
        "Run check_estimability() to diagnose non-estimable parameters."
      ))
    }
  )

  cor_mat <- stats::cov2cor(cov_mat)

  # --- Extract upper-triangle pairs -----------------------------------------

  idx  <- which(upper.tri(cor_mat), arr.ind = TRUE)
  cors <- cor_mat[idx]

  pairs <- data.frame(
    par1        = par_sub[idx[, 1L]],
    par2        = par_sub[idx[, 2L]],
    correlation = cors,
    stringsAsFactors = FALSE
  )

  # --- Filter and sort ------------------------------------------------------

  pairs <- pairs[abs(pairs$correlation) > threshold, , drop = FALSE]

  if (nrow(pairs) == 0L) {
    return(sprintf(
      "No pairs show strong positive or negative correlation relative to %.2f threshold.",
      threshold
    ))
  }

  pairs <- pairs[order(abs(pairs$correlation), decreasing = TRUE), , drop = FALSE]

  if (!is.null(digits)) {
    pairs$correlation <- signif(pairs$correlation, digits = digits)
  }

  rownames(pairs) <- NULL
  pairs
}
