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
