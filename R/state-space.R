#' Numbers-at-age state process contribution
#'
#' Calculates the negative log-density for a latent start-of-year
#' numbers-at-age state around the corresponding deterministic population
#' transition. States are parameterised on the log scale so they remain
#' positive during optimisation.
#'
#' @param log_number_state_a Numeric vector containing a latent log
#'   numbers-at-age state.
#' @param number_pred_a Numeric vector containing the deterministic prediction
#'   for the same state.
#' @param log_sigma_state Numeric scalar or age-specific vector of log process
#'   standard deviations.
#' @param bias_correct Logical. If \code{TRUE}, subtract
#'   \code{0.5 * sigma_state^2} from the predicted log state so the transition
#'   is mean-unbiased on the arithmetic scale.
#' @param state_floor Small positive constant added before taking logarithms.
#' @return A list containing the scalar negative log-likelihood \code{nll}, the
#'   process \code{residual_a}, the log-scale transition mean
#'   \code{mean_log_number_a}, and \code{sigma_state_a}.
#' @importFrom RTMB ADoverload dnorm
#' @export
get_state_process_nll <- function(log_number_state_a, number_pred_a,
                                  log_sigma_state, bias_correct = TRUE,
                                  state_floor = 1e-12) {
  "[<-" <- ADoverload("[<-")

  n_age <- length(number_pred_a)
  if (length(log_number_state_a) != n_age) {
    stop("'log_number_state_a' and 'number_pred_a' must have equal length.")
  }
  if (!length(log_sigma_state) %in% c(1L, n_age)) {
    stop("'log_sigma_state' must have length 1 or match the number of ages.")
  }
  if (length(state_floor) != 1L || !is.finite(state_floor) || state_floor <= 0) {
    stop("'state_floor' must be a finite positive scalar.")
  }

  sigma_state <- exp(log_sigma_state)
  sigma_state_a <- number_pred_a * 0 + sigma_state[1L]
  if (length(sigma_state) == n_age) sigma_state_a[] <- sigma_state

  mean_log_number_a <- log(number_pred_a + state_floor)
  if (isTRUE(bias_correct)) {
    mean_log_number_a <- mean_log_number_a - 0.5 * sigma_state_a^2
  }
  residual_a <- (log_number_state_a - mean_log_number_a) / sigma_state_a
  nll <- -sum(dnorm(
    x = log_number_state_a,
    mean = mean_log_number_a,
    sd = sigma_state_a,
    log = TRUE
  ))

  list(
    nll = nll,
    residual_a = residual_a,
    mean_log_number_a = mean_log_number_a,
    sigma_state_a = sigma_state_a
  )
}

#' Initialise latent numbers-at-age parameters
#'
#' Adds state-space parameters to an existing parameter list using the annual
#' states reported by a deterministic \code{do_dynamics()} or
#' \code{opal_model()} run. The first year is conditioned on the existing
#' equilibrium initialisation; rows 2 through \code{n_year + 1} initialise the
#' latent transition states.
#'
#' @param parameters Existing model parameter list.
#' @param number_ysa Numeric array with dimensions
#'   \code{[n_year + 1, n_season, n_age]}.
#' @param process_sigma Positive initial process standard deviation.
#' @param state_floor Small positive floor applied before taking logarithms.
#' @return \code{parameters} with \code{log_number_state_ya} and
#'   \code{log_sigma_state} added.
#' @export
initialize_state_space_parameters <- function(parameters, number_ysa,
                                              process_sigma = 0.1,
                                              state_floor = 1e-12) {
  state_dims <- dim(number_ysa)
  if (is.null(state_dims) || length(state_dims) != 3L || state_dims[1L] < 2L) {
    stop(paste0(
      "'number_ysa' must have dimensions [n_year + 1, n_season, n_age] ",
      "with at least two annual states."
    ))
  }
  if (length(process_sigma) != 1L || !is.finite(process_sigma) || process_sigma <= 0) {
    stop("'process_sigma' must be a finite positive scalar.")
  }
  if (length(state_floor) != 1L || !is.finite(state_floor) || state_floor <= 0) {
    stop("'state_floor' must be a finite positive scalar.")
  }

  n_year <- state_dims[1L] - 1L
  n_age <- state_dims[3L]
  state_ya <- matrix(
    number_ysa[seq.int(2L, state_dims[1L]), 1L, ],
    nrow = n_year,
    ncol = n_age
  )
  if (any(!is.finite(state_ya)) || any(state_ya < 0)) {
    stop("'number_ysa' must contain finite, non-negative annual states.")
  }

  parameters$log_number_state_ya <- log(pmax(state_ya, state_floor))
  parameters$log_sigma_state <- log(process_sigma)
  parameters
}
