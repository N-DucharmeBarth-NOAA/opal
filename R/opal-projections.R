# Resolve posterior selection once, before entering legacy projection code.
.opal_projection_input <- function(data, object = NULL, mcmc = NULL,
                                    uncertainty = NULL, dynamics = TRUE) {
  if (inherits(data, "opal_fit")) data <- opal_from_fit(data)
  if (!inherits(data, "opal_obj")) {
    if (!is.null(uncertainty)) stop("Use uncertainty with an opal_obj.", call. = FALSE)
    return(c(.opal_model_input(data, object), list(mcmc = mcmc)))
  }
  if (!is.null(object) || !is.null(mcmc)) {
    stop("Use the objective and posterior stored in the opal_obj.", call. = FALSE)
  }
  data <- opal_build(data)
  alternatives <- if (dynamics) c("mvn", "mcmc") else c("fit", "mcmc")
  if (is.null(uncertainty)) {
    if (!is.null(data$mcmc)) stop("Choose uncertainty explicitly: ", paste(alternatives, collapse = " or "), ".", call. = FALSE)
    uncertainty <- alternatives[1L]
  }
  uncertainty <- match.arg(uncertainty, alternatives)
  if (uncertainty == "mcmc") {
    mcmc <- opal_as_tmbfit(data)
    if (length(data$random) && data$mcmc$parameter_scope != "complete") {
      stop("Random-effect projections require complete joint posterior draws.", call. = FALSE)
    }
  } else {
    if (is.null(data$fit$opt)) stop("This uncertainty source requires a fitted model.", call. = FALSE)
    if (dynamics && length(data$random)) {
      stop("MVN projections for random effects are unsupported; use complete joint MCMC draws.", call. = FALSE)
    }
  }
  c(.opal_model_input(data), list(mcmc = mcmc, uncertainty = uncertainty))
}

#' Project an Opal object and retain the result
#'
#' Runs [project_dynamics()] and stores the result, settings, and the identity
#' of its source fit or posterior. Inputs for future recruitment, selectivity,
#' and catch remain explicit scientific choices. Updating the model or
#' replacing the source fit/posterior clears stored projections.
#' @param x A fitted or sampled `opal_obj`.
#' @param uncertainty Either `"mvn"` or `"mcmc"`; required when a posterior exists.
#' @param name Name of the result in `x$derived`.
#' @param seed Optional random seed, restored after the projection.
#' @param ... Arguments passed to [project_dynamics()], including future inputs.
#' @return An updated `opal_obj` with a stored projection.
#' @export
opal_project <- function(x, uncertainty = NULL, name = "projection", seed = NULL, ...) {
  validate_opal_obj(x)
  if (!is.character(name) || length(name) != 1L || is.na(name) || !nzchar(name)) stop("Supply one result name.")
  args <- list(...)
  .opal_validate_named_list(args, "projection arguments")
  if (length(intersect(names(args), c("data", "object", "mcmc")))) stop("Projection inputs come from x.")
  .opal_obj_portable(args, "projection arguments")
  source <- .opal_projection_input(x, uncertainty = uncertainty)
  if (!is.null(seed)) {
    had_seed <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
    old_seed <- if (had_seed) get(".Random.seed", .GlobalEnv) else NULL
    on.exit(if (had_seed) assign(".Random.seed", old_seed, .GlobalEnv) else
      if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
    set.seed(seed)
  }
  result <- do.call(project_dynamics, c(source[c("data", "object", "mcmc")], args))
  x$derived[[name]] <- list(result = result, uncertainty = source$uncertainty,
    identity = if (source$uncertainty == "mcmc") .opal_check_identity(x, "mcmc") else x$identity,
    settings = c(list(seed = seed), args))
  .opal_obj_seal(x)
}
