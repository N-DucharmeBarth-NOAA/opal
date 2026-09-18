.opal_test_cache <- new.env(parent = emptyenv())

opaka_inputs <- function() {
  if (is.null(.opal_test_cache$opaka)) {
    .opal_test_cache$opaka <- opaka_quickstart_inputs()
  }
  .opal_test_cache$opaka
}

opaka_obj <- function(inputs = opaka_inputs(), data = inputs$data,
                      parameters = inputs$parameters, map = inputs$map, ...) {
  suppressWarnings(RTMB::MakeADFun(
    func = cmb(opal_model, data), parameters = parameters,
    map = map, silent = TRUE, ...
  ))
}

fit_quickstart <- function(obj, parameters, start = obj$par,
                           control = list(eval.max = 10000, iter.max = 10000)) {
  bounds <- get_bounds(obj, parameters)
  opt <- nlminb(start, obj$fn, obj$gr,
                lower = bounds$lower, upper = bounds$upper, control = control)
  opt <- nlminb(opt$par, obj$fn, obj$gr,
                lower = bounds$lower, upper = bounds$upper, control = control)
  invisible(obj$fn(opt$par))
  opt
}

free_in_map <- function(map, parameters, free) {
  for (nm in intersect(names(free), names(parameters))) {
    if (is.null(map[[nm]])) next
    n <- length(parameters[[nm]])
    idx <- if (isTRUE(free[[nm]])) seq_len(n) else as.integer(free[[nm]])
    lv <- as.integer(as.character(map[[nm]]))
    new <- idx[is.na(lv[idx])]
    lv[new] <- max(c(0L, lv), na.rm = TRUE) + seq_along(new)
    map[[nm]] <- factor(lv)
  }
  map
}