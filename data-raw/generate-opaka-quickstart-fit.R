devtools::load_all()

inputs <- opaka_quickstart_inputs()
# Historical objects are configuration inputs only, never compatible runtimes.
previous <- Sys.getenv("OPAL_REBUILD_FROM")
if (nzchar(previous)) {
  old <- readRDS(previous)
  inputs <- old[c("data", "parameters", "map")]
}
control <- list(eval.max = 10000, iter.max = 10000)
obj <- RTMB::MakeADFun(
  func = cmb(opal_model, inputs$data),
  parameters = inputs$parameters,
  map = inputs$map,
  silent = TRUE
)
bounds <- get_bounds(obj, inputs$parameters)
# OPAL_REBUILD_FROM reattaches a verified optimum for report-only changes.
if (nzchar(previous)) {
  stopifnot(identical(names(obj$par), names(old$fit$opt$par)))
  opt <- old$fit$opt
  stopifnot(abs(obj$fn(opt$par) - opt$objective) < 1e-8,
            all(is.finite(obj$gr(opt$par))))
} else {
  opt <- nlminb(
    obj$par, obj$fn, obj$gr,
    lower = bounds$lower, upper = bounds$upper, control = control
  )
  opt <- nlminb(
    opt$par, obj$fn, obj$gr,
    lower = bounds$lower, upper = bounds$upper, control = control
  )
}
invisible(obj$fn(opt$par))
fit <- opal_fit(
  data = inputs$data,
  obj = obj,
  opt = opt,
  bounds = list(lower = bounds$lower, upper = bounds$upper),
  control = control,
  estimability = tryCatch(check_estimability(obj), error = identity),
  diagnostics = list(max_gradient = max(abs(obj$gr(opt$par)))),
  metadata = list(stock = "Opakapaka", model = "Fixed-effects quickstart")
)
output <- "inst/extdata/opaka_quickstart_fit.rds"
if (file.exists(output)) unlink(output)
save_opal_fit(fit, output)
read_opal_fit(output, strict = TRUE)