library(opal)

source("inst/scripts/opaka-quickstart.R")

inputs <- opaka_quickstart_inputs()
control <- list(eval.max = 10000, iter.max = 10000)
obj <- RTMB::MakeADFun(
	func = cmb(opal_model, inputs$data),
	parameters = inputs$parameters,
	map = inputs$map,
	silent = TRUE
)
bounds <- get_bounds(obj, inputs$parameters)
opt <- nlminb(
	obj$par, obj$fn, obj$gr,
	lower = bounds$lower, upper = bounds$upper, control = control
)
opt <- nlminb(
	opt$par, obj$fn, obj$gr,
	lower = bounds$lower, upper = bounds$upper, control = control
)
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