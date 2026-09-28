# Staged assessment workflow; run deliberately from an installed development opal.
# No production assessment or MCMC is launched by sourcing this file.
library(opal)
inputs <- opaka_quickstart_inputs()
assessment <- opal_obj(inputs$data, inputs$parameters, inputs$map,
                       metadata = list(stock = "Opakapaka"))
assessment <- opal_fit(assessment)
print(summary(assessment))
print(plot_cpue(assessment))

file <- tempfile(fileext = ".rds")
opal_save(assessment, file)
restored <- opal_read(file, rebuild = TRUE)
stopifnot(isTRUE(all.equal(opal_report(restored), opal_report(assessment))))
unlink(file)

# Run these steps explicitly with suitable sampling effort and scientific review:
# assessment <- opal_mcmc(assessment, seed = 42)
# assessment <- opal_check(assessment, scope = "mcmc")
# assessment <- opal_project(assessment, uncertainty = "mcmc", seed = 42,
#   n_proj = 5, n_iter = n_iter, rdev_y = future_recruitment,
#   sel_fya = future_selectivity, catch_ysf = future_catch)
# opal_save(assessment, "assessment.rds")
