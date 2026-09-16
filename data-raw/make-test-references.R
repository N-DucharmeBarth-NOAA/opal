# Generate the reviewed whole-model golden reference for the Opakapaka quickstart.
# Run from the package root with:
#   Rscript data-raw/make-test-references.R
#
# Regeneration protocol for intentional numerical changes:
# 1. Show the golden test failing on the parent reference and record moved labels.
# 2. Explain the change in NEWS.md and the pull request description.
# 3. Bump .opal_model_scientific_version in R/opal-fit.R.
# 4. Update meta$reason and regenerate this reference.
# 5. Regenerate inst/extdata/opaka_quickstart_fit.rds with
#    data-raw/generate-opaka-quickstart-fit.R.
# 6. Commit the model change, contract bump, reference, and bundled fit together.

devtools::load_all()

inputs <- opaka_quickstart_inputs()
obj <- RTMB::MakeADFun(
  func = cmb(opal_model, inputs$data),
  parameters = inputs$parameters,
  map = inputs$map,
  silent = TRUE
)
bounds <- get_bounds(obj, inputs$parameters)
control <- list(eval.max = 10000, iter.max = 10000)
opt <- nlminb(
  obj$par, obj$fn, obj$gr,
  lower = bounds$lower, upper = bounds$upper, control = control
)
opt <- nlminb(
  opt$par, obj$fn, obj$gr,
  lower = bounds$lower, upper = bounds$upper, control = control
)

snap_point <- function(par, with_gr = TRUE) {
  list(
    par = par,
    nll = obj$fn(par),
    gr = if (with_gr) as.vector(obj$gr(par)) else NULL,
    max_gr = max(abs(obj$gr(par))),
    report = obj$report(par)
  )
}

ref <- list(
  meta = list(
    created = format(Sys.time(), tz = "UTC"),
    git_sha = system("git rev-parse HEAD", intern = TRUE),
    opal = as.character(packageVersion("opal")),
    RTMB = as.character(packageVersion("RTMB")),
    RTMBdist = as.character(packageVersion("RTMBdist")),
    R = R.version.string,
    platform = R.version$platform,
    scientific_version = .opal_model_scientific_version,
    opal_model_signature = .opal_model_metadata()$signature,
    reason = "Initial reference prior to dev -> main merge"
  ),
  par_names = names(obj$par),
  start = snap_point(obj$par),
  mle = snap_point(opt$par, with_gr = FALSE)
)

output <- "tests/testthat/_reference/opaka-quickstart.rds"
dir.create(dirname(output), showWarnings = FALSE, recursive = TRUE)
saveRDS(ref, output, version = 3)
message("Wrote ", output)
