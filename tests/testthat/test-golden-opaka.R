ref_path <- test_path("_reference", "opaka-quickstart.rds")

test_that("opaka quickstart reference exists", {
  expect_true(file.exists(ref_path))
})

test_that("opaka quickstart reference matches the current scientific-model contract", {
  ref <- readRDS(ref_path)
  expect_identical(
    ref$meta$scientific_version,
    .opal_model_scientific_version,
    label = "bump .opal_model_scientific_version and regenerate together"
  )
})

test_that("opaka quickstart active-parameter layout matches reference", {
  ref <- readRDS(ref_path)
  expect_identical(names(opaka_obj()$par), ref$par_names)
})

for (point in c("start", "mle")) {
  test_that(
    sprintf("opaka quickstart objective and report match reference (%s)", point),
    {
      skip_on_cran()
      ref <- readRDS(ref_path)[[point]]
      obj <- opaka_obj()
      tol <- reference_tolerance()

      expect_equal(obj$fn(ref$par), ref$nll, tolerance = tol, label = "nll")
      if (!is.null(ref$gr)) {
        expect_equal(
          as.vector(obj$gr(ref$par)), ref$gr,
          tolerance = tol, label = "gradient"
        )
      }

      report <- obj$report(ref$par)
      missing <- setdiff(names(ref$report), names(report))
      expect_identical(
        missing, character(0),
        label = "REPORT() elements removed or renamed since reference"
      )
      for (nm in intersect(names(ref$report), names(report))) {
        expect_equal(
          report[[nm]], ref$report[[nm]],
          tolerance = tol, label = paste0("report$", nm)
        )
      }
    }
  )
}

test_that("bundled opaka quickstart fit is compatible and reproduces its objective", {
  path <- system.file("extdata", "opaka_quickstart_fit.rds", package = "opal")
  expect_true(nzchar(path))
  fit <- suppressWarnings(read_opal_fit(path, strict = TRUE))
  expect_true(opal_fit_compatibility(fit)$compatible)
  expect_identical(names(fit$fit$opt$par), names(opaka_obj()$par))
})
