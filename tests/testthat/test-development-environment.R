test_that("local installs retain the locked GitHub provenance in both metadata files", {
  script <- testthat::test_path("..", "..", ".github", "scripts", "setup-r-env.R")
  skip_if_not(file.exists(script), "Development setup is not included in package builds")
  setup <- new.env(parent = globalenv())
  sys.source(script, envir = setup)

  library <- tempfile("opal-setup-test-")
  withr::defer(unlink(library, recursive = TRUE))
  package <- "opalSetupFixture"
  path <- file.path(library, package)
  dir.create(file.path(path, "Meta"), recursive = TRUE)
  description <- c(
    Package = package, Version = "1.0", Title = "Setup fixture",
    Description = "Temporary package metadata for development setup tests.",
    RemoteType = "local", RemoteUrl = "/tmp/source"
  )
  write.dcf(as.data.frame(as.list(description)), file.path(path, "DESCRIPTION"))
  saveRDS(list(DESCRIPTION = description), file.path(path, "Meta", "package.rds"))
  record <- list(
    RemoteType = "github", RemoteUsername = "fixture", RemoteRepo = "fixture",
    RemoteSha = "locked-commit", RemoteSubdir = "package"
  )

  expect_null(utils::packageDescription(package, lib.loc = library)$RemoteSha)
  setup$set_locked_provenance(path, record)
  installed <- utils::packageDescription(package, lib.loc = library)
  expect_identical(unclass(installed)[names(record)], record)
  expect_null(installed$RemoteUrl)
  expect_identical(unname(read.dcf(file.path(path, "DESCRIPTION"))[1L, "RemoteSha"]),
                   record$RemoteSha)

  setup$set_locked_provenance(path, record)
  expect_identical(utils::packageDescription(package, lib.loc = library), installed)
})
