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

test_that("Codespaces downloads locked GitHub commits without forwarding credentials", {
  script <- testthat::test_path("..", "..", ".github", "scripts", "setup-r-env.R")
  skip_if_not(file.exists(script), "Development setup is not included in package builds")
  setup <- new.env(parent = globalenv())
  sys.source(script, envir = setup)
  record <- list(Source = "GitHub", RemoteUsername = "fixture",
                 RemoteRepo = "package", RemoteSha = "locked-commit")
  api <- "https://api.github.com/repos/fixture/package/tarball/locked-commit"
  archive <- "https://codeload.github.com/fixture/package/tar.gz/locked-commit"
  calls <- list()
  status <- 0L
  download <- function(url, destfile, quiet, mode, headers = NULL) {
    calls[[length(calls) + 1L]] <<- list(url = url, headers = headers)
    if (identical(url, archive)) status else 0L
  }
  downloader <- setup$locked_github_downloader(list(record), download, download)
  headers <- c(Authorization = "test-only-placeholder")

  expect_identical(downloader(api, tempfile(), headers = headers), 0L)
  expect_identical(calls, list(list(url = archive, headers = NULL)))

  calls <- list()
  status <- 1L
  expect_identical(downloader(api, tempfile(), headers = headers), 0L)
  expect_identical(calls, list(list(url = archive, headers = NULL),
                               list(url = api, headers = headers)))

  calls <- list()
  other <- "https://example.org/source.tar.gz"
  downloader(other, tempfile(), headers = headers)
  expect_identical(calls, list(list(url = other, headers = headers)))

  calls <- list()
  enterprise <- modifyList(record, list(RemoteHost = "github.example.org"))
  setup$locked_github_downloader(list(enterprise), download, download)(api, tempfile())
  expect_identical(calls, list(list(url = api, headers = NULL)))
})

test_that("API fallback uses normal renv headers without recursively invoking the hook", {
  script <- testthat::test_path("..", "..", ".github", "scripts", "setup-r-env.R")
  skip_if_not(file.exists(script), "Development setup is not included in package builds")
  setup <- new.env(parent = globalenv())
  sys.source(script, envir = setup)
  override <- function(...) stop("Recursive override")
  withr::local_options(list(renv.download.override = override))
  setup$getFromNamespace <- function(name, namespace) {
    expect_identical(name, "download")
    expect_identical(namespace, "renv")
    function(url, destfile, quiet, headers) {
      expect_null(getOption("renv.download.override"))
      expect_identical(url, "https://api.github.com/example")
    }
  }
  expect_identical(setup$download_with_renv("https://api.github.com/example",
                                           tempfile(), TRUE, "wb", NULL), 0L)
  expect_identical(getOption("renv.download.override"), override)

  setup$getFromNamespace <- function(...) function(...) stop("Download failed")
  expect_error(setup$download_with_renv("https://api.github.com/example",
                                       tempfile(), TRUE, "wb", NULL), "Download failed")
  expect_identical(getOption("renv.download.override"), override)
})
