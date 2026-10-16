# Run from the repository root; .Rprofile bootstraps renv on a fresh checkout.
set_locked_provenance <- function(package_path, record) {
  path <- file.path(package_path, "DESCRIPTION")
  desc <- as.list(read.dcf(path)[1L, ])
  desc[grep("^Remote", names(desc))] <- NULL
  fields <- grep("^Remote", names(record), value = TRUE)
  desc[fields] <- record[fields]
  write.dcf(as.data.frame(desc, stringsAsFactors = FALSE), path)
  metadata_path <- file.path(package_path, "Meta", "package.rds")
  metadata <- readRDS(metadata_path)
  metadata$DESCRIPTION <- unlist(desc, use.names = TRUE)
  saveRDS(metadata, metadata_path)
}

setup_r_env <- function() {
  timed <- function(label, code) {
    message("\n== ", label, " ==")
    start <- proc.time()[["elapsed"]]
    on.exit(message(sprintf("== %s: %.1f seconds ==", label,
                            proc.time()[["elapsed"]] - start)))
    force(code)
  }

  timed("Bootstrap locked renv", {
    if (!requireNamespace("renv", quietly = TRUE))
      source("renv/activate.R")
    lock <- renv::lockfile_read("renv.lock")
    version <- lock$Packages$renv$Version
    stopifnot(length(version) == 1L, nzchar(version))
    if (as.character(utils::packageVersion("renv")) != version) {
      renv::install(paste0("renv@", version), prompt = FALSE)
      # A loaded namespace cannot be replaced safely. Restart with the newly
      # installed renv, without re-running the older activate.R bootstrap.
      status <- system2(file.path(R.home("bin"), "Rscript"),
                        c("--vanilla", ".github/scripts/setup-r-env.R",
                          commandArgs(trailingOnly = TRUE)),
                        env = paste0("R_LIBS=", shQuote(paste(.libPaths(),
                                                           collapse = .Platform$path.sep))))
      quit(status = status)
    }
    renv::load(project = getwd())
    if (tolower(Sys.getenv("RENV_CONFIG_PAK_ENABLED", "FALSE")) %in%
        c("true", "t", "1")) {
      renv::restore(packages = "pak", prompt = FALSE)
    }
  })

  phases <- c("tmb", "composition", "restore", "matrix", "verify")
  phase <- commandArgs(trailingOnly = TRUE)
  if (!length(phase)) {
    # The old workflow used separate R sessions. Preserve that isolation:
    # restoring/rebuilding a package cannot replace an already loaded DLL.
    for (phase in phases) {
      status <- system2(file.path(R.home("bin"), "Rscript"),
                        c("--vanilla", ".github/scripts/setup-r-env.R", phase),
                        env = paste0("R_LIBS=", shQuote(paste(.libPaths(),
                                                           collapse = .Platform$path.sep))))
      if (status != 0L)
        quit(status = status)
    }
    return(invisible())
  }
  stopifnot(length(phase) == 1L, phase %in% phases)

  graft_osa <- function() {
    header <- "include/contrib/OSA_multivariate_dists-main/distr.hpp"
    if (!file.exists(system.file(header, package = "TMB"))) {
      # Required by compResidual; these contrib headers are not shipped by TMB.
      TMB:::install.contrib(
        "https://github.com/vtrijoulet/OSA_multivariate_dists/archive/main.zip"
      )
    }
    stopifnot(file.exists(system.file(header, package = "TMB")))
  }

  if (phase == "tmb") timed("Restore TMB and graft OSA headers", {
    # TMB first, in its own transaction so it lands in the real library.
    # This is the short list of TMB build dependencies, not a second lockfile.
    build <- c("Rcpp", "Matrix", "RcppEigen", "cpp11",
               "BH", "RcppParallel", "rapidjsonr", "TMB")
    renv::restore(packages = build, prompt = FALSE)
    graft_osa()
  })

  if (phase == "composition") timed("Install compResidual from sources", {
    # It lives in a sub-directory of its GitHub repo, and the repo and the
    # sub-directory share a name. renv 1.2.x's graph installer collides on
    # that while staging sources:
    #   Failed to prepare 'compResidual': target file
    #   '/tmp/.../renv-graph-staging-*/compResidual' already exists
    # Only visible on a cold cache — with a warm renv cache the package is
    # linked and this code path never runs.
    #
    # So fetch the sources ourselves, at the SHA recorded in renv.lock, and
    # install from the unpacked directory.
    rec <- lock$Packages$compResidual
    stopifnot(!is.null(rec), !is.null(rec$RemoteSha))
    user <- rec$RemoteUsername
    repo <- rec$RemoteRepo
    subdir <- if (!is.null(rec$RemoteSubdir)) rec$RemoteSubdir else ""
    ref <- rec$RemoteSha
    installed <- tryCatch(suppressWarnings(utils::packageDescription(rec$Package)),
                          error = function(e) NULL)
    if (!is.list(installed) || !identical(installed$Version, rec$Version) ||
        !identical(installed$RemoteSha, ref)) {
      message(rec$Package, ": ", user, "/", repo, " @ ", ref,
              if (nzchar(subdir)) paste0(" [", subdir, "]") else "")
      tarball <- tempfile(fileext = ".tar.gz")
      exdir <- tempfile("compResidual-src-")
      dir.create(exdir)
      on.exit(unlink(c(tarball, exdir), recursive = TRUE), add = TRUE)
      download.file(
        sprintf("https://github.com/%s/%s/archive/%s.tar.gz", user, repo, ref),
        tarball, mode = "wb"
      )
      untar(tarball, exdir = exdir)
      root <- list.dirs(exdir, recursive = FALSE)
      stopifnot(length(root) == 1L)
      pkgdir <- if (nzchar(subdir)) file.path(root, subdir) else root
      stopifnot(file.exists(file.path(pkgdir, "DESCRIPTION")))

      # renv handles a plain local directory fine; keep pak out of this one
      # call so the path is deterministic regardless of the global setting.
      pak_enabled <- getOption("renv.config.pak.enabled")
      options(renv.config.pak.enabled = FALSE)
      tryCatch(renv::install(pkgdir, prompt = FALSE),
               finally = options(renv.config.pak.enabled = pak_enabled))

      # Local installs lose their GitHub identity. Retain the locked provenance
      # so status()/snapshot() and subsequent setup runs recognise the package.
      set_locked_provenance(find.package(rec$Package), rec)
      stopifnot(identical(utils::packageDescription(rec$Package)$RemoteSha, ref))
    }
  })

  if (phase == "restore") timed("Restore remaining locked packages", {
    # Everything else, from a lockfile with the compResidual record removed
    # so restore() does not undo the out-of-band installation.
    trimmed <- tempfile(fileext = ".lock")
    on.exit(unlink(trimmed), add = TRUE)
    lock$Packages$compResidual <- NULL
    renv::lockfile_write(lock, file = trimmed)
    renv::restore(lockfile = trimmed, prompt = FALSE)
  })

  if (phase == "matrix") timed("Align native matrix dependencies with R 4.6", {
    # Matrix 1.7-1 in the R 4.5 lockfile predates the R 4.6 C-API changes.
    # Its sparse Cholesky path can segfault on R 4.6. Keep the locked R 4.5
    # job, but use an explicit compatible Matrix on newer R, and rebuild
    # its native consumers against the installed headers.
    if (getRversion() >= "4.6.0") {
      renv::install("Matrix@1.7-6", prompt = FALSE)
      renv::rebuild(c("TMB", "RTMB"), recursive = FALSE,
                    type = "source", prompt = FALSE)
    }
  })

  if (phase == "verify") timed("Verify environment", {
    # Cheap insurance: if anything above reinstalled TMB, distr.hpp is gone
    # and downstream compiles fail with confusing errors much later.
    graft_osa()
    a <- Matrix::sparseMatrix(i = c(1L, 2L), j = c(1L, 2L),
                             x = c(2, 3), symmetric = TRUE)
    factor <- Matrix::Cholesky(a, LDL = FALSE, super = FALSE)
    stopifnot(isTRUE(all.equal(as.numeric(Matrix::solve(factor, c(2, 3))),
                              c(1, 1))))
    required <- c("TMB", "RTMB", "RTMBdist", "compResidual", "SparseNUTS")
    stopifnot(all(vapply(required, requireNamespace, logical(1), quietly = TRUE)))
    message("R ", getRversion(), "; renv ", utils::packageVersion("renv"),
            "; TMB ", utils::packageVersion("TMB"),
            "; RTMB ", utils::packageVersion("RTMB"),
            "; Matrix ", utils::packageVersion("Matrix"))
  })
}

if (sys.nframe() == 0L)
  setup_r_env()
