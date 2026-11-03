# Opt-in performance measurements. Nothing in this helper runs at source time.
benchmark_root <- function() {
  root <- Sys.getenv("OPAL_BENCHMARK_ROOT", "")
  if (nzchar(root)) return(normalizePath(root, mustWork = TRUE))
  root <- normalizePath(getwd(), mustWork = TRUE)
  repeat {
    if (file.exists(file.path(root, "DESCRIPTION")) &&
        file.exists(file.path(root, "tests", "testthat", "helper-benchmark.R"))) {
      return(root)
    }
    parent <- dirname(root)
    if (identical(parent, root)) stop("Set OPAL_BENCHMARK_ROOT to the opal checkout.")
    root <- parent
  }
}

benchmark_cases <- function() {
  cases <- c(
    paste0("opaka_", c("makeadfun", "fn", "gr", "fit", "report", "sdreport")),
    paste0("bet_", c("makeadfun", "fn", "gr", "trajectory", "fit", "sdreport")),
    paste0("io_yft_", c("makeadfun", "fn", "gr", "trajectory", "fit_evals")),
    "opal_model_plain"
  )
  selection <- Sys.getenv("OPAL_BENCHMARK_CASES", "")
  if (!nzchar(selection)) {
    if (identical(Sys.getenv("OPAL_BENCHMARK_FULL_FITS"), "false")) {
      cases <- setdiff(cases, "io_yft_fit_evals")
    }
    return(cases)
  }
  selected <- trimws(strsplit(selection, ",", fixed = TRUE)[[1L]])
  if (any(!selected %in% cases)) {
    stop("Unknown benchmark case(s): ", paste(setdiff(selected, cases), collapse = ", "))
  }
  cases[cases %in% selected]
}

benchmark_timing_cases <- function(profile = Sys.getenv("OPAL_BENCHMARK_TIMING_PROFILE", "full")) {
  if (!is.character(profile) || length(profile) != 1L || is.na(profile) ||
      !profile %in% c("smoke", "full")) {
    stop("Timing profile must be smoke or full.")
  }
  cases <- setdiff(benchmark_cases(), "io_yft_fit_evals")
  if (profile == "smoke") {
    cases <- intersect(cases, c("opaka_fn", "opaka_gr", "bet_fn", "bet_gr",
                                "io_yft_fn", "io_yft_gr"))
  }
  if (!length(cases)) stop("No timing cases selected for this profile.")
  cases
}

benchmark_inputs <- function(model) {
  if (identical(model, "opaka")) return(opal::opaka_quickstart_inputs())
  if (identical(model, "io_yft")) {
    env <- new.env()
    for (name in c("io_yft_data", "io_yft_parameters", "io_yft_map")) {
      utils::data(list = name, package = "opal", envir = env)
    }
    return(list(data = env$io_yft_data, parameters = env$io_yft_parameters,
                map = env$io_yft_map))
  }
  if (!identical(model, "bet")) stop("Unknown benchmark model: ", model)
  env <- new.env()
  for (name in c("opal_baseline_data", "opal_baseline_parameters", "opal_baseline_map")) {
    utils::data(list = name, package = "opal", envir = env)
  }
  list(data = env$opal_baseline_data, parameters = env$opal_baseline_parameters,
       map = env$opal_baseline_map)
}

benchmark_objective <- function(inputs) {
  suppressWarnings(RTMB::MakeADFun(
    opal::cmb(opal::opal_model, inputs$data), inputs$parameters,
    map = inputs$map, silent = TRUE
  ))
}

benchmark_fit <- function(obj, parameters, record = FALSE) {
  bounds <- opal::get_bounds(obj, parameters)
  counts <- c(fn_evals = 0L, gr_evals = 0L, iterations = 0L)
  events <- list()
  wrap <- function(kind) {
    force(kind)
    function(par) {
      key <- paste0(kind, "_evals")
      counts[key] <<- counts[key] + 1L
      value <- obj[[kind]](par)
      if (record) {
        events[[length(events) + 1L]] <<- list(kind = kind, par = par, value = value)
        if (length(events) %% 500L == 0L) message("Recorded ", length(events), " fn/gr calls")
      }
      value
    }
  }
  fn <- wrap("fn")
  gr <- wrap("gr")
  start <- obj$par
  passes <- vector("list", 2L)
  for (pass in seq_len(2L)) {
    opt <- stats::nlminb(start, fn, gr, lower = bounds$lower, upper = bounds$upper,
                         control = list(eval.max = 10000L, iter.max = 10000L))
    if (opt$convergence != 0L) {
      stop("Benchmark fit pass ", pass, " did not converge: ", opt$message)
    }
    counts["iterations"] <- counts["iterations"] + opt$iterations
    passes[[pass]] <- opt
    start <- opt$par
  }
  if (!is.finite(opt$objective) || any(!is.finite(gr(opt$par)))) {
    stop("Non-finite objective or gradient at the benchmark solution.")
  }
  list(opt = opt, counts = counts, passes = passes, events = events)
}

benchmark_read_trajectory <- function(model, trajectory_dir, inputs) {
  path <- file.path(trajectory_dir, paste0(model, ".rds"))
  if (!file.exists(path)) {
    stop("Missing benchmark trajectory: ", path,
         ". Run Rscript dev/update-benchmark-trajectories.R explicitly.")
  }
  trajectory <- readRDS(path)
  if (!identical(trajectory$inputs, inputs)) {
    stop("Benchmark trajectory inputs changed; explicitly regenerate ", path)
  }
  if (!length(trajectory$events) ||
      any(!vapply(trajectory$events, function(x) {
        x$kind %in% c("fn", "gr") && is.numeric(x$par) && all(is.finite(x$par))
      }, logical(1)))) stop("Invalid benchmark trajectory: ", path)
  trajectory
}

benchmark_sample_trajectory <- function(trajectory, per_kind = 50L) {
  if (length(per_kind) != 1L || !is.finite(per_kind) || per_kind < 1L ||
      per_kind != as.integer(per_kind)) stop("per_kind must be a positive integer.")
  events <- trajectory$reference_events
  if (is.null(events)) events <- trajectory$events
  kinds <- vapply(events, `[[`, "", "kind")
  indices <- sort(unique(unlist(lapply(c("fn", "gr"), function(kind) {
    index <- which(kinds == kind)
    if (!length(index)) return(integer())
    index[unique(round(seq(1, length(index), length.out = min(per_kind, length(index)))))]
  }))))
  trajectory$reference_events <- events
  trajectory$events <- events[indices]
  trajectory$replay_indices <- indices
  trajectory
}

benchmark_validate_trajectory <- function(trajectory) {
  sequences <- list(trajectory$events)
  if (!is.null(trajectory$reference_events)) {
    sequences <- c(list(trajectory$reference_events), sequences)
  }
  for (events in sequences) {
    obj <- benchmark_objective(trajectory$inputs)
    for (i in seq_along(events)) {
      event <- events[[i]]
      value <- obj[[event$kind]](event$par)
      if (!identical(value, event$value)) {
        stop("Reference trajectory replay differs at event ", i, " (", event$kind, ").")
      }
    }
  }
  invisible(TRUE)
}

benchmark_operation <- function(case, trajectory_dir =
                                  file.path(benchmark_root(), "tests/testthat/_benchmark/trajectories")) {
  plain <- identical(case, "opal_model_plain")
  model <- if (plain) "opaka" else sub("_(makeadfun|fn|gr|trajectory|fit|report|sdreport|fit_evals)$", "", case)
  operation <- if (plain) "plain" else substring(case, nchar(model) + 2L)
  if (!operation %in% c("makeadfun", "fn", "gr", "trajectory", "fit",
                        "report", "sdreport", "fit_evals", "plain")) {
    stop("Unknown benchmark case: ", case)
  }
  inputs <- benchmark_inputs(model)
  obj <- NULL
  fit <- NULL
  trajectory <- NULL
  if (!operation %in% c("makeadfun", "plain")) obj <- benchmark_objective(inputs)
  if (operation %in% c("trajectory", "sdreport", "report")) {
    trajectory <- benchmark_read_trajectory(model, trajectory_dir, inputs)
    if (operation != "trajectory") invisible(obj$fn(trajectory$opt$par))
  }
  reference_par <- if (!is.null(obj)) obj$par else NULL
  direction <- -1
  next_par <- function() {
    direction <<- -direction
    reference_par + direction * 1e-10
  }
  run <- switch(operation,
    makeadfun = function() benchmark_objective(inputs),
    fn = function() obj$fn(next_par()),
    gr = function() obj$gr(next_par()),
    fit = function() {
      fit <<- benchmark_fit(obj, inputs$parameters)
      fit$opt
    },
    fit_evals = function() {
      fit <<- benchmark_fit(obj, inputs$parameters)
      fit$counts
    },
    report = function() obj$report(trajectory$opt$par),
    sdreport = function() RTMB::sdreport(obj, par.fixed = trajectory$opt$par),
    trajectory = function() {
      value <- NULL
      for (event in trajectory$events) value <- obj[[event$kind]](event$par)
      value
    },
    plain = function() opal::opal_model(inputs$parameters, inputs$data)
  )
  attr(run, "benchmark_state") <- function() list(obj = obj, fit = fit, inputs = inputs)
  run
}

benchmark_hwm <- function() {
  if (!file.exists("/proc/self/status")) return(NA_real_)
  line <- grep("^VmHWM:", readLines("/proc/self/status", warn = FALSE), value = TRUE)
  if (!length(line)) return(NA_real_)
  as.numeric(strsplit(trimws(sub("^VmHWM:", "", line)), "[[:space:]]+")[[1L]][1L]) * 1024
}

benchmark_load <- function(library = NULL) {
  if (!is.null(library)) {
    library <- normalizePath(library, mustWork = TRUE)
    if (!file.exists(file.path(library, "opal", "DESCRIPTION"))) {
      stop("No installed opal package in library: ", library)
    }
    .libPaths(c(library, .libPaths()))
  }
  if (requireNamespace("opal", quietly = TRUE)) return(invisible("installed"))
  if (!is.null(library)) stop("Cannot load opal from the supplied library.")
  if (!requireNamespace("pkgload", quietly = TRUE)) stop("Install opal or pkgload first.")
  pkgload::load_all(benchmark_root(), quiet = TRUE)
  invisible("pkgload")
}

benchmark_child <- function(case = NULL, library = NULL, iterations = 20L, idle = FALSE) {
  mode <- benchmark_load(library)
  if (!requireNamespace("bench", quietly = TRUE)) stop("Install the suggested bench package.")
  gc()
  if (idle) {
    hwm <- benchmark_hwm()
    return(list(hwm = hwm, metadata = benchmark_metadata(), load_mode = mode))
  }
  run <- benchmark_operation(case)
  result <- NULL
  measured_run <- function() {
    result <<- run()
    invisible(result)
  }
  # The allocation probe already warms long fits and tape construction.
  if (!grepl("_(makeadfun|fit|fit_evals)$", case)) measured_run()
  gc()
  if (grepl("_fit_evals$", case)) {
    # Count-only fits run once, not once for allocation plus N timing repeats.
    # Retain their one-shot timing and memory for the weekly trend artifacts.
    allocation <- profmem::profmem(timing <- bench::system_time(measured_run()))
    measured <- list(median = timing[["real"]],
                     mem_alloc = sum(allocation$bytes, na.rm = TRUE))
  } else {
    # bench runs a memory allocation probe and timed iterations separately.
    fixed_iterations <- if (grepl("_(makeadfun|fit|trajectory|sdreport)$", case)) {
      min(iterations, 3L)
    } else {
      NULL
    }
    measured <- bench::mark(measured_run(), iterations = fixed_iterations,
                            min_iterations = max(20L, iterations), min_time = 1,
                            check = FALSE, filter_gc = TRUE, memory = TRUE)
  }
  peak <- benchmark_hwm()
  state <- attr(run, "benchmark_state")()
  obj <- state$obj
  if (grepl("_makeadfun$", case)) obj <- result
  values <- c(median_seconds = as.numeric(measured$median),
              mem_alloc_bytes = as.numeric(measured$mem_alloc),
              structural_bytes = as.numeric(utils::object.size(result)),
              peak_rss_bytes = peak)
  if (!is.null(obj)) {
    # RTMB exposes these native counters; object.size omits external tape memory.
    info <- .Call(get("InfoADFunObject", asNamespace("RTMB")), obj$env$ADFun$ptr)
    values <- c(values, tape_opstack_size = info$opstack_size,
                tape_values_size = info$values_size, tape_inputs_size = info$inputs_size)
  }
  if (!is.null(state$fit)) values <- c(values, state$fit$counts)
  units <- ifelse(grepl("_bytes$", names(values)), "bytes",
                  ifelse(names(values) == "median_seconds", "seconds", "count"))
  data.frame(case = case, metric = names(values), value = unname(values),
             units = units, stringsAsFactors = FALSE)
}

run_benchmarks <- function(cases = benchmark_cases(), library = NULL, iterations = 20L,
                           idle = NULL) {
  if (length(iterations) != 1L || !is.finite(iterations) || iterations < 1 ||
      iterations != as.integer(iterations)) stop("iterations must be a positive integer.")
  if (!length(cases) || anyDuplicated(cases)) stop("Supply unique benchmark cases.")
  if (!requireNamespace("callr", quietly = TRUE)) stop("Install the suggested callr package.")
  root <- benchmark_root()
  scratch <- tempdir()
  child <- function(case, idle) {
    callr::r(function(root, case, library, iterations, idle) {
      source(file.path(root, "tests/testthat/helper-benchmark.R"), local = globalenv())
      benchmark_child(case, library, iterations, idle)
    }, args = list(root, case, library, as.integer(iterations), idle),
    libpath = .libPaths(),
    env = c(OPAL_BENCHMARK_ROOT = root, TMPDIR = scratch, OMP_NUM_THREADS = "1",
            OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1"),
    system_profile = FALSE, user_profile = FALSE)
  }
  if (is.null(idle)) idle <- child(NULL, TRUE)
  rows <- lapply(cases, function(case) {
    message("Benchmark: ", case)
    result <- child(case, FALSE)
    peak <- result$metric == "peak_rss_bytes"
    result$value[peak] <- pmax(0, result$value[peak] - idle$hwm)
    result
  })
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  attr(result, "idle_hwm_bytes") <- idle$hwm
  attr(result, "benchmark_idle") <- idle
  attr(result, "benchmark_metadata") <- c(idle$metadata,
                                         list(load_mode = idle$load_mode,
                                              iterations = iterations, cases = cases))
  result
}

benchmark_tolerance <- function() {
  value <- suppressWarnings(as.numeric(Sys.getenv("OPAL_BENCHMARK_TOL", "0.10")))
  if (length(value) != 1L || !is.finite(value) || value < 0) {
    stop("OPAL_BENCHMARK_TOL must be a finite non-negative number")
  }
  value
}

compare_to_baseline <- function(observed, baseline, tolerance = benchmark_tolerance()) {
  if (is.list(baseline) && !is.data.frame(baseline)) baseline <- baseline$results
  required <- c("case", "metric", "value", "units")
  for (x in list(observed, baseline)) {
    if (!all(required %in% names(x))) stop("Benchmark results require case, metric, value, units.")
    if (anyDuplicated(x[c("case", "metric")])) stop("Duplicate benchmark case/metric.")
  }
  if (length(tolerance) != 1L || !is.finite(tolerance) || tolerance < 0) stop("Invalid tolerance.")
  names(observed)[names(observed) == "value"] <- "observed"
  names(baseline)[names(baseline) == "value"] <- "baseline"
  # Compare every observed row. Missing stored metrics cannot silently pass.
  result <- merge(observed, baseline, by = c("case", "metric"), all.x = TRUE,
                  suffixes = c(".observed", ".baseline"), sort = FALSE)
  result$ratio <- ifelse(result$baseline == 0,
                         ifelse(result$observed == 0, 1, Inf),
                         result$observed / result$baseline)
  valid <- is.finite(result$baseline) & is.finite(result$observed) &
    result$baseline >= 0 & result$observed >= 0 &
    result$units.observed == result$units.baseline
  # Wall time is deliberately recorded, not gated in a single measurement.
  peak_noise <- result$metric == "peak_rss_bytes" &
    result$observed - result$baseline < 5 * 1024^2
  result$pass <- valid & (result$metric == "median_seconds" |
                           peak_noise |
                           result$observed <= result$baseline * (1 + tolerance))
  result$pass[is.na(result$pass)] <- FALSE
  improved <- valid & result$ratio < 1 - tolerance
  if (any(improved, na.rm = TRUE)) {
    message("Performance improved beyond tolerance; refresh the benchmark baseline: ",
            paste(paste(result$case[which(improved)], result$metric[which(improved)]),
                  collapse = ", "))
  }
  result[c("case", "metric", "baseline", "observed", "ratio", "pass")]
}

benchmark_metadata <- function() {
  root <- benchmark_root()
  benchmark_load()
  git <- function(...) {
    tryCatch(system2("git", c("-C", shQuote(root), ...), stdout = TRUE,
                     stderr = FALSE), error = function(e) NA_character_)
  }
  packages <- c("opal", "RTMB", "TMB", "bench", "callr", "pkgload")
  versions <- vapply(packages, function(p) {
    if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p)) else NA_character_
  }, character(1))
  list(schema_version = 1L, timestamp_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
       git_sha = git("rev-parse", "HEAD"), git_dirty = git("status", "--porcelain"),
       R = R.version.string, platform = R.version$platform, system = Sys.info(),
       scientific_version = get(".opal_model_scientific_version", asNamespace("opal")),
       packages = versions, opal_path = if (requireNamespace("opal", quietly = TRUE)) find.package("opal") else NA_character_,
       runner = Sys.getenv(c("RUNNER_OS", "RUNNER_ARCH", "RUNNER_NAME", "GITHUB_RUN_ID")),
       workflow = list(
         run_url = if (nzchar(Sys.getenv("GITHUB_RUN_ID"))) paste0(
           Sys.getenv("GITHUB_SERVER_URL"), "/", Sys.getenv("GITHUB_REPOSITORY"),
           "/actions/runs/", Sys.getenv("GITHUB_RUN_ID")) else NA_character_,
         attempt = Sys.getenv("GITHUB_RUN_ATTEMPT"), job = Sys.getenv("GITHUB_JOB"),
         timing_profile = Sys.getenv("OPAL_BENCHMARK_TIMING_PROFILE"),
         rounds = Sys.getenv("OPAL_BENCHMARK_ROUNDS")),
       threads = Sys.getenv(c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")),
       cpu = if (file.exists("/proc/cpuinfo")) unique(grep("^model name", readLines("/proc/cpuinfo"), value = TRUE)) else NA_character_,
       peak_rss_method = "Linux VmHWM minus a fresh idle process loading opal and bench; 5 MiB comparison floor",
       fit_evals_method = "One two-pass nlminb fit; fn/gr wrappers count actual calls; iterations summed across passes",
       source_md5 = tools::md5sum(c(file.path(root, "tests/testthat/helper-benchmark.R"),
                                  list.files(file.path(root, "data"), pattern = "\\.rda$",
                                             full.names = TRUE))))
}

write_benchmark_results <- function(results, out) {
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(results, file.path(out, "results.csv"), row.names = FALSE)
  metadata <- attr(results, "benchmark_metadata")
  if (is.null(metadata)) metadata <- benchmark_metadata()
  metadata$idle_hwm_bytes <- attr(results, "idle_hwm_bytes")
  saveRDS(list(results = results, metadata = metadata), file.path(out, "results.rds"))
  formatted <- as.data.frame(lapply(results, function(x) {
    if (is.numeric(x)) format(x, digits = 9, trim = TRUE) else as.character(x)
  }), stringsAsFactors = FALSE)
  rows <- if (nrow(formatted)) apply(formatted, 1L, paste, collapse = " | ") else character()
  summary <- c("# Benchmark results", "", paste("Measured:", metadata$timestamp_utc),
               paste("Commit:", paste(metadata$git_sha, collapse = " ")), "",
               paste0("| ", paste(names(results), collapse = " | "), " |"),
               paste0("| ", paste(rep("---", ncol(results)), collapse = " | "), " |"),
               paste0("| ", rows, " |"))
  writeLines(summary, file.path(out, "summary.md"))
  invisible(results)
}
