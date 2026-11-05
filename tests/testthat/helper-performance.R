benchmark_ratio_ci <- function(base, head, replicates = 2000L, seed = 163L) {
  if (length(base) != length(head) || length(base) < 5L ||
      any(!is.finite(c(base, head))) || any(c(base, head) <= 0)) {
    stop("A/B comparison needs at least five paired, positive finite rounds")
  }
  old_seed <- if (exists(".Random.seed", .GlobalEnv)) .GlobalEnv$.Random.seed else NULL
  on.exit({
    if (is.null(old_seed)) {
      if (exists(".Random.seed", .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv)
    } else {
      assign(".Random.seed", old_seed, .GlobalEnv)
    }
  })
  set.seed(seed)
  ratios <- replicate(replicates, {
    idx <- sample.int(length(base), replace = TRUE)
    median(head[idx]) / median(base[idx])
  })
  c(ratio = median(head) / median(base),
    lower = unname(quantile(ratios, 0.025)),
    upper = unname(quantile(ratios, 0.975)))
}

benchmark_instruction_count <- function(path) {
  lines <- readLines(path, warn = FALSE)
  events <- grep("^events:", lines, value = TRUE)
  summary <- grep("^summary:", lines, value = TRUE)
  if (length(events) != 1L || length(summary) != 1L) {
    stop("Missing Cachegrind events/summary in ", path)
  }
  names <- strsplit(trimws(sub("^events:", "", events)), "\\s+")[[1]]
  values <- as.numeric(strsplit(trimws(sub("^summary:", "", summary)), "\\s+")[[1]])
  value <- values[match("Ir", names)]
  if (length(value) != 1L || !is.finite(value) || value <= 0) {
    stop("Invalid Cachegrind instruction count")
  }
  value
}
