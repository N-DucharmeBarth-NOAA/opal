.opal_composition_input <- function(x, upper, model_edges, tol, default_width = NULL) {
  meta <- c("fishery", "year", "month", "ts")
  out <- list()
  add <- function(info, counts, edges) {
    if (!is.numeric(counts) || any(!is.finite(counts)) || any(counts < 0))
      stop("Composition counts must be finite and nonnegative.", call. = FALSE)
    if (length(edges) != length(counts) + 1L ||
        any(!is.finite(edges)) || any(diff(edges) <= 0))
      stop("Composition edges must be finite and strictly increasing.", call. = FALSE)
    out[[length(out) + 1L]] <<- list(info = info, counts = counts, edges = edges,
                                   input_row = length(out) + 1L)
  }
  frames <- if (is.data.frame(x)) list(x) else x
  if (!is.list(frames) || !all(vapply(frames, is.data.frame, logical(1))))
    stop("Composition input must be a data frame or list of data frames.", call. = FALSE)
  if (!is.null(upper) && length(upper) != 1L && length(upper) != length(frames))
    stop("upper_edge must have length one or one value per input data frame.",
         call. = FALSE)
  native <- is.data.frame(x) && !all(c("bin_lower", "bin_upper", "count") %in% names(x))
  for (j in seq_along(frames)) {
    d <- as.data.frame(frames[[j]])
    if (!all(meta %in% names(d)))
      stop("Composition input requires fishery, year, month, and ts.", call. = FALSE)
    if (all(c("bin_lower", "bin_upper", "count") %in% names(d))) {
      keys <- intersect(c("fishery", "ts", "source"), names(d))
      if (!nrow(d)) next
      group <- do.call(paste, c(lapply(d[keys], as.character), sep = "\r"))
      for (key in unique(group)) {
        rows <- which(group == key)
        rows <- rows[order(d$bin_lower[rows])]
        lo <- d$bin_lower[rows]
        hi <- d$bin_upper[rows]
        if (!is.numeric(lo) || !is.numeric(hi) || any(!is.finite(c(lo, hi))) ||
            any(hi <= lo) || (length(lo) > 1L &&
              any(abs(hi[-length(hi)] - lo[-1L]) > tol)))
          stop("Long composition edges must be finite, increasing and contiguous.", call. = FALSE)
        add(d[rows[1L], intersect(c(meta, "source"), names(d)), drop = FALSE],
            d$count[rows], c(lo, tail(hi, 1L)))
      }
      native <- FALSE
    } else {
      cols <- setdiff(names(d), c(meta, "source"))
      lower <- suppressWarnings(as.numeric(cols))
      if (!length(lower) || any(!is.finite(lower)) || anyDuplicated(lower))
        stop("Wide composition bin names must be unique finite lower edges.", call. = FALSE)
      ord <- order(lower)
      lower <- lower[ord]
      end <- if (is.null(upper)) attr(d, "upper_edge") else
        if (length(upper) == 1L) upper[[1L]] else upper[[j]]
      matches <- length(lower) == length(model_edges) - 1L &&
        all(abs(lower - head(model_edges, -1L)) <= tol)
      if (is.null(end)) {
        end <- if (is.null(default_width)) tail(model_edges, 1L) else
          tail(lower, 1L) + default_width
      }
      if (!is.numeric(end) || length(end) != 1L || !is.finite(end))
        stop("upper_edge must be one finite number per wide data frame.", call. = FALSE)
      edges <- c(lower, end)
      native <- native && matches && abs(end - tail(model_edges, 1L)) <= tol
      counts <- as.matrix(d[, cols[ord], drop = FALSE])
      if ((length(counts) && !is.numeric(counts)) ||
          any(!is.finite(counts)) || any(counts < 0))
        stop("Composition counts must be finite and nonnegative.", call. = FALSE)
      if (any(diff(edges) <= 0))
        stop("Composition edges must be finite and strictly increasing.", call. = FALSE)
      for (i in seq_len(nrow(d)))
        add(d[i, intersect(c(meta, "source"), names(d)), drop = FALSE],
            as.numeric(counts[i, ]), edges)
    }
  }
  list(rows = out, native = native)
}

.opal_pool_composition <- function(x, lo, hi) {
  if (lo > 1L) x[lo] <- sum(x[seq_len(lo)])
  if (hi < length(x)) x[hi] <- sum(x[hi:length(x)])
  x[lo:hi]
}

.opal_composition_grid_width <- function(edges, source_edges, tol) {
  offsets <- abs(c(edges, tail(source_edges, 1L)) - source_edges[1L])
  eps <- max(tol, .Machine$double.eps * max(1, offsets) * 16)
  offsets <- offsets[offsets > eps]
  if (!length(offsets)) return(min(diff(source_edges)))
  grid <- offsets[1L]
  for (value in offsets[-1L]) {
    a <- max(grid, value)
    b <- min(grid, value)
    while (b > eps) {
      remainder <- a %% b
      if (remainder <= eps || b - remainder <= eps) {
        a <- b
        break
      }
      a <- b
      b <- remainder
    }
    grid <- a
  }
  width <- min(diff(source_edges))
  grid / max(1, ceiling(grid / width - eps / width))
}

.opal_composition_registry <- function(rows, source_edges, model_edges,
                                       minbin, maxbin, tails, tol,
                                       prefix, allow_misaligned = FALSE,
                                       native_matrix = NULL) {
  registry <- matrices <- list()
  native_mapping <- logical()
  scheme <- integer(length(rows))
  width <- diff(source_edges)
  for (i in seq_along(rows)) {
    row <- rows[[i]]
    edges <- row$edges
    f <- row$info$fishery
    native <- length(edges) == length(model_edges) &&
      all(abs(edges - model_edges) <= tol)
    if (native) edges <- edges[minbin[f]:(maxbin[f] + 1L)]
    # Snap model-aligned boundaries without changing the recorded observed edges.
    matrix_edges <- edges
    if (prefix == "lf") {
      nearest <- vapply(edges, function(e) source_edges[which.min(abs(source_edges - e))],
                        numeric(1))
      aligned <- abs(edges - nearest) <= tol
      matrix_edges[aligned] <- nearest[aligned]
    }
    same <- which(vapply(registry, function(e)
      length(e) == length(edges) && all(abs(e - edges) <= tol), logical(1)))
    # Identical displayed edges can represent different weight tail operators:
    # native observations retain the old finite intermediate weight grid.
    if (prefix == "wf") same <- same[native_mapping[same] == native]
    if (length(same)) {
      scheme[i] <- same[1L]
      next
    }
    s <- length(registry) + 1L
    if (prefix == "lf") {
      bad <- !aligned
      if (any(bad)) {
        grid <- .opal_composition_grid_width(edges, source_edges, tol)
        split_edges <- edges[bad & edges > source_edges[1L] &
                               edges < tail(source_edges, 1L)]
        share <- if (length(split_edges)) max(vapply(split_edges, function(e) {
          p <- ((e - source_edges[1L]) %% width[1L]) / width[1L]
          min(p, 1 - p)
        }, numeric(1))) else 0
        affected <- which(vapply(rows, function(r)
          length(r$edges) == length(row$edges) &&
            all(abs(r$edges - row$edges) <= tol), logical(1)))
        fisheries <- unique(vapply(rows[affected], function(r) r$info$fishery, numeric(1)))
        input_rows <- vapply(rows[affected], function(r) r$input_row, integer(1))
        msg <- paste0("LF scheme ", s, ", fisheries ", paste(fisheries, collapse = ","),
                      ", rows ", paste(input_rows, collapse = ","), ": misaligned edges ",
                      paste(signif(edges[bad], 10), collapse = ","),
                      "; maximum split share ", signif(share, 5),
                      "; suggested largest grid width dividing offsets: ", signif(grid, 10),
                      ". Check bin labeling/units or use a finer model grid",
                      if (any(edges < source_edges[1L] - tol |
                              edges > tail(source_edges, 1L) + tol))
                        " and extend its range to contain all observed edges" else "",
                      ". Changing mean/sd len_mid selectivity scaling means par_sel is not transferable.")
        if (!allow_misaligned) stop(msg, call. = FALSE)
        warning(msg, call. = FALSE)
      }
      if (any(diff(edges) < min(width) - tol))
        warning("LF scheme ", s, " has observed bins narrower than the model grid.",
                call. = FALSE)
    } else {
      if (tail(edges, 1L) <= source_edges[1L] ||
          edges[1L] >= tail(source_edges, 1L)) {
        warning("WF scheme ", s, " is entirely outside the model weight range.",
                call. = FALSE)
      }
      narrower <- vapply(seq_len(length(edges) - 1L), function(k) {
        overlap <- head(source_edges, -1L) < edges[k + 1L] &
          tail(source_edges, -1L) > edges[k]
        any(diff(edges)[k] < width[overlap] - tol)
      }, logical(1))
      if (any(narrower))
        warning("WF scheme ", s, " has observed bins narrower than overlapping source bins.",
                call. = FALSE)
    }
    if (native && !is.null(native_matrix)) {
      m <- native_matrix
      lo <- minbin[f]
      hi <- maxbin[f]
      if (tails == "accumulate") {
        if (lo > 1L) m[lo, ] <- colSums(m[seq_len(lo), , drop = FALSE])
        if (hi < nrow(m)) m[hi, ] <- colSums(m[hi:nrow(m), , drop = FALSE])
      }
      m <- m[lo:hi, , drop = FALSE]
    } else {
      if (tails == "accumulate") {
        matrix_edges[1L] <- min(matrix_edges[1L], source_edges[1L])
        matrix_edges[length(matrix_edges)] <- max(tail(matrix_edges, 1L),
                                                  tail(source_edges, 1L))
      }
      m <- rebin_matrix(source_edges, matrix_edges)
    }
    registry[[s]] <- edges
    matrices[[s]] <- m
    native_mapping[s] <- native
    scheme[i] <- s
  }
  list(scheme_edges = registry, rebin_s = matrices, scheme = scheme,
       native_i = vapply(rows, function(r) length(r$edges) == length(model_edges) &&
         all(abs(r$edges - model_edges) <= tol), logical(1)),
       nbins_i = as.integer(vapply(rows, function(r) {
         f <- r$info$fishery
         as.integer(if (length(r$edges) == length(model_edges) &&
             all(abs(r$edges - model_edges) <= tol)) maxbin[f] - minbin[f] + 1L
         else length(r$counts))
       }, integer(1))))
}

.opal_prepare_compositions <- function(data, input, prefix, keep, switch, minbin,
                                      maxbin, adjust, cap, add, tails, tol,
                                      allow_misaligned = FALSE) {
  is_lf <- prefix == "lf"
  lower <- seq(if (is_lf) data$len_bin_start else data$wt_bin_start,
               by = if (is_lf) data$len_bin_width else data$wt_bin_width,
               length.out = if (is_lf) data$n_len else data$n_wt)
  width <- if (is_lf) data$len_bin_width else data$wt_bin_width
  model_edges <- c(lower, tail(lower, 1L) + width)
  len_lower <- seq(data$len_bin_start, by = data$len_bin_width, length.out = data$n_len)
  len_edges <- c(len_lower, tail(len_lower, 1L) + data$len_bin_width)
  source_edges <- if (is_lf) len_edges else data$lw_a * len_edges^data$lw_b
  if (any(!is.finite(source_edges)) || any(diff(source_edges) <= 0))
    stop("Model composition edges must be finite and increasing.", call. = FALSE)
  if (is.null(minbin)) minbin <- rep(1L, data$n_fishery)
  if (is.null(maxbin)) maxbin <- rep(as.integer(length(lower)), data$n_fishery)
  if (is.null(adjust)) adjust <- rep(1, data$n_fishery)
  if (length(minbin) != data$n_fishery || length(maxbin) != data$n_fishery ||
      anyNA(c(minbin, maxbin)) || any(minbin < 1 | minbin != floor(minbin)) ||
      any(maxbin > length(lower) | maxbin != floor(maxbin)) || any(maxbin < minbin))
    stop("minbin/maxbin must be valid per-fishery bin indices.", call. = FALSE)
  minbin <- as.integer(minbin)
  maxbin <- as.integer(maxbin)
  if (length(adjust) != data$n_fishery || any(!is.finite(adjust)) || any(adjust <= 0))
    stop("Variance adjustments must be finite and positive.", call. = FALSE)
  if (length(add) != 1L || !is.finite(add) || add < 0)
    stop("addtocomp must be finite and nonnegative.", call. = FALSE)
  if (!is.null(cap) && (length(cap) != 1L || !is.finite(cap) || cap <= 0))
    stop("Composition cap must be positive and finite.", call. = FALSE)
  rows <- input$rows
  fishery <- vapply(rows, function(r) r$info$fishery, numeric(1))
  if (any(!is.finite(fishery)) || any(fishery != floor(fishery)) ||
      any(fishery < 1 | fishery > data$n_fishery))
    stop("Composition fisheries must be valid model indices.", call. = FALSE)
  use <- vapply(rows, function(r) sum(r$counts) > 0, logical(1))
  if (!is.null(keep)) use <- use & fishery %in% keep
  rows <- rows[use]
  fishery <- as.integer(fishery[use])
  n <- vapply(rows, function(r) sum(r$counts), numeric(1))
  raw_prop <- Map(function(r, size) r$counts / size, rows, n)
  if (!is.null(cap)) n <- pmin(n, cap)
  n <- n / adjust[fishery]
  n_int <- as.integer(floor(n + 0.5))
  fisheries <- sort(unique(fishery))
  group <- factor(fishery, levels = fisheries)
  row_fi <- split(seq_along(rows), group)
  flat_order <- unlist(row_fi, use.names = FALSE)
  legacy_matrix <- if (is_lf) NULL else rebin_matrix(source_edges, model_edges)
  registry <- .opal_composition_registry(rows[flat_order], source_edges, model_edges,
    minbin, maxbin, tails, tol, prefix, allow_misaligned, legacy_matrix)
  processed <- lapply(seq_along(rows), function(i) {
    r <- rows[[i]]
    p <- raw_prop[[i]]
    native <- length(r$edges) == length(model_edges) &&
      all(abs(r$edges - model_edges) <= tol)
    if (native) p <- .opal_pool_composition(p, minbin[fishery[i]], maxbin[fishery[i]])
    p <- p + add
    p / sum(p) * n[i]
  })
  obs <- lapply(row_fi, function(indices) {
    z <- processed[indices]
    if (length(unique(lengths(z))) == 1L) {
      m <- do.call(rbind, z)
      dimnames(m) <- NULL
      m
    } else z
  })
  names(obs) <- NULL
  set <- function(name, value) data[[paste0(prefix, "_", name)]] <<- value
  set("switch", as.integer(switch))
  raw <- if (!length(raw_prop)) matrix(numeric(), 0L, length(lower)) else
    if (length(unique(lengths(raw_prop))) == 1L) do.call(rbind, raw_prop) else raw_prop
  if (is.matrix(raw)) {
    dimnames(raw) <- NULL
    if (input$native) colnames(raw) <- as.character(lower)
  }
  set("obs_in", raw)
  set("n", n)
  set("n_int", n_int)
  set("fishery", fishery)
  set("fishery_f", fisheries)
  set("n_f", as.integer(table(group)))
  year <- as.integer(vapply(rows, function(r) r$info$ts, numeric(1)))
  set("year", year)
  set("season", rep(1L, length(rows)))
  set("n_fi", split(n, group))
  set("n_int_fi", split(n_int, group))
  set("year_fi", split(year, group))
  set("row_fi", row_fi)
  set("minbin", minbin)
  set("maxbin", maxbin)
  set("var_adjust", adjust)
  set("addtocomp", add)
  set("tails", tails)
  set("edge_tol", tol)
  for (field in names(registry)) set(field, registry[[field]])
  set("obs_data", list(obs = obs))
  set("obs_flat", as.numeric(unlist(processed[flat_order], use.names = FALSE)))
  set("obs_ints", as.integer(unlist(lapply(flat_order, function(i)
    .opal_integerize_composition(processed[[i]], n_int[i])$counts), use.names = FALSE)))
  set("obs_prop", as.numeric(unlist(lapply(processed[flat_order], function(x) {
    p <- x / sum(x)
    p <- p + add
    p / sum(p)
  }), use.names = FALSE)))
  set("nbins", if (length(registry$nbins_i)) registry$nbins_i[1L] else 0L)
  zero <- n_int == 0L
  set("dm_zero_n_f", tabulate(fishery[zero], data$n_fishery))
  set("dm_zero_weight_f", vapply(seq_len(data$n_fishery), function(f)
    sum(n[zero & fishery == f]), numeric(1)))
  if (any(zero)) message(toupper(prefix), " Dirichlet-multinomial integerization skipped ",
                        sum(zero), " observations with total effective sample size ",
                        signif(sum(n[zero]), 6), ".")
  data[[paste0("n_", prefix)]] <- as.integer(length(rows))
  if (is_lf) {
    data$len_lower <- lower
    data$len_upper <- lower + width
    data$len_mid <- lower + width / 2
    data$removal_switch_f <- rep(0L, data$n_fishery)
    data$lf_obs <- data$lf_obs_vec <- data$lf_nbins_f <- NULL
  } else {
    data$wt_lower <- lower
    data$wt_upper <- lower + width
    data$wt_mid <- lower + width / 2
    data$wt_bin_edges <- model_edges
    data$wf_rebin_matrix <- legacy_matrix
  }
  data
}
