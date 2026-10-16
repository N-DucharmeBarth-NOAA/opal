# RTMBdist's ordinary density is retained verbatim for fitting and simulation.
# During OSA, factor the composition into sequential beta-binomial conditionals.
# This mirrors RTMB's multinomial OSA factorisation and honours its observation
# permutation. It supports the density indicators used by oneStepGeneric;
# CDF indicators require a different integration algorithm and are rejected.
.opal_ddirmult <- function(x, size, alpha, log = FALSE) {
  if (!inherits(x, "osa")) return(RTMBdist::ddirmult(x, size, alpha, log = log))
  if (!isTRUE(log)) stop("OSA requires log densities.")
  if (ncol(x@keep) != 1L) stop("Dirichlet-multinomial OSA requires method = 'oneStepGeneric'.")
  order <- order(attr(x@keep, "ord"))
  x <- x[order]
  alpha <- alpha[order]
  n <- length(x)
  ans <- 0
  remaining <- size
  if (n > 1L) for (i in seq_len(n - 1L)) {
    # Remaining trial counts are data during the OSA density integration.
    total <- remaining
    count <- x@x[i]
    a <- alpha[i]
    b <- sum(alpha[seq.int(i + 1L, n)])
    density <- lgamma(total + 1) - lgamma(count + 1) - lgamma(total - count + 1) +
      lgamma(count + a) + lgamma(total - count + b) - lgamma(total + a + b) +
      lgamma(a + b) - lgamma(a) - lgamma(b)
    ans <- ans + density * x@keep[i, 1L]
    remaining <- remaining - count
  }
  ans
}
