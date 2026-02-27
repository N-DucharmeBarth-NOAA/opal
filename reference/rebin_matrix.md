# Compute rebinning weight matrix

Returns an M x N matrix W where `W[i,j]` is the fraction of source bin j
that falls within destination bin i. Applying `W %*% src_counts` gives
the same result as `rebin_counts(src_edges, src_counts, dest_edges)`.

## Usage

``` r
rebin_matrix(src_edges, dest_edges)
```

## Arguments

- src_edges:

  numeric vector (length N+1). Source bin boundaries. Must be
  monotonically increasing.

- dest_edges:

  numeric vector (length M+1). Destination bin boundaries. Must be
  monotonically increasing.

## Value

numeric matrix (M x N) of overlap fractions.

## Details

Precomputing W is efficient when the same bin edges are reused for many
observations (e.g., all weight composition observations share the same
length-to-weight bin mapping).

## Examples

``` r
src_edges  <- 0:4
dest_edges <- c(0, 2, 4)
W <- rebin_matrix(src_edges, dest_edges)
# W %*% c(10, 20, 15, 5) gives the coarsened counts
```
