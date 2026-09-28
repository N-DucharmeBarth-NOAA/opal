# Linear area rebinning of frequency data

Redistributes counts from source bins to destination bins using
proportional overlap. For each destination bin, the function finds all
source bins that overlap it and assigns a fraction of each source bin's
count proportional to the overlap width relative to the source bin
width.

## Usage

``` r
rebin_counts(src_edges, src_counts, dest_edges)
```

## Arguments

- src_edges:

  numeric vector (length N+1). Boundaries of source bins. Must be
  monotonically increasing.

- src_counts:

  numeric vector (length N). Counts in each source bin. May be an AD
  vector from RTMB.

- dest_edges:

  numeric vector (length M+1). Boundaries of destination bins. Must be
  monotonically increasing.

## Value

numeric vector (length M). Rebinned counts in destination bins.

## Details

AD-safe: src_edges and dest_edges are data; only src_counts may carry AD
values. The operation is linear in src_counts.

## Examples

``` r
# Identity rebin
edges <- 0:5
counts <- c(10, 20, 15, 5, 30)
rebin_counts(edges, counts, edges)
#> [1] 10 20 15  5 30

# Coarsen: merge pairs of bins
src_edges  <- 0:4
src_counts <- c(10, 20, 15, 5)
dest_edges <- c(0, 2, 4)
rebin_counts(src_edges, src_counts, dest_edges)
#> [1] 30 20
```
