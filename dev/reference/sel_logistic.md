# Logistic selectivity as a function of length

Logistic selectivity as a function of length

## Usage

``` r
sel_logistic(len, par)
```

## Arguments

- len:

  Numeric vector of length-bin midpoints.

- par:

  Numeric vector of length 6 containing selectivity parameters. Only the
  first two entries are used: `par[1]` (a) and `par[2]` (b). `a` is the
  inflection point on the real line, transformed via
  `mean(len) + a * sd(len)`, so `a = 0` gives inflection at `mean(len)`.
  `b` is log-scale 95% width, transformed via `exp(b) * sd(len)`, so
  `b = 0` gives a 95% width equal to `sd(len)`.

## Value

Numeric vector of selectivity values in (0, 1\].
