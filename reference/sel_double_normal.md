# Double-normal selectivity as a function of length (SS3 pattern 24, full form)

All parameters are on the real line, making this parameterization
suitable for unconstrained optimization and AD. The length-bin vector
`x` is used to define natural scales for position and width parameters.

## Usage

``` r
sel_double_normal(x, par)
```

## Arguments

- x:

  Numeric vector of length-bin midpoints.

- par:

  Numeric vector of length 6 containing selectivity parameters.

  `par[1]` (a)

  :   Peak location (real line). Transformed via `mean(x) + a * sd(x)`,
      so `a = 0` places the peak at `mean(x)`.

  `par[2]` (b)

  :   Plateau width (real line). Controls the distance from peak to the
      start of the descending limb via logistic transform of the
      available range:
      `peak + bin_width + (0.99 * max(x) - peak - bin_width) / (1 + exp(-b))`.

  `par[3]` (c)

  :   Ascending width (real line, log-space). Actual width =
      `exp(c) * sd(x)`. `c = 0` gives an ascending width equal to
      `sd(x)`.

  `par[4]` (d)

  :   Descending width (real line, log-space). Actual width =
      `exp(d) * sd(x)`. `d = 0` gives a descending width equal to
      `sd(x)`.

  `par[5]` (e)

  :   Initial selectivity (real line, logit-space). Transformed via
      `1 / (1 + exp(-e))`, so `e = 0` gives initial selectivity of 0.5,
      large negative values give ~0, large positive values give ~1.

  `par[6]` (f)

  :   Final selectivity (real line, logit-space). Same transform as `e`.

## Value

Numeric vector of selectivity values in \[0, 1\].

## Details

Only the full form is implemented (e and f both active). The simplified
forms where e \<= -999 or f \<= -999 are not supported.
