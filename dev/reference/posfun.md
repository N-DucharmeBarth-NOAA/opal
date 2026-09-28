# Positive Constraint Penalty Function

Ensures a value remains above a small threshold using a smooth
approximation and penalty. See https://github.com/kaskr/adcomp/issues/7
for discussion.

## Usage

``` r
posfun(x, eps = 0.001)
```

## Arguments

- x:

  Numeric value to constrain.

- eps:

  Minimum allowable value (default 0.001).

## Value

A list with:

- new:

  Transformed value.

- penalty:

  Penalty applied if `x < eps`.
