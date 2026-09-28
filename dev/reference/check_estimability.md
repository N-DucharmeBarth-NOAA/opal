# Check for identifiability of fixed effects

Calculates the matrix of second-derivatives of the marginal likelihood
with respect to fixed effects, to see if any linear combinations are not
estimable (i.e. cannot be uniquely estimated conditional upon model
structure and available data, e.g., resulting in a likelihood ridge and
singular, non-invertable Hessian matrix)

## Usage

``` r
check_estimability(obj, h)
```

## Arguments

- obj:

  The compiled object

- h:

  optional argument containing pre-computed Hessian matrix

## Value

A tagged list of the hessian and the message
