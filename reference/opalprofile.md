# Adaptive likelihood profiling

Calculate 1D likelihood profiles with respect to single parameters or
more generally, with respect to arbitrary linear combinations of
parameters (e.g. contrasts).

## Usage

``` r
opalprofile(
  obj,
  name,
  lincomb,
  h = 1e-04,
  ytol = 2,
  ystep = 0.1,
  maxit = ceiling(5 * ytol/ystep),
  parm.range = c(-Inf, Inf),
  slice = FALSE,
  adaptive = TRUE,
  trace = TRUE
)
```

## Arguments

- obj:

  Object from MakeADFun that has been optimized.

- name:

  Name or index of a parameter to profile.

- lincomb:

  Optional linear combination of parameters to profile. By default a
  unit vector corresponding to name.

- h:

  Initial adaptive stepsize on parameter axis.

- ytol:

  Adjusts the range of the likelihood values.

- ystep:

  Adjusts the resolution of the likelihood profile.

- maxit:

  Max number of iterations for adaptive algorithm.

- parm.range:

  Valid parameter range.

- slice:

  Do slicing rather than profiling?

- adaptive:

  Logical; Use adaptive step size?

- trace:

  Trace progress? (TRUE, or a numeric value of 1, gives basic tracing:
  numeric values \> 1 give more information).

## Value

a `vector` penalty
