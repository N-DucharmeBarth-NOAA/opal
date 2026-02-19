# Obtain samples from the posterior distribution of reported quantities

When samples from the posterior distribution are obtained using
`tmbstan`, only the model parameters can be accessed directly from the
`stanfit` object (i.e., derived quantities such as `par_B0` or `M_a`
cannot be accessed directly). Accessing derived quantities must be done
using the this function.

## Usage

``` r
get_posterior(
  object,
  posterior,
  pars = "par_B0",
  iters = NULL,
  option = 2,
  type = "df"
)
```

## Arguments

- object:

  A `list` specifying the AD object created using `MakeADFun`.

- posterior:

  An `rstan` objected created using the `tmbstan` function.

- pars:

  The parameter(s) to be extracted.

- iters:

  The number of iterations to be extracted.

- option:

  The parallel option to use.

- type:

  default data frame "df" or "array"

## Value

A `data.frame` or `array`
