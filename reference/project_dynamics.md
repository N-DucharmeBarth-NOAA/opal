# Project dynamics

Project dynamics

## Usage

``` r
project_dynamics(
  data,
  object,
  mcmc = NULL,
  n_proj = 5,
  n_iter = 1,
  rdev_y,
  sel_fya,
  catch_ysf
)
```

## Arguments

- data:

  a `list` of parameter values.

- object:

  a `list` of parameter values.

- mcmc:

  a `list` of parameter values.

- n_proj:

  the number of projection years.

- n_iter:

  a `list` of inputs.

- rdev_y:

  a `list` of inputs.

- sel_fya:

  a `list` of inputs.

- catch_ysf:

  a `list` of inputs.

## Value

the negative log-likelihood (NLL) value.
