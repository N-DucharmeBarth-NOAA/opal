# Set a prior lambda

Sets the contribution weight for all elements of a prior block or for
one selected element. For an `opal_obj`, the updated object is returned
through
[`opal_update()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_update.md)
so dependent results are invalidated.

## Usage

``` r
opal_set_prior_lambda(x, parameter, value, element = NULL)
```

## Arguments

- x:

  An `opal_obj` or a named priors list.

- parameter:

  Name of the parameter block.

- value:

  Non-negative finite lambda value(s).

- element:

  Optional one-based element index within the block, in column-major
  order for matrices. The block length must be known: from resolved
  parameters on an `opal_obj`, or from a vector `par1`, `par2`, or
  `lambda`.

## Value

The updated `opal_obj` or priors list.
