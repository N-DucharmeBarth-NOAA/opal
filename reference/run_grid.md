# Run a grid

Fit a model for each grid cell defined in `grid` where each grid cell
has a different combination of fixed parameter values.

## Usage

``` r
run_grid(
  data,
  grid_parameters,
  bounds,
  map = list(),
  random = c(),
  control = list(eval.max = 10000, iter.max = 10000)
)
```

## Arguments

- data:

  a `list` containing the data created using the `get_data` function.

- grid_parameters:

  a `list` of parameter inputs for each grid cell. This can be created
  using `get_grid`.

- bounds:

  the lower and upper parameter bounds created using the `get_bounds`
  function.

- map:

  a `list` defining how to optionally collect and fix parameters.

- random:

  a character `vector` defining the random effect parameters.

- control:

  a `list` defining the controls passed to `nlminb`.

## Value

a `list` of grid cells.
