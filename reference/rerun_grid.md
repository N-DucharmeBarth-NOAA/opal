# Run a grid again

Fit a model for each grid cell defined in `grid` where each grid cell
has a different combination of fixed parameter values.

## Usage

``` r
rerun_grid(grid, bounds, control = list(eval.max = 10000, iter.max = 10000))
```

## Arguments

- grid:

  a `list` of parameter inputs for each grid cell. This can be created
  using `get_grid`.

- bounds:

  the lower and upper parameter bounds created using the `get_bounds`
  function.

- control:

  a `list` defining the controls passed to `nlminb`.

## Value

a `list` of grid cells.
