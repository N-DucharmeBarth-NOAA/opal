# Set up a grid

Set up a grid

## Usage

``` r
get_grid(
  parameters,
  m0 = c(0.5, 0.4, 0.45),
  m10 = c(0.065, 0.085, 0.105),
  h = c(0.55, 0.63, 0.72, 0.8),
  psi = c(1.75, 1.5, 2)
)
```

## Arguments

- parameters:

  a `list` containing the parameter inputs.

- m0:

  the M0 values to be included in the grid.

- m10:

  the M10 values to be included in the grid.

- h:

  the h values to be included in the grid.

- psi:

  the psi values to be included in the grid.

## Value

a `list` of parameter inputs ready to be passed to `MakeADFun`.

## Examples

``` r
if (FALSE) { # \dontrun{
  #parameters <- get_parameters(data = data)
  #grid_parameters <- get_grid(parameters = parameters)
} # }
```
