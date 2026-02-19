# Obtain the length at age array

Obtain the mean length (cm) at age for each year and season. This
function is used within `get_data` so is generally not needed directly.
It simply structures the input `data.frame` into an `array`.

## Usage

``` r
get_length_at_age(length_mean)
```

## Arguments

- length_mean:

  a `data.frame` containing the mean length at age (a) for each year (y)
  and season (s).

## Value

an `array`.

## Details

The input `data.frame` must contain the mean length (cm) at age arranged
into the columns `Year`, `Season`, and ages 0 to 30:

|      |        |     |      |     |       |       |
|------|--------|-----|------|-----|-------|-------|
| Year | Season | 0   | 1    | ... | 29    | 30    |
| 1931 | 1      | 45  | 57.1 | ... | 74.1  | 88.8  |
| 1932 | 1      | 45  | 57.1 | ... | 74.1  | 88.8  |
| ...  | ...    | ... | ...  | ... | ...   | ...   |
| 2021 | 2      | 50  | 68.9 | ... | 184.0 | 184.0 |
| 2022 | 2      | 50  | 68.9 | ... | 184.0 | 184.0 |

## Examples

``` r
length_mu_ysa <- get_length_at_age(length_mean = opal::length_mean)
#> Error: 'length_mean' is not an exported object from 'namespace:opal'
```
