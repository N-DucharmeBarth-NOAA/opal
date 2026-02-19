# Resolve a biology vector to age-basis

If the vector length matches `n_age`, it is returned unchanged (assumed
to already be on an age basis). If the vector length matches `n_len`
(and `n_len != n_age`), it is converted to age using the
probability-of-length-at-age matrix: `vec_a = t(pla) %*% vec_l`.
Otherwise an informative error is thrown.

## Usage

``` r
resolve_bio_vector(vec, n_age, n_len, pla, name = "vector")
```

## Arguments

- vec:

  Numeric vector (length `n_age` or `n_len`). May be an AD type inside
  [`bet_model()`](https://n-ducharmebarth-noaa.github.io/opal/reference/bet_model.md).

- n_age:

  Integer. Number of age classes.

- n_len:

  Integer. Number of length bins.

- pla:

  Matrix (`n_len` x `n_age`). Probability of length at age (columns sum
  to 1), as returned by
  [`get_pla`](https://n-ducharmebarth-noaa.github.io/opal/reference/get_pla.md).

- name:

  Character. Name of the vector used in warning/error messages.

## Value

Numeric vector of length `n_age` on an age basis.

## Details

**Note:** When `n_age == n_len` the vector is always treated as
age-basis. If you need to pass a length-basis vector in that situation
you must convert it to age externally before calling the model.
