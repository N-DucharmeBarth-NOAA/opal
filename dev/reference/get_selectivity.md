# Compute selectivity-at-age from parametric selectivity curves

Defines selectivity for each fishery as a parametric curve (logistic,
double-normal, or double Richards) evaluated either on length or
directly on age, chosen per fishery by `data$sel_basis_f`:

- 1 = length (default):

  The curve is evaluated at `len_mid` and converted to
  selectivity-at-age by matrix-multiplying with the
  probability-of-length-at-age (`t(pla) %*% sel_at_length`).

- 2 = age:

  The curve is evaluated directly at the ages `age_a`, and `pla` is not
  used, so the selectivity does not depend on growth. Parameters are
  scaled by `mean(age_a)` and `sd(age_a)` rather than by the length
  bins.

If `data$sel_basis_f` is absent or `NULL`, every fishery is
length-based.

## Usage

``` r
get_selectivity(data, par_sel, pla, len_mid, age_a = NULL)
```

## Arguments

- data:

  A list containing model data. Required elements: n_fishery, n_year,
  n_age, sel_type_f. Optional: `sel_basis_f`, an integer vector of
  length n_fishery (1 = length, 2 = age), and `min_age` (default 1),
  used to build `age_a` when it is not supplied.

- par_sel:

  Numeric matrix of dimensions `[n_fishery, 6]`. Each row is a real-line
  parameter vector. For logistic (sel_type_f == 1), only columns 1:2 are
  used. For double-normal (sel_type_f == 2) and double Richards
  (sel_type_f == 3), all 6 are used.

- pla:

  Numeric matrix of dimensions `[n_len, n_age]` containing
  probability-of-length-at-age (age-length key). Computed via
  [`get_pla()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/get_pla.md).
  Used only for length-based fisheries.

- len_mid:

  Numeric vector (length n_len) of length-bin midpoints.

- age_a:

  Numeric vector (length n_age) of ages at which age-based curves are
  evaluated. Defaults to `seq(min_age, by = 1, length.out = n_age)`, as
  in
  [`get_growth()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/get_growth.md).

## Value

3D array `sel_fya` of dimensions `[n_fishery, n_year, n_age]`.

## Details

Selectivity is currently time-invariant within each fishery (constant
across years).
