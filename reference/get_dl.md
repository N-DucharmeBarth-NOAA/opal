# Obtain dl

Obtain dl for use later in calculating phi. This function is used within
the
[`get_data()`](https://n-ducharmebarth-noaa.github.io/opal/reference/get_data.md)
function so it is generally not needed directly.

## Usage

``` r
get_dl(length_mu_ysa, length_sd_a)
```

## Arguments

- length_mu_ysa:

  an `array` containing the mean length at age (a) for each year (y) and
  season (s).

- length_sd_a:

  a `vector` containing the standard deviation of the mean length at
  age.

## Value

an `array`.

## Details

Obtain dl for use later in calculating phi(a, y). From the
PRELIMINARY_CALCS_SECTION of `opalmod.tpl`. More fine-scale than
lenage_dist_syal. Needs a bit more detail to get it right. Integrates
over length-at-age distribution to get phi(age,year).

## Examples

``` r
if (FALSE) { # \dontrun{
  length_mu_ysa <- get_length_at_age(length_mean = opal::length_mean)
  dl_yal <- get_dl(length_mu_ysa = length_mu_ysa, length_sd_a = opal::length_sd$SD)
} # }
```
