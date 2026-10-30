# Compute SD of length-at-age

Computes the SD of length at each age from two growth-variability
parameters at the reference lengths `L1` and `L2`, under one of two
forms selected by `sd_model`.

## Usage

``` r
get_sd_at_age(
  mu_a,
  L1,
  L2,
  log_CV1 = NULL,
  log_CV2 = NULL,
  sd_model = "cv",
  log_SD1 = NULL,
  log_SD2 = NULL,
  A1 = NULL,
  A2 = NULL
)
```

## Arguments

- mu_a:

  Numeric vector. Mean length at age (from
  [`get_growth`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/get_growth.md)).
  May be AD.

- L1:

  Numeric. Lower reference length, at age A1 (may be AD).

- L2:

  Numeric. Upper reference length, at age A2 (may be AD).

- log_CV1:

  Numeric. Log CV at `L1` (may be AD). Used by `"cv"`.

- log_CV2:

  Numeric. Log CV at `L2` (may be AD). Used by `"cv"`.

- sd_model:

  Character (data). `"cv"` (default) or `"sd"`. See Details.

- log_SD1:

  Numeric. Log SD at `L1` (may be AD). Used by `"sd"`.

- log_SD2:

  Numeric. Log SD at `L2` (may be AD). Used by `"sd"`.

- A1, A2:

  Integer (data). Optional positions in `mu_a` of the reference ages.
  When both are supplied, the SD at positions below `A1` is held at its
  value at `L1`, and the SD at positions `A2` and above at its value at
  `L2`. When `NULL` (default), the linear form is extrapolated to all
  ages.

## Value

Numeric vector of SD at each age.

## Details

- `"cv"`:

  CV is linear in mean length, \$\$CV_a = CV_1 + \frac{\mu_a -
  L_1}{L_2 - L_1}(CV_2 - CV_1),\quad SD_a = \mu_a CV_a,\$\$ which is
  SS3's `CV_Growth_Pattern = 0` (CV = f(LAA)). SD is then quadratic in
  mean length.

- `"sd"`:

  SD is linear in mean length, \$\$SD_a = SD_1 + \frac{\mu_a -
  L_1}{L_2 - L_1}(SD_2 - SD_1),\$\$ which is SS3's
  `CV_Growth_Pattern = 2` (SD = f(LAA)).

Under both patterns, SS3 holds the CV (or SD) at its young value for
ages below `A1` and at its old value for ages at or above `A2`
(`SS_ALK.tpl`, SS_Label_Info_16.5.3). Supply `A1` and `A2` to do the
same.
[`opal_model`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_model.md)
supplies them under `"sd"` only, so `"cv"` objectives are unchanged.

When SS3's `Growth_Age_for_L2 = 999`, its upper reference length is
Linf, but the SD is still held at the old value from the oldest age.
With `A2` set to the oldest age, `"sd"` then matches SS3 at `A2` and
below `A1`, and differs slightly in between because the slope uses
`mu_a[A2]` instead of Linf.

## References

[Stock Synthesis
source](https://github.com/nmfs-ost/ss3-source-code/blob/main/SS_ALK.tpl),
`SS_ALK.tpl`, SS_Label_Info_16.5.2–16.5.4.
