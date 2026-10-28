# Compute mean length-at-age

Computes mean length at each modeled age under one of three growth
curves, selected by `growth_model`: the Schnute parameterization of the
von Bertalanffy (VB) curve, the Richards curve, or the growth cessation
model.

## Usage

``` r
get_growth(
  n_age,
  A1,
  A2,
  L1 = NULL,
  L2 = NULL,
  log_k = NULL,
  min_age = 1L,
  log_k_mult = NULL,
  k_mult_ages = NULL,
  plus_group_adjust = FALSE,
  growth_model = "vonB",
  b = NULL,
  L0 = NULL,
  log_rmax = NULL,
  log_k_gc = NULL,
  t50 = NULL,
  age_offset = 0
)
```

## Arguments

- n_age:

  Integer. Number of age classes.

- A1:

  Integer. Reference age for L1 (data).

- A2:

  Integer. Reference age for L2 (data).

- L1:

  Numeric. Length at age A1 (may be AD). Used by `"vonB"` and
  `"richards"`.

- L2:

  Numeric. Length at age A2 (may be AD). Used by `"vonB"` and
  `"richards"`.

- log_k:

  Numeric. Log growth coefficient per model age step (may be AD). Used
  by `"vonB"` and `"richards"`.

- min_age:

  Integer. Minimum age (default 1L).

- log_k_mult:

  Numeric vector (may be AD). Log K multipliers; `NULL` or length zero
  retains the single-K curve. `"vonB"` only.

- k_mult_ages:

  Integer vector (data). Ages at which the multipliers apply, strictly
  increasing within the modeled ages and matching `log_k_mult`. `"vonB"`
  only.

- plus_group_adjust:

  Logical (data). Apply SS3's oldest-age length approximation,
  independently of the multipliers (default `FALSE`). `"vonB"` only.

- growth_model:

  Character (data). Growth curve: `"vonB"` (default), `"richards"` or
  `"cessation"`. See Details.

- b:

  Numeric. Richards shape parameter on the natural scale (may be AD).
  Must be non-zero. Used by `"richards"`.

- L0:

  Numeric. Length at absolute age 0, natural scale (may be AD). Used by
  `"cessation"`.

- log_rmax:

  Numeric. Log maximum growth rate, in length per model age step (may be
  AD). Used by `"cessation"`.

- log_k_gc:

  Numeric. Log steepness of the logistic decline in growth rate, per
  model age step (may be AD). Used by `"cessation"`.

- t50:

  Numeric. Absolute age at the logistic midpoint (may be AD). Used by
  `"cessation"`.

- age_offset:

  Numeric (data). Added to the modeled ages to give the absolute age `t`
  used by `"cessation"` (default 0).

## Value

Numeric vector of length `n_age`: mean length at each age
`a = min_age, ..., min_age + n_age - 1`.

## Details

### Growth models

|  |  |  |
|----|----|----|
| `growth_model` | Curve | Arguments used |
| `"vonB"` (default) | Schnute VB | `L1`, `L2`, `log_k` (optionally `log_k_mult`, `k_mult_ages`, `plus_group_adjust`) |
| `"richards"` | Richards | `L1`, `L2`, `log_k`, `b` |
| `"cessation"` | Growth cessation | `L0`, `log_rmax`, `log_k_gc`, `t50`, `age_offset` |

Stock Synthesis (SS3) GrowthModel 3 (age-specific K) is reproduced by
`"vonB"` with `log_k_mult` and `k_mult_ages`.

All rates (`log_k`, `log_rmax` and `log_k_gc`) are per model age step.

### von Bertalanffy (growth_model = "vonB")

Uses the Schnute parameterization of the von Bertalanffy growth curve,
\$\$L_a = L_1 + (L_2 - L_1)\frac{1 - e^{-k(a - A_1)}}{1 - e^{-k(A_2 -
A_1)}}.\$\$ Optionally reproduces SS3 GrowthModel 3 (age-specific K),
including its cumulative multipliers and plus-group approximation.

Multipliers apply cumulatively: K at age a is the base K times the
product of all multipliers at ages less than or equal to a. The final K
carries forward. Linf is derived from L1, L2, A1 and A2 using the
oldest-age K; growth from age a to a+1 uses K at a. Consequently, with
multipliers the curve need not pass through L2 at A2. A1 must be a
modeled age for this option. Below A1, opal retains base-K exponential
growth toward this Linf; SS3 instead uses a linear segment.

The plus-group approximation averages a linear increase from the oldest
length to Linf over a further `max_age` ages, with weights
`exp(-0.2)^j`, j = 0, ..., max_age. The fixed 0.2 is not mortality; this
option requires a positive integer oldest age.

`log_k` is the log rate per model age unit, whereas SS3 reads an annual
rate and scales by season duration. For a quarterly age step, use
`log(K_SS3 * 3/12)`; multipliers are dimensionless. This
reimplementation is pinned to SS3 v3.24AB beginning-of-season lengths;
later SS3 versions may differ. CV interpolation is handled separately by
[`get_sd_at_age`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/get_sd_at_age.md).

### Richards (growth_model = "richards")

Uses the Richards (1959) curve as parameterized by Schnute (1981),
\$\$L_a = \left\[L_1^b + (L_2^b - L_1^b)\frac{1 - e^{-k(a - A_1)}}{1 -
e^{-k(A_2 - A_1)}}\right\]^{1/b}.\$\$ The parameterization and the
guidance in this subsection follow the Stock Synthesis user manual
(section 8.5.2.2.2, Richards growth function). See the Richards growth
vignette on the Stock Synthesis website for further technical guidance.

`b = 1` gives the VB curve. As `b` approaches 0, the curve approaches
the Gompertz curve; Gompertz growth can be obtained by fixing `b` at a
small value such as `1e-4`. `b` is a floating-point parameter on the
natural scale that may be positive or negative but not exactly 0. If its
initial value is non-zero, there is effectively no risk of an estimated
`b` landing exactly on 0.

When `A1` is greater than the youngest modeled age, some parameter
combinations give undefined (`NaN`) predicted lengths at the younger
ages. Below `A1` the fraction in the equation is negative, so with
`b > 0` the bracket can become negative. With `b < 0` the asymptotic
bracket can be zero or negative, in which case Linf is undefined. The
choice of `A1` and `A2` also constrains the possible curve shapes. opal
does not clamp the bracket.

### Growth cessation (growth_model = "cessation")

Uses the growth cessation model of Maunder et al. (2018), Eq. 2, \$\$L_t
= L_0 + r\_{max}\left\[\frac{\ln(e^{-k t\_{50}} + 1) - \ln(e^{k(t -
t\_{50})} + 1)}{k} + t\right\],\$\$ with `rmax = exp(log_rmax)` and
`k = exp(log_k_gc)`. The growth rate declines logistically from `rmax`,
with midpoint age `t50`. `L(0) = L0`, and the asymptotic length is
\\L_0 + r\_{max}\[\ln(e^{-k t\_{50}} + 1)/k + t\_{50}\]\\. As `k`
approaches 0, the curve approaches `L0 + rmax * t / 2`. Both log terms
are computed as an overflow-safe \\\log(1 + e^x)\\
([`RTMB::logspace_add`](https://rdrr.io/pkg/RTMB/man/Distributions.html)
for AD inputs), which avoids overflow for large `k * (t - t50)`.

Unlike the VB and Richards curves, which depend only on `a - A1`, this
curve depends on absolute age `t`, so the age origin matters. `t` is the
modeled age plus `age_offset`. For example, if model age 1 corresponds
to biological age 0, use `age_offset = -1`.

## References

Maunder, M.N., Deriso, R.B., Schaefer, K.M., Fuller, D.W.,
Aires-da-Silva, A.M., Minte-Vera, C.V. and Campana, S.E. (2018). The
growth cessation model: a growth model for species showing a near
cessation in growth with application to bigeye tuna (*Thunnus obesus*).
Marine Biology 165: 76.
[doi:10.1007/s00227-018-3336-9](https://doi.org/10.1007/s00227-018-3336-9)
.

Methot, R.D. Jr. and Wetzel, C.R. (2013). Stock synthesis: A biological
and statistical framework for fish stock assessment and fishery
management. Fisheries Research 142: 86–99.
[doi:10.1016/j.fishres.2012.10.012](https://doi.org/10.1016/j.fishres.2012.10.012)
.

Methot, R.D. Jr., Wetzel, C.R., Taylor, I.G., Doering, K.L., Perl, E.F.
and Johnson, K.F. Stock Synthesis User Manual, Version 3.30.25.1. NOAA
Fisheries, Seattle, WA. Section 8.5.2.2.2, Richards growth function.
<https://nmfs-ost.github.io/ss3-doc/SS330_User_Manual_release.html>.

Richards, F.J. (1959). A flexible growth function for empirical use.
Journal of Experimental Botany 10: 290–301.

Schnute, J. (1981). A versatile growth model with statistically stable
parameters. Canadian Journal of Fisheries and Aquatic Sciences 38:
1128–1140.

[Stock Synthesis
source](https://github.com/nmfs-ost/ss3-source-code/blob/d0da9a5/SS3trial.tpl),
sections SS_Label_Info_16.2.2–16.2.4.1.4.

## Examples

``` r
if (FALSE) { # \dontrun{
  # Starting values from the IO yellowfin inputs (quarterly ages, SS3
  # GrowthModel 3 solution). In opal_model() the same options come from
  # data (growth_model, k_mult_ages, growth_plus_group_adjust,
  # growth_age_offset) and parameters (log_L1, log_L2, log_k, log_k_mult,
  # richards_b, log_L0, log_rmax, log_k_gc, t50_gc).
  data(io_yft_data)
  data(io_yft_parameters)
  d <- io_yft_data
  p <- io_yft_parameters
  L1 <- exp(p$log_L1)
  L2 <- exp(p$log_L2)

  # von Bertalanffy (Schnute), the default: single K
  L_vb <- get_growth(d$n_age, d$A1, d$A2, L1, L2, p$log_k)

  # SS3 GrowthModel 3: age-specific K multipliers and plus-group adjustment
  L_ss3 <- get_growth(d$n_age, d$A1, d$A2, L1, L2, p$log_k,
                      log_k_mult = p$log_k_mult,
                      k_mult_ages = d$k_mult_ages,
                      plus_group_adjust = TRUE)

  # Richards (b = 1 recovers VB; b near 0 approximates Gompertz)
  L_rich <- get_growth(d$n_age, d$A1, d$A2, L1, L2, log_k = log(0.25),
                       growth_model = "richards", b = -1)

  # Growth cessation: absolute age t = model age + age_offset
  L_gc <- get_growth(d$n_age, d$A1, d$A2, growth_model = "cessation",
                     L0 = 15, log_rmax = log(7), log_k_gc = log(0.6),
                     t50 = 18, age_offset = 0)

  matplot(cbind(L_vb, L_ss3, L_rich, L_gc), type = "l", lty = 1,
          xlab = "Model age (quarters)", ylab = "Mean length (cm)")
  legend("bottomright", c("vonB", "SS3 GM3", "Richards", "Cessation"),
         col = 1:4, lty = 1)
} # }
```
