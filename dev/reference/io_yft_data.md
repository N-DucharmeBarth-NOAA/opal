# Indian Ocean yellowfin tuna (simulated) assessment data

Model inputs for a quarterly, single-region assessment of simulated
Indian Ocean yellowfin tuna, converted from the IATTC analyst team's
Stock Synthesis 3 (SS3) model for the Spatial Stock Assessment
Simulation Experiment (replicate 4). Used by the "Comparing opal with
Stock Synthesis" article.

## Usage

``` r
io_yft_data
```

## Format

A named list of 62 elements in the form
[`opal_obj()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/opal_obj.md)
expects: dimensions, catch (`catch_obs_ysf`), the index (`cpue_data`),
biology, growth settings, length compositions (`lf_*`), the recruitment
bias adjustment (`bias_adj_y`) and priors.

## Source

Converted from SS3 v3.24 files with the SS3-to-opal converter in
<https://github.com/N-DucharmeBarth/wcpfc-bet-single-region>
(`code/rtmb/scripts/04-export-io-yft-for-opal.r`). SS3 model by M.N.
Maunder, H. Xu and C. Minte-Vera (IATTC); data simulated from the
operating model described in Goethel et al. (2024)
[doi:10.1111/faf.12819](https://doi.org/10.1111/faf.12819) .

## Details

The model runs 256 quarterly time steps (`n_year = 256`, `n_season = 1`)
and SS3 ages 1–28, so the opal age index equals the SS3 age. There are
20 fisheries: 16 extraction fleets with catch in numbers, and the
longline CPUE index split into four quarterly series (fisheries 17–20),
each with its own catchability. Length compositions (5 cm bins) are
multinomial for fisheries 1–14 and 16, with SS3's `mult_by_lencomp_N`
carried as `lf_var_adjust`;
[`prep_lf_data()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/prep_lf_data.md)
has already been applied.

Biology is SS3's age-based end-of-run biology: M, weight, maturity and
spawning potential at age, with M scaled to the quarter. Growth uses
SS3's `GrowthModel 3` (age-specific K): `k_mult_ages` gives the ages
with K multipliers, which are in `io_yft_parameters$log_k_mult`, and
`growth_plus_group_adjust` applies SS3's plus-group length
approximation. Selectivity is length-based (`sel_type_f`: 2 = double
normal for fisheries 1–16, 1 = logistic for the index series).
`sel_basis_f` is not set, so every curve is evaluated on length; adding
`sel_basis_f = rep(2L, 20)` evaluates them on age instead, as SS3 does
(see the "Length- or age-based selectivity" article). `priors` holds the
default priors from
[`get_priors()`](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/get_priors.md).

## See also

[io_yft_parameters](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/io_yft_parameters.md),
[io_yft_map](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/io_yft_map.md),
[io_yft_reference](https://n-ducharmebarth-noaa.github.io/opal/dev/reference/io_yft_reference.md)

## Examples

``` r
data(io_yft_data)
io_yft_data$n_year
#> [1] 256
head(io_yft_data$cpue_data)
#>   year month  ts fishery index metric units    value        se
#> 1   81     2  81      17     1   cpue     2 1.399491 0.0628665
#> 2   85     2  85      17     1   cpue     2 1.244067 0.1169320
#> 3   89     2  89      17     1   cpue     2 1.235876 0.0581591
#> 4   93     2  93      17     1   cpue     2 1.201179 0.0749399
#> 5   97     2  97      17     1   cpue     2 1.194558 0.0583808
#> 6  101     2 101      17     1   cpue     2 1.261683 0.1402180
```
