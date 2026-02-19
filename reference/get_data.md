# Set up the data input file

Set up the data input file to be passed to `MakeADFun`. This function
runs data cross validation tests and appends several inputs to the data
list including model dimensions and processed inputs:

## Usage

``` r
get_data(data_in)
```

## Arguments

- data_in:

  a `list` containing the data inputs.

## Value

a `list` ready to be passed to `MakeADFun`.

## Details

- `n_year`: dervied from `first_yr` and `last_yr`

- `n_season`: set to 2

- `n_length`: not in use

- `n_age`: derived from `min_age` and `max_age`

- `n_fishery`: set to 6

- `age_a`: sequence of modeled ages derived from `min_age` and `max_age`

- `length_mu_ysa`: derived from the `length_mean` input

- `length_sd_a`: derived from the `length_sd` input

- `dl_yal`: derived from `length_mu_ysa` and `length_sd_a`

- `weight_fya`: derived from `length_mu_ysa` and `length_sd_a`

- `catch_obs_ysf`: derived from `catch`, `catch_UA`, `scenarios_LL1`,
  and `scenarios_surf`

- `sel_change_year_fy`: derived from `sel_change_sd_fy`

This function produces the data input file to be passed to `MakeADFun`.
