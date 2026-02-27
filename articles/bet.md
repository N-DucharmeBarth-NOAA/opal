# The BET model

## Introduction

`opal` is an open, modular R package for fisheries stock assessment
built on [RTMB](https://github.com/kaskr/RTMB) for automatic
differentiation. It includes bundled data from the western and central
Pacific Ocean (WCPO) bigeye tuna (BET) stock assessment. This vignette
demonstrates model setup, fitting, and diagnostics using those data.

## Load inputs

Load the `opal` package and the `RTMB` dependency. The `ggplot2` package
is used for plotting.

``` r
# library(ggplot2)
# library(dplyr)
library(tidyverse)
library(opal)
library(RTMB)

theme_set(theme_bw())
```

The bundled data object `wcpo_bet_data` contains all biological
parameters, catch and CPUE observations, length structure, and prior
specifications needed for the BET assessment model:

``` r
data(wcpo_bet_data)
data <- wcpo_bet_data

names(data)
#>  [1] "age_a"          "n_age"          "n_season"       "n_fishery"     
#>  [5] "len_bin_start"  "len_bin_width"  "n_len"          "first_yr"      
#>  [9] "last_yr"        "years"          "n_year"         "first_yr_catch"
#> [13] "catch_units_f"  "cpue_switch"    "cpue_data"      "A1"            
#> [17] "A2"             "lw_a"           "lw_b"           "maturity"      
#> [21] "fecundity"      "M"              "catch_obs_ysf"  "sel_type_f"    
#> [25] "priors"
```

Key dimensions of the data:

``` r
cat("Number of ages:", data$n_age, "\n")
#> Number of ages: 40
cat("Number of years:", data$n_year, "\n")
#> Number of years: 268
cat("Number of fisheries:", data$n_fishery, "\n")
#> Number of fisheries: 15
cat("Number of length bins:", data$n_len, "\n")
#> Number of length bins: 95
```

### Biological inputs

Natural mortality at age, maturity at length, and length-weight
parameters are all contained in the data object.

``` r
ages <- data$age_a
real_age <- ages / 4

# Natural mortality
ggplot(data.frame(age = real_age, M = data$M), aes(x = age, y = M)) +
  geom_line() + geom_point() +
  labs(x = "Age (years)", y = "M (quarterly rate)",
       title = "Natural Mortality-at-Age")
```

![](bet_files/figure-html/plot-biology-1.png)

``` r

# Maturity at length
len_mid <- seq(data$len_bin_start + data$len_bin_width / 2,
               by = data$len_bin_width, length.out = data$n_len)
ggplot(data.frame(length = len_mid, maturity = data$maturity),
       aes(x = length, y = maturity)) +
  geom_line() + geom_point() +
  labs(x = "Length (cm)", y = "Maturity",
       title = "Maturity-at-Length")
```

![](bet_files/figure-html/plot-biology-2.png)

### Length composition data

Length composition observations are loaded from the bundled
`wcpo_bet_lf` long-format data object and transformed into model-ready
arrays.

``` r
data(wcpo_bet_lf)

# Pivot to wide format: one row per fishery x timestep, bins as columns
# Note: LF data exist for fisheries 8-14, but we initially use only fisheries 8 & 9
lf_wide <- wcpo_bet_lf %>%
  tidyr::pivot_wider(
    id_cols = c(fishery, year, month, ts),
    names_from = bin,
    values_from = value,
    values_fill = 0
  ) %>%
  arrange(fishery, ts)

# define variance adjustment scalars for each fishery if using lf_switch = 1 (multinomial)
var_adjust_scalars <- 1/rep(1,data$n_fishery)
# var_adjust_scalars <- 1/rep(20000,data$n_fishery)
# var_adjust_scalars[c(1,4,5,6,15)] <- 1/40000

data <- prep_lf_data(data, lf_wide, lf_keep_fisheries = c(8, 9),
                     lf_var_adjust = var_adjust_scalars)
```

### Weight composition data

Weight composition observations are loaded from the bundled
`wcpo_bet_wf` long-format data object and transformed into model-ready
arrays.

``` r
data(wcpo_bet_wf)

# Weight bin scalars (1 kg bins, 1–200 kg)
data$wt_bin_start <- 1
data$wt_bin_width <- 1
data$n_wt         <- 200L

# Pivot to wide format: one row per fishery x timestep, bins as columns
wf_wide <- wcpo_bet_wf |>
  tidyr::pivot_wider(
    id_cols     = c(fishery, year, month, ts),
    names_from  = bin,
    values_from = value,
    values_fill = 0
  ) |>
  dplyr::arrange(fishery, ts)

data <- prep_wf_data(data, wf_wide,
                     wf_keep_fisheries = c(2),
                     wf_switch = 1L)
```

## Model setup

### Parameters

Define the initial parameter values. Growth parameters (`log_L1`,
`log_L2`, `log_k`, `log_CV1`, `log_CV2`) and selectivity (`par_sel`) are
initialised at reasonable starting values. Load these from the bundled
`wcpo_bet_parameters` data object:

``` r
data(wcpo_bet_parameters)

parameters <- list(
  log_B0 = 20,
  log_h = as.numeric(wcpo_bet_parameters$log_h),
  log_sigma_r = as.numeric(wcpo_bet_parameters$log_sigma_r),
  log_cpue_q = as.numeric(wcpo_bet_parameters$log_cpue_q),
  cpue_creep = as.numeric(wcpo_bet_parameters$cpue_creep),
  log_cpue_tau = as.numeric(wcpo_bet_parameters$log_cpue_tau),
  log_cpue_omega = as.numeric(wcpo_bet_parameters$log_cpue_omega),
  log_lf_tau = as.numeric(log(rep(0.1, data$n_fishery))),
  log_wf_tau = rep(0, data$n_fishery),
  log_L1 = as.numeric(wcpo_bet_parameters$log_L1),
  log_L2 = as.numeric(wcpo_bet_parameters$log_L2),
  log_k = as.numeric(wcpo_bet_parameters$log_k),
  log_CV1 = as.numeric(wcpo_bet_parameters$log_CV1),
  log_CV2 = as.numeric(wcpo_bet_parameters$log_CV2),
  par_sel = as.matrix(wcpo_bet_parameters$par_sel),
  rdev_y = as.numeric(wcpo_bet_parameters$rdev_y)
)
```

### Priors

Priors are specified using
[`get_priors()`](https://n-ducharmebarth-noaa.github.io/opal/reference/get_priors.md).
The data object already contains prior center values for growth
parameters:

``` r
data$priors <- get_priors(parameters = parameters, data = data)
evaluate_priors(parameters = parameters, priors = data$priors)
#> [1] 92427.24
```

### Parameter map

Use RTMB’s `map` option to turn parameters on/off. Parameters mapped to
`factor(NA)` are fixed at their initial values:

``` r
map_sel <- matrix(NA, nrow(parameters$par_sel), ncol(parameters$par_sel))
map_sel[8, 1] <- NA
map_sel[9, 1] <- NA
map_sel[8, 3] <- NA
map_sel[9, 3] <- NA
map_sel[8, 4] <- NA
map_sel[9, 4] <- NA
map_rdev <- rep(NA, length(parameters$rdev_y))

map_lf_tau <- rep(NA, length(parameters$log_lf_tau))
map_lf_tau[c(8, 9)] <- NA

map <- list(
  # log_B0 = factor(NA),
  log_h = factor(NA),
  log_sigma_r = factor(NA),
  # log_cpue_q = factor(NA),
  cpue_creep = factor(NA),
  log_cpue_tau = factor(NA),
  log_cpue_omega = factor(NA),
  log_lf_tau = factor(map_lf_tau),
  log_wf_tau = factor(rep(NA, data$n_fishery)),  # fixed initially
  log_L1  = factor(NA),
  log_L2  = factor(NA),
  log_k   = factor(NA),
  log_CV1 = factor(NA),
  log_CV2 = factor(NA),
  par_sel = factor(map_sel)
  # rdev_y = as.factor(map_rdev),
)
```

### Build the AD object

Using the `data`, the `parameters`, the parameter `map`, and the model
(`opal_model`), the AD object is created using RTMB’s `MakeADFun`
function:

``` r
# data$lf_switch <- 1 # multinomial likelihood on flat counts (default)
# data$lf_switch <- 2 # fails while fitting, need to sort out log_lf_tau for this, use simulate to tune and find
# data$lf_switch <- 3 # fails - wants integers - think this should be an issue to RTMBdist guys
data$lf_switch <- 0 # skip length comps (removal only)

# Note: wf_switch was already set by prep_wf_data(); it is re-stated here for
# clarity alongside the equivalent lf_switch assignment.
# data$wf_switch <- 2 # Dirichlet
# data$wf_switch <- 3 # Dirichlet-multinomial
data$wf_switch <- 0 # skip weight comps (removal only)

# Optionally, we could also specify random effects here (e.g. `random = "rdev_y"`), but we'll start with a simpler fixed-effects model to check everything is working first.
obj <- MakeADFun(func = cmb(opal_model, data), parameters = parameters, map = map)
unique(names(obj$par))
#> [1] "log_B0"     "log_cpue_q" "rdev_y"
obj$fn()
#> [1] 696.6957
obj$gr()
#> outer mgc:  847.8735
#>         [,1]     [,2]      [,3]      [,4]      [,5]      [,6]      [,7]
#> [1,] 1.36223 847.8735 -14.42373 -17.79091 -11.96625 -12.08163 -23.62672
#>           [,8]      [,9]   [,10]     [,11]     [,12]     [,13]     [,14]
#> [1,] -17.25554 -17.79412 -22.425 -21.47014 -9.492821 -9.487044 -34.25533
#>          [,15]     [,16]     [,17]     [,18]     [,19]     [,20]     [,21]
#> [1,] -13.94379 -24.41818 -21.90769 -13.48098 -7.581473 -20.25592 -26.55708
#>          [,22]     [,23]     [,24]     [,25]     [,26]     [,27]     [,28]
#> [1,] -11.59815 -17.71142 -14.77482 -11.33643 -8.327975 -8.960074 -11.67251
#>          [,29]     [,30]     [,31]     [,32]     [,33]     [,34]     [,35]
#> [1,] -14.40482 -18.79487 -10.89254 -12.52937 -13.77129 -16.75682 -7.612175
#>          [,36]     [,37]     [,38]     [,39]     [,40]     [,41]     [,42]
#> [1,] -19.80659 -13.20922 -10.29418 -14.54024 -12.02668 -28.96289 -10.74707
#>          [,43]     [,44]     [,45]     [,46]     [,47]     [,48]     [,49]
#> [1,] -10.28152 -14.14101 -22.33511 -9.906668 -8.281704 -11.64379 -9.064911
#>         [,50]     [,51]     [,52]     [,53]     [,54]     [,55]     [,56]
#> [1,] -10.2941 -9.664533 -13.90955 -8.889439 -10.53537 -21.16136 -6.615962
#>         [,57]     [,58]    [,59]     [,60]     [,61]     [,62]     [,63]
#> [1,] -7.25638 -8.705425 -6.17196 -5.831573 -6.230935 -4.644877 -12.92899
#>          [,64]     [,65]     [,66]     [,67]     [,68]     [,69]     [,70]
#> [1,] -4.819591 -4.237939 -4.831417 -4.725637 -3.043372 -2.803211 -2.736426
#>          [,71]     [,72]     [,73]     [,74]     [,75]     [,76]     [,77]
#> [1,] -3.042287 -2.957644 -2.832817 -2.951742 -3.327376 -3.373636 -2.718246
#>          [,78]     [,79]     [,80]     [,81]      [,82]     [,83]     [,84]
#> [1,] -1.905599 -1.534254 -1.112567 -2.980091 -0.7524688 -3.513066 -1.777633
#>         [,85]     [,86]     [,87]     [,88]     [,89]     [,90]     [,91]
#> [1,] -2.77622 -2.898418 -3.203531 -3.533987 -3.949122 -3.965853 -3.891109
#>          [,92]     [,93]     [,94]     [,95]     [,96]     [,97]     [,98]
#> [1,] -3.869417 -3.657707 -3.873963 -3.423115 -2.867199 -2.524664 -2.428816
#>          [,99]    [,100]    [,101]   [,102]  [,103]    [,104]   [,105]   [,106]
#> [1,] -1.881653 0.1099419 0.8935999 4.313273 3.53943 -1.380832 3.849519 1.001104
#>        [,107]    [,108]   [,109]   [,110]   [,111]   [,112]     [,113]   [,114]
#> [1,] 1.260817 -1.860571 2.235185 1.621885 3.384444 3.651069 -0.4034103 10.06855
#>        [,115]   [,116]   [,117]   [,118]   [,119]   [,120]   [,121]   [,122]
#> [1,] 4.423438 7.950092 9.743676 2.905259 4.412664 10.42574 6.050362 5.262348
#>        [,123]   [,124]    [,125]    [,126]   [,127]   [,128]   [,129]   [,130]
#> [1,] 19.63047 1.499814 0.0709136 0.8013582 10.39791 13.36318 4.934542 17.81396
#>          [,131]   [,132]   [,133]   [,134]   [,135]   [,136]   [,137]   [,138]
#> [1,] -0.6327937 1.060031 2.913669 3.210989 2.645825 1.973106 8.814456 1.905641
#>          [,139]    [,140]   [,141]   [,142]   [,143]   [,144]  [,145]  [,146]
#> [1,] -0.5361419 -1.843739 16.19188 3.008929 2.066632 5.887302 1.33293 2.34599
#>        [,147]   [,148]   [,149]   [,150]    [,151]   [,152]   [,153]   [,154]
#> [1,] 2.694726 6.489514 1.825036 1.483763 0.3097918 7.556589 3.024082 9.470925
#>        [,155]   [,156]   [,157]   [,158]   [,159]  [,160]   [,161]   [,162]
#> [1,] 4.270987 5.030244 3.461459 5.489262 7.919619 4.12505 4.281111 1.847724
#>        [,163]   [,164]   [,165]    [,166]   [,167]   [,168]    [,169]   [,170]
#> [1,] 4.132976 2.679274 2.881558 0.5435661 2.279492 8.243004 -1.876656 8.524661
#>        [,171]   [,172]  [,173]   [,174]   [,175]   [,176]  [,177]   [,178]
#> [1,] 3.978502 1.049552 5.94895 4.913562 4.518005 4.382705 5.11398 4.263811
#>        [,179]   [,180]    [,181]   [,182]   [,183]   [,184]   [,185]   [,186]
#> [1,] 3.804904 1.952369 0.5251225 21.01821 7.529496 11.68518 10.81984 2.003909
#>        [,187]   [,188]   [,189]   [,190]   [,191]   [,192]  [,193]     [,194]
#> [1,] 25.51814 13.05216 6.110153 5.844232 12.74616 1.187579 9.65867 -0.8051976
#>        [,195]   [,196]   [,197]   [,198]   [,199]   [,200]   [,201]   [,202]
#> [1,] 5.905875 3.543131 3.104756 7.223417 4.720595 12.58776 12.51461 20.60459
#>        [,203]   [,204]   [,205]   [,206]   [,207]   [,208]    [,209]   [,210]
#> [1,] 15.85585 13.19304 14.98456 5.014382 7.704384 4.966446 0.9461077 8.565709
#>        [,211]   [,212]   [,213]   [,214]  [,215]   [,216]   [,217]   [,218]
#> [1,] 28.13899 31.12039 30.31993 11.78485 24.7352 6.036837 17.65149 22.38168
#>        [,219]   [,220]   [,221]   [,222]   [,223]  [,224]  [,225]   [,226]
#> [1,] 14.94333 14.78684 4.432771 16.58045 2.901582 11.1743 9.21042 11.13793
#>        [,227]   [,228]   [,229]   [,230]   [,231]   [,232]   [,233]     [,234]
#> [1,] 14.87578 14.09604 13.00995 18.83471 6.427051 13.80929 9.082483 -0.2643707
#>        [,235]   [,236]   [,237]   [,238]   [,239]   [,240]  [,241]   [,242]
#> [1,] 10.46214 31.25587 8.392217 9.580797 1.915908 6.143048 8.12524 4.091881
#>        [,243]   [,244]   [,245]   [,246]   [,247]  [,248]   [,249]   [,250]
#> [1,] 20.24496 9.960366 3.628648 2.183022 5.307158 9.07459 5.761786 5.801557
#>        [,251]       [,252]   [,253]   [,254]   [,255]   [,256]   [,257] [,258]
#> [1,] 6.351286 -0.001576827 8.981838 6.469136 5.821842 2.562364 12.84008 3.0483
#>        [,259]   [,260]   [,261]   [,262]    [,263]   [,264]   [,265]    [,266]
#> [1,] 5.601956 4.766498 2.654791 5.269718 0.8332309 1.495876 2.504979 -1.093281
#>          [,267]    [,268]    [,269]    [,270]
#> [1,] -0.7009227 -2.285167 -2.101524 0.6256006
```

Inspect initial model outputs:

``` r
plot(obj$report()$spawning_biomass_y, type = "l",
     xlab = "Time step", ylab = "Spawning biomass (mt)",
     main = "Initial spawning biomass trajectory")
```

![](bet_files/figure-html/init-checks-1.png)

``` r

plot_catch(data = data, obj = obj)
#> [1] "The maximum catch difference was: 1.10394466901198e-08"
```

![](bet_files/figure-html/init-checks-2.png)

### Parameter bounds

``` r
Lwr <- rep(-Inf, length(obj$par))
Upr <- rep(Inf, length(obj$par))
Lwr[grep("log_B0", names(obj$par))] <- log(1)
Upr[grep("log_B0", names(obj$par))] <- log(exp(22))
Lwr[grep("log_cpue_q", names(obj$par))] <- log(0.1)
Upr[grep("log_cpue_q", names(obj$par))] <- log(10)
Lwr[grep("log_lf_tau", names(obj$par))] <- rep(-9, length(grep("log_lf_tau", names(obj$par))))
Upr[grep("log_lf_tau", names(obj$par))] <- rep(9, length(grep("log_lf_tau", names(obj$par))))
Lwr[grep("rdev_y", names(obj$par))] <- rep(-5, length(grep("rdev_y", names(obj$par))))
Upr[grep("rdev_y", names(obj$par))] <- rep(5, length(grep("rdev_y", names(obj$par))))
# Lwr[grep("par_sel", names(obj$par))] <- rep(-7, length(grep("par_sel", names(obj$par))))
# Upr[grep("par_sel", names(obj$par))] <- rep(7, length(grep("par_sel", names(obj$par))))
bounds <- data.frame(par = names(obj$par), lower = Lwr, upper = Upr)
```

## Optimisation

Optimise using the `nlminb` function, do it twice to be sure to be sure
(said with Irish accent):

``` r
control <- list(eval.max = 10000, iter.max = 10000)
opt <- nlminb(start = obj$par, objective = obj$fn, gradient = obj$gr,
              hessian = obj$he, lower = Lwr, upper = Upr, control = control)
#> Warning in nlminb(start = obj$par, objective = obj$fn, gradient = obj$gr, :
#> NA/NaN function evaluation
#> Warning in nlminb(start = obj$par, objective = obj$fn, gradient = obj$gr, :
#> NA/NaN function evaluation
#> Warning in nlminb(start = obj$par, objective = obj$fn, gradient = obj$gr, :
#> NA/NaN function evaluation
opt <- nlminb(start = opt$par, objective = obj$fn, gradient = obj$gr,
              hessian = obj$he, lower = Lwr, upper = Upr, control = control)
obj$fn()
max(abs(obj$gr()))
obj$report()$lp_prior
obj$report()$lp_penalty
obj$report()$lp_rec
sum(obj$report()$lp_cpue)
sum(obj$report()$lp_lf)
sum(obj$report()$lp_wf)
```

Compare initial and estimated parameter values:

``` r
get_par_table(obj, parameters, map, lower = Lwr, upper = Upr, grad_tol = 1e-2, digits = 2L)
#> outer mgc:  4.072735e-12
#>          par init   est  lwr  upr       gr gr_chk bd_chk
#> 1     log_B0   20 14.00  0.0 22.0 -4.1e-12     OK     OK
#> 2 log_cpue_q    0 -0.02 -2.3  2.3 -2.3e-12     OK     OK
```

### Diagnostics

Check that all parameters are estimable using the `check_estimability`
function:

``` r
check_estimability(obj = obj)
#> outer mgc:  4.072735e-12 
#> outer mgc:  7.301959 
#> outer mgc:  7.417775 
#> outer mgc:  7.425739 
#> outer mgc:  7.425739 
#> outer mgc:  0.03705092 
#> outer mgc:  0.03701259 
#> outer mgc:  0.03777039 
#> outer mgc:  0.03773221 
#> outer mgc:  0.03855943 
#> outer mgc:  0.03852048 
#> outer mgc:  0.04029141 
#> outer mgc:  0.04025065 
#> outer mgc:  0.04058785 
#> outer mgc:  0.0405457 
#> outer mgc:  0.04280243 
#> outer mgc:  0.04275888 
#> outer mgc:  0.04312595 
#> outer mgc:  0.04308214 
#> outer mgc:  0.04344454 
#> outer mgc:  0.0434004 
#> outer mgc:  0.04228394 
#> outer mgc:  0.04224013 
#> outer mgc:  0.04352656 
#> outer mgc:  0.04348221 
#> outer mgc:  0.04347253 
#> outer mgc:  0.04342828 
#> outer mgc:  0.04371646 
#> outer mgc:  0.04367193 
#> outer mgc:  0.04074535 
#> outer mgc:  0.0407036 
#> outer mgc:  0.04109553 
#> outer mgc:  0.04105475 
#> outer mgc:  0.03835994 
#> outer mgc:  0.03832198 
#> outer mgc:  0.03608784 
#> outer mgc:  0.03605214 
#> outer mgc:  0.03193525 
#> outer mgc:  0.03190177 
#> outer mgc:  0.03167213 
#> outer mgc:  0.0316406 
#> outer mgc:  0.02989747 
#> outer mgc:  0.02986787 
#> outer mgc:  0.02850442 
#> outer mgc:  0.02847622 
#> outer mgc:  0.02533635 
#> outer mgc:  0.02530866 
#> outer mgc:  0.02648068 
#> outer mgc:  0.02645423 
#> outer mgc:  0.02645692 
#> outer mgc:  0.02643076 
#> outer mgc:  0.02673562 
#> outer mgc:  0.02670921 
#> outer mgc:  0.02531674 
#> outer mgc:  0.02528939 
#> outer mgc:  0.02713217 
#> outer mgc:  0.02710517 
#> outer mgc:  0.02757828 
#> outer mgc:  0.02755108 
#> outer mgc:  0.02797858 
#> outer mgc:  0.027951 
#> outer mgc:  0.02763535 
#> outer mgc:  0.02760712 
#> outer mgc:  0.02861779 
#> outer mgc:  0.02858954 
#> outer mgc:  0.02907004 
#> outer mgc:  0.02904146 
#> outer mgc:  0.02974009 
#> outer mgc:  0.02971088 
#> outer mgc:  0.03029216 
#> outer mgc:  0.03026184 
#> outer mgc:  0.03149609 
#> outer mgc:  0.03146479 
#> outer mgc:  0.03272599 
#> outer mgc:  0.03269393 
#> outer mgc:  0.03335334 
#> outer mgc:  0.03332072 
#> outer mgc:  0.03287251 
#> outer mgc:  0.03283961 
#> outer mgc:  0.03248496 
#> outer mgc:  0.03245265 
#> outer mgc:  0.03205766 
#> outer mgc:  0.03202629 
#> outer mgc:  0.03148134 
#> outer mgc:  0.03145055 
#> outer mgc:  0.03102964 
#> outer mgc:  0.03099895 
#> outer mgc:  0.03116143 
#> outer mgc:  0.03113071 
#> outer mgc:  0.03166157 
#> outer mgc:  0.03163062 
#> outer mgc:  0.03253915 
#> outer mgc:  0.03250738 
#> outer mgc:  0.033387 
#> outer mgc:  0.03335367 
#> outer mgc:  0.03490109 
#> outer mgc:  0.03486654 
#> outer mgc:  0.03591642 
#> outer mgc:  0.03588142 
#> outer mgc:  0.03592828 
#> outer mgc:  0.03589328 
#> outer mgc:  0.03531244 
#> outer mgc:  0.03527759 
#> outer mgc:  0.03465394 
#> outer mgc:  0.03461984 
#> outer mgc:  0.03363056 
#> outer mgc:  0.03359775 
#> outer mgc:  0.03253897 
#> outer mgc:  0.03250723 
#> outer mgc:  0.03138182 
#> outer mgc:  0.03135066 
#> outer mgc:  0.03059518 
#> outer mgc:  0.03056492 
#> outer mgc:  0.02984814 
#> outer mgc:  0.02981896 
#> outer mgc:  0.02925111 
#> outer mgc:  0.02922252 
#> outer mgc:  0.02893164 
#> outer mgc:  0.02890287 
#> outer mgc:  0.0294495 
#> outer mgc:  0.02942037 
#> outer mgc:  0.03039898 
#> outer mgc:  0.03036926 
#> outer mgc:  0.03146594 
#> outer mgc:  0.03143522 
#> outer mgc:  0.03138598 
#> outer mgc:  0.03135417 
#> outer mgc:  0.03095806 
#> outer mgc:  0.03092699 
#> outer mgc:  0.02945591 
#> outer mgc:  0.02942706 
#> outer mgc:  0.02728879 
#> outer mgc:  0.02726206 
#> outer mgc:  0.02505717 
#> outer mgc:  0.02503155 
#> outer mgc:  0.02413791 
#> outer mgc:  0.02411345 
#> outer mgc:  0.02418304 
#> outer mgc:  0.02415924 
#> outer mgc:  0.02470402 
#> outer mgc:  0.0246794 
#> outer mgc:  0.02598373 
#> outer mgc:  0.02595717 
#> outer mgc:  0.02844988 
#> outer mgc:  0.02842121 
#> outer mgc:  0.0309806 
#> outer mgc:  0.03094995 
#> outer mgc:  0.03263572 
#> outer mgc:  0.03260341 
#> outer mgc:  0.03329664 
#> outer mgc:  0.03326311 
#> outer mgc:  0.03274996 
#> outer mgc:  0.03271711 
#> outer mgc:  0.03064222 
#> outer mgc:  0.03061179 
#> outer mgc:  0.02796335 
#> outer mgc:  0.02793564 
#> outer mgc:  0.02544004 
#> outer mgc:  0.02541439 
#> outer mgc:  0.02367066 
#> outer mgc:  0.02364684 
#> outer mgc:  0.02248174 
#> outer mgc:  0.02245933 
#> outer mgc:  0.02191011 
#> outer mgc:  0.02188835 
#> outer mgc:  0.02213394 
#> outer mgc:  0.02211179 
#> outer mgc:  0.02310288 
#> outer mgc:  0.02307987 
#> outer mgc:  0.02443416 
#> outer mgc:  0.02440989 
#> outer mgc:  0.02617322 
#> outer mgc:  0.02614725 
#> outer mgc:  0.02782609 
#> outer mgc:  0.02779825 
#> outer mgc:  0.02938929 
#> outer mgc:  0.02936007 
#> outer mgc:  0.03028732 
#> outer mgc:  0.03025733 
#> outer mgc:  0.03047173 
#> outer mgc:  0.03044151 
#> outer mgc:  0.03004622 
#> outer mgc:  0.03001622 
#> outer mgc:  0.02916754 
#> outer mgc:  0.02913847 
#> outer mgc:  0.0281353 
#> outer mgc:  0.0281073 
#> outer mgc:  0.0276073 
#> outer mgc:  0.02757975 
#> outer mgc:  0.02782265 
#> outer mgc:  0.02779461 
#> outer mgc:  0.02868982 
#> outer mgc:  0.02866097 
#> outer mgc:  0.02929985 
#> outer mgc:  0.02927042 
#> outer mgc:  0.02884304 
#> outer mgc:  0.028814 
#> outer mgc:  0.02642563 
#> outer mgc:  0.02639812 
#> outer mgc:  0.02449987 
#> outer mgc:  0.0244748 
#> outer mgc:  0.02287367 
#> outer mgc:  0.02285052 
#> outer mgc:  0.02184823 
#> outer mgc:  0.02182604 
#> outer mgc:  0.02103326 
#> outer mgc:  0.02101124 
#> outer mgc:  0.020828 
#> outer mgc:  0.02080656 
#> outer mgc:  0.02015474 
#> outer mgc:  0.02013409 
#> outer mgc:  0.01918406 
#> outer mgc:  0.01916441 
#> outer mgc:  0.01724616 
#> outer mgc:  0.01722754 
#> outer mgc:  0.01613423 
#> outer mgc:  0.01611744 
#> outer mgc:  0.01507217 
#> outer mgc:  0.01505673 
#> outer mgc:  0.01431762 
#> outer mgc:  0.01430289 
#> outer mgc:  0.01372884 
#> outer mgc:  0.01371407 
#> outer mgc:  0.01405542 
#> outer mgc:  0.01404072 
#> outer mgc:  0.01454303 
#> outer mgc:  0.01452784 
#> outer mgc:  0.01530579 
#> outer mgc:  0.01528968 
#> outer mgc:  0.01606874 
#> outer mgc:  0.01605145 
#> outer mgc:  0.01686147 
#> outer mgc:  0.01684333 
#> outer mgc:  0.01710547 
#> outer mgc:  0.01708688 
#> outer mgc:  0.01700595 
#> outer mgc:  0.01698723 
#> outer mgc:  0.01647916 
#> outer mgc:  0.01646061 
#> outer mgc:  0.01565848 
#> outer mgc:  0.01564079 
#> outer mgc:  0.01462556 
#> outer mgc:  0.01460893 
#> outer mgc:  0.01370259 
#> outer mgc:  0.01368679 
#> outer mgc:  0.01283882 
#> outer mgc:  0.01282361 
#> outer mgc:  0.01236047 
#> outer mgc:  0.01234579 
#> outer mgc:  0.01190606 
#> outer mgc:  0.01189174 
#> outer mgc:  0.01141919 
#> outer mgc:  0.01140522 
#> outer mgc:  0.01082422 
#> outer mgc:  0.01081057 
#> outer mgc:  0.01058847 
#> outer mgc:  0.010575 
#> outer mgc:  0.01043105 
#> outer mgc:  0.01041747 
#> outer mgc:  0.01014746 
#> outer mgc:  0.01013383 
#> outer mgc:  0.009685618 
#> outer mgc:  0.009672127 
#> outer mgc:  0.009312784 
#> outer mgc:  0.009299263 
#> outer mgc:  0.008912726 
#> outer mgc:  0.008899084 
#> outer mgc:  0.008184648 
#> outer mgc:  0.008171483 
#> outer mgc:  0.007197504 
#> outer mgc:  0.007185436 
#> outer mgc:  0.006214566 
#> outer mgc:  0.00620384 
#> outer mgc:  0.005658703 
#> outer mgc:  0.00565157 
#> outer mgc:  0.005565039 
#> outer mgc:  0.005558498 
#> outer mgc:  0.005464202 
#> outer mgc:  0.005458053 
#> outer mgc:  0.005459641 
#> outer mgc:  0.005453687 
#> outer mgc:  0.005512257 
#> outer mgc:  0.005506202 
#> outer mgc:  0.005951965 
#> outer mgc:  0.005948283 
#> outer mgc:  0.006568977 
#> outer mgc:  0.006564679 
#> outer mgc:  0.007521996 
#> outer mgc:  0.007517236 
#> outer mgc:  0.008950055 
#> outer mgc:  0.008944883 
#> outer mgc:  0.01043485 
#> outer mgc:  0.01042924 
#> outer mgc:  0.01186668 
#> outer mgc:  0.01186071 
#> outer mgc:  0.01261979 
#> outer mgc:  0.01261368 
#> outer mgc:  0.01239226 
#> outer mgc:  0.012386 
#> outer mgc:  0.01199413 
#> outer mgc:  0.01198746 
#> outer mgc:  0.01185217 
#> outer mgc:  0.01184475 
#> outer mgc:  0.01225857 
#> outer mgc:  0.01225065 
#> outer mgc:  0.01298552 
#> outer mgc:  0.012977 
#> outer mgc:  0.01400607 
#> outer mgc:  0.01399681 
#> outer mgc:  0.01512083 
#> outer mgc:  0.01511059 
#> outer mgc:  0.0161317 
#> outer mgc:  0.01612097 
#> outer mgc:  0.01823648 
#> outer mgc:  0.01823396 
#> outer mgc:  0.01972271 
#> outer mgc:  0.01971924 
#> outer mgc:  0.02004934 
#> outer mgc:  0.020044 
#> outer mgc:  0.02000971 
#> outer mgc:  0.02000286 
#> outer mgc:  0.0198004 
#> outer mgc:  0.01979208 
#> outer mgc:  0.01998454 
#> outer mgc:  0.01997522 
#> outer mgc:  0.02086931 
#> outer mgc:  0.02085933 
#> outer mgc:  0.0219632 
#> outer mgc:  0.02195245 
#> outer mgc:  0.02321409 
#> outer mgc:  0.02320266 
#> outer mgc:  0.02416807 
#> outer mgc:  0.02415585 
#> outer mgc:  0.02420947 
#> outer mgc:  0.02419606 
#> outer mgc:  0.02415429 
#> outer mgc:  0.02414035 
#> outer mgc:  0.02429449 
#> outer mgc:  0.02428002 
#> outer mgc:  0.02527904 
#> outer mgc:  0.02526371 
#> outer mgc:  0.02726969 
#> outer mgc:  0.02725253 
#> outer mgc:  0.03057398 
#> outer mgc:  0.03055479 
#> outer mgc:  0.03640508 
#> outer mgc:  0.03638394 
#> outer mgc:  0.04475474 
#> outer mgc:  0.04473133 
#> outer mgc:  0.05563076 
#> outer mgc:  0.05560552 
#> outer mgc:  0.06843682 
#> outer mgc:  0.06841202 
#> outer mgc:  0.07521148 
#> outer mgc:  0.07518597 
#> outer mgc:  0.07407738 
#> outer mgc:  0.07404968 
#> outer mgc:  0.0674127 
#> outer mgc:  0.06738145 
#> outer mgc:  0.06054421 
#> outer mgc:  0.06051058 
#> outer mgc:  0.05713168 
#> outer mgc:  0.05709774 
#> outer mgc:  0.0577185 
#> outer mgc:  0.0576851 
#> outer mgc:  0.0607472 
#> outer mgc:  0.0607114 
#> outer mgc:  0.06511649 
#> outer mgc:  0.06507893 
#> outer mgc:  0.07434123 
#> outer mgc:  0.07430351 
#> outer mgc:  0.08572162 
#> outer mgc:  0.08568225 
#> outer mgc:  0.09931663 
#> outer mgc:  0.09927607 
#> outer mgc:  0.1118563 
#> outer mgc:  0.1118172 
#> outer mgc:  0.1171744 
#> outer mgc:  0.1171334 
#> outer mgc:  0.1119395 
#> outer mgc:  0.1118911 
#> outer mgc:  0.1083053 
#> outer mgc:  0.1082523 
#> outer mgc:  0.1085015 
#> outer mgc:  0.1084476 
#> outer mgc:  0.1131364 
#> outer mgc:  0.1130823 
#> outer mgc:  0.1168027 
#> outer mgc:  0.1167441 
#> outer mgc:  0.12824 
#> outer mgc:  0.1281782 
#> outer mgc:  0.1448289 
#> outer mgc:  0.1447682 
#> outer mgc:  0.1488663 
#> outer mgc:  0.148807 
#> outer mgc:  0.1331433 
#> outer mgc:  0.1330816 
#> outer mgc:  0.1160581 
#> outer mgc:  0.115995 
#> outer mgc:  0.1032122 
#> outer mgc:  0.1031486 
#> outer mgc:  0.09707163 
#> outer mgc:  0.09700854 
#> outer mgc:  0.1017144 
#> outer mgc:  0.1016516 
#> outer mgc:  0.1115651 
#> outer mgc:  0.111499 
#> outer mgc:  0.1248447 
#> outer mgc:  0.1247761 
#> outer mgc:  0.1382938 
#> outer mgc:  0.1382258 
#> outer mgc:  0.1434027 
#> outer mgc:  0.1433369 
#> outer mgc:  0.1357267 
#> outer mgc:  0.1356579 
#> outer mgc:  0.1304035 
#> outer mgc:  0.1303344 
#> outer mgc:  0.1326064 
#> outer mgc:  0.1325377 
#> outer mgc:  0.1399934 
#> outer mgc:  0.1399237 
#> outer mgc:  0.1478892 
#> outer mgc:  0.1478136 
#> outer mgc:  0.1607946 
#> outer mgc:  0.160719 
#> outer mgc:  0.1747331 
#> outer mgc:  0.1746599 
#> outer mgc:  0.1797385 
#> outer mgc:  0.1796667 
#> outer mgc:  0.1718467 
#> outer mgc:  0.1717704 
#> outer mgc:  0.1622535 
#> outer mgc:  0.1621772 
#> outer mgc:  0.1580853 
#> outer mgc:  0.1580109 
#> outer mgc:  0.1579289 
#> outer mgc:  0.1578556 
#> outer mgc:  0.1591844 
#> outer mgc:  0.1591053 
#> outer mgc:  0.1640026 
#> outer mgc:  0.1639217 
#> outer mgc:  0.1801885 
#> outer mgc:  0.180115 
#> outer mgc:  0.196168 
#> outer mgc:  0.1961033 
#> outer mgc:  0.1963476 
#> outer mgc:  0.1962794 
#> outer mgc:  0.1890751 
#> outer mgc:  0.1890018 
#> outer mgc:  0.1847872 
#> outer mgc:  0.1847152 
#> outer mgc:  0.1826456 
#> outer mgc:  0.1825769 
#> outer mgc:  0.1778515 
#> outer mgc:  0.1777826 
#> outer mgc:  0.1690758 
#> outer mgc:  0.1690048 
#> outer mgc:  0.1631956 
#> outer mgc:  0.1631231 
#> outer mgc:  0.163937 
#> outer mgc:  0.1638659 
#> outer mgc:  0.1739192 
#> outer mgc:  0.1738504 
#> outer mgc:  0.1800206 
#> outer mgc:  0.1799495 
#> outer mgc:  0.1942705 
#> outer mgc:  0.1942062 
#> outer mgc:  0.2003827 
#> outer mgc:  0.2003342 
#> outer mgc:  0.1823035 
#> outer mgc:  0.1822476 
#> outer mgc:  0.1519063 
#> outer mgc:  0.151837 
#> outer mgc:  0.1329732 
#> outer mgc:  0.1329013 
#> outer mgc:  0.1273073 
#> outer mgc:  0.1272425 
#> outer mgc:  0.1270056 
#> outer mgc:  0.1269448 
#> outer mgc:  0.1181345 
#> outer mgc:  0.1180691 
#> outer mgc:  0.1215816 
#> outer mgc:  0.1215165 
#> outer mgc:  0.1341553 
#> outer mgc:  0.1340888 
#> outer mgc:  0.1763906 
#> outer mgc:  0.1763399 
#> outer mgc:  0.2123527 
#> outer mgc:  0.2123155 
#> outer mgc:  0.2587935 
#> outer mgc:  0.2588034 
#> outer mgc:  0.1943547 
#> outer mgc:  0.1943242 
#> outer mgc:  0.1489883 
#> outer mgc:  0.148937 
#> outer mgc:  0.1166996 
#> outer mgc:  0.116642 
#> outer mgc:  0.1013811 
#> outer mgc:  0.1013279 
#> outer mgc:  0.09537271 
#> outer mgc:  0.09532768 
#> outer mgc:  0.09008129 
#> outer mgc:  0.09003705 
#> outer mgc:  0.08350529 
#> outer mgc:  0.08346109 
#> outer mgc:  0.08242265 
#> outer mgc:  0.08238205 
#> outer mgc:  0.08302936 
#> outer mgc:  0.08299201 
#> outer mgc:  0.08797039 
#> outer mgc:  0.08794074 
#> outer mgc:  0.08446968 
#> outer mgc:  0.08443955 
#> outer mgc:  0.08047273 
#> outer mgc:  0.08044473 
#> outer mgc:  0.07381704 
#> outer mgc:  0.07379129 
#> outer mgc:  0.06453203 
#> outer mgc:  0.06450665 
#> outer mgc:  0.04994719 
#> outer mgc:  0.0499207 
#> outer mgc:  0.03847386 
#> outer mgc:  0.03845168 
#> outer mgc:  0.02705258 
#> outer mgc:  0.02703592 
#> outer mgc:  0.01685738 
#> outer mgc:  0.01684572 
#> outer mgc:  0.009547628 
#> outer mgc:  0.009540424 
#> outer mgc:  0.00457244 
#> outer mgc:  0.004568829 
#> outer mgc:  0.002775735 
#> outer mgc:  0.002775718 
#> outer mgc:  0.002774406 
#> outer mgc:  0.002774408 
#> outer mgc:  0.002776984 
#> outer mgc:  0.002776985 
#> outer mgc:  0.002777671 
#> outer mgc:  0.002777671 
#> outer mgc:  0.002777778 
#> outer mgc:  0.002777778
```

Look for any parameters with high correlations (absolute correlation \>
0.95) using the `get_cor_pairs` function:

``` r
get_cor_pairs(obj = obj, threshold = 0.95)
#> [1] "No pairs show strong positive or negative correlation relative to 0.95 threshold."
```

Calculate standard deviations of all model parameters:

``` r
Report <- sdreport(obj)
#> outer mgc:  4.072735e-12 
#> outer mgc:  7.301959 
#> outer mgc:  7.417775 
#> outer mgc:  7.425739 
#> outer mgc:  7.425739 
#> outer mgc:  0.03705092 
#> outer mgc:  0.03701259 
#> outer mgc:  0.03777039 
#> outer mgc:  0.03773221 
#> outer mgc:  0.03855943 
#> outer mgc:  0.03852048 
#> outer mgc:  0.04029141 
#> outer mgc:  0.04025065 
#> outer mgc:  0.04058785 
#> outer mgc:  0.0405457 
#> outer mgc:  0.04280243 
#> outer mgc:  0.04275888 
#> outer mgc:  0.04312595 
#> outer mgc:  0.04308214 
#> outer mgc:  0.04344454 
#> outer mgc:  0.0434004 
#> outer mgc:  0.04228394 
#> outer mgc:  0.04224013 
#> outer mgc:  0.04352656 
#> outer mgc:  0.04348221 
#> outer mgc:  0.04347253 
#> outer mgc:  0.04342828 
#> outer mgc:  0.04371646 
#> outer mgc:  0.04367193 
#> outer mgc:  0.04074535 
#> outer mgc:  0.0407036 
#> outer mgc:  0.04109553 
#> outer mgc:  0.04105475 
#> outer mgc:  0.03835994 
#> outer mgc:  0.03832198 
#> outer mgc:  0.03608784 
#> outer mgc:  0.03605214 
#> outer mgc:  0.03193525 
#> outer mgc:  0.03190177 
#> outer mgc:  0.03167213 
#> outer mgc:  0.0316406 
#> outer mgc:  0.02989747 
#> outer mgc:  0.02986787 
#> outer mgc:  0.02850442 
#> outer mgc:  0.02847622 
#> outer mgc:  0.02533635 
#> outer mgc:  0.02530866 
#> outer mgc:  0.02648068 
#> outer mgc:  0.02645423 
#> outer mgc:  0.02645692 
#> outer mgc:  0.02643076 
#> outer mgc:  0.02673562 
#> outer mgc:  0.02670921 
#> outer mgc:  0.02531674 
#> outer mgc:  0.02528939 
#> outer mgc:  0.02713217 
#> outer mgc:  0.02710517 
#> outer mgc:  0.02757828 
#> outer mgc:  0.02755108 
#> outer mgc:  0.02797858 
#> outer mgc:  0.027951 
#> outer mgc:  0.02763535 
#> outer mgc:  0.02760712 
#> outer mgc:  0.02861779 
#> outer mgc:  0.02858954 
#> outer mgc:  0.02907004 
#> outer mgc:  0.02904146 
#> outer mgc:  0.02974009 
#> outer mgc:  0.02971088 
#> outer mgc:  0.03029216 
#> outer mgc:  0.03026184 
#> outer mgc:  0.03149609 
#> outer mgc:  0.03146479 
#> outer mgc:  0.03272599 
#> outer mgc:  0.03269393 
#> outer mgc:  0.03335334 
#> outer mgc:  0.03332072 
#> outer mgc:  0.03287251 
#> outer mgc:  0.03283961 
#> outer mgc:  0.03248496 
#> outer mgc:  0.03245265 
#> outer mgc:  0.03205766 
#> outer mgc:  0.03202629 
#> outer mgc:  0.03148134 
#> outer mgc:  0.03145055 
#> outer mgc:  0.03102964 
#> outer mgc:  0.03099895 
#> outer mgc:  0.03116143 
#> outer mgc:  0.03113071 
#> outer mgc:  0.03166157 
#> outer mgc:  0.03163062 
#> outer mgc:  0.03253915 
#> outer mgc:  0.03250738 
#> outer mgc:  0.033387 
#> outer mgc:  0.03335367 
#> outer mgc:  0.03490109 
#> outer mgc:  0.03486654 
#> outer mgc:  0.03591642 
#> outer mgc:  0.03588142 
#> outer mgc:  0.03592828 
#> outer mgc:  0.03589328 
#> outer mgc:  0.03531244 
#> outer mgc:  0.03527759 
#> outer mgc:  0.03465394 
#> outer mgc:  0.03461984 
#> outer mgc:  0.03363056 
#> outer mgc:  0.03359775 
#> outer mgc:  0.03253897 
#> outer mgc:  0.03250723 
#> outer mgc:  0.03138182 
#> outer mgc:  0.03135066 
#> outer mgc:  0.03059518 
#> outer mgc:  0.03056492 
#> outer mgc:  0.02984814 
#> outer mgc:  0.02981896 
#> outer mgc:  0.02925111 
#> outer mgc:  0.02922252 
#> outer mgc:  0.02893164 
#> outer mgc:  0.02890287 
#> outer mgc:  0.0294495 
#> outer mgc:  0.02942037 
#> outer mgc:  0.03039898 
#> outer mgc:  0.03036926 
#> outer mgc:  0.03146594 
#> outer mgc:  0.03143522 
#> outer mgc:  0.03138598 
#> outer mgc:  0.03135417 
#> outer mgc:  0.03095806 
#> outer mgc:  0.03092699 
#> outer mgc:  0.02945591 
#> outer mgc:  0.02942706 
#> outer mgc:  0.02728879 
#> outer mgc:  0.02726206 
#> outer mgc:  0.02505717 
#> outer mgc:  0.02503155 
#> outer mgc:  0.02413791 
#> outer mgc:  0.02411345 
#> outer mgc:  0.02418304 
#> outer mgc:  0.02415924 
#> outer mgc:  0.02470402 
#> outer mgc:  0.0246794 
#> outer mgc:  0.02598373 
#> outer mgc:  0.02595717 
#> outer mgc:  0.02844988 
#> outer mgc:  0.02842121 
#> outer mgc:  0.0309806 
#> outer mgc:  0.03094995 
#> outer mgc:  0.03263572 
#> outer mgc:  0.03260341 
#> outer mgc:  0.03329664 
#> outer mgc:  0.03326311 
#> outer mgc:  0.03274996 
#> outer mgc:  0.03271711 
#> outer mgc:  0.03064222 
#> outer mgc:  0.03061179 
#> outer mgc:  0.02796335 
#> outer mgc:  0.02793564 
#> outer mgc:  0.02544004 
#> outer mgc:  0.02541439 
#> outer mgc:  0.02367066 
#> outer mgc:  0.02364684 
#> outer mgc:  0.02248174 
#> outer mgc:  0.02245933 
#> outer mgc:  0.02191011 
#> outer mgc:  0.02188835 
#> outer mgc:  0.02213394 
#> outer mgc:  0.02211179 
#> outer mgc:  0.02310288 
#> outer mgc:  0.02307987 
#> outer mgc:  0.02443416 
#> outer mgc:  0.02440989 
#> outer mgc:  0.02617322 
#> outer mgc:  0.02614725 
#> outer mgc:  0.02782609 
#> outer mgc:  0.02779825 
#> outer mgc:  0.02938929 
#> outer mgc:  0.02936007 
#> outer mgc:  0.03028732 
#> outer mgc:  0.03025733 
#> outer mgc:  0.03047173 
#> outer mgc:  0.03044151 
#> outer mgc:  0.03004622 
#> outer mgc:  0.03001622 
#> outer mgc:  0.02916754 
#> outer mgc:  0.02913847 
#> outer mgc:  0.0281353 
#> outer mgc:  0.0281073 
#> outer mgc:  0.0276073 
#> outer mgc:  0.02757975 
#> outer mgc:  0.02782265 
#> outer mgc:  0.02779461 
#> outer mgc:  0.02868982 
#> outer mgc:  0.02866097 
#> outer mgc:  0.02929985 
#> outer mgc:  0.02927042 
#> outer mgc:  0.02884304 
#> outer mgc:  0.028814 
#> outer mgc:  0.02642563 
#> outer mgc:  0.02639812 
#> outer mgc:  0.02449987 
#> outer mgc:  0.0244748 
#> outer mgc:  0.02287367 
#> outer mgc:  0.02285052 
#> outer mgc:  0.02184823 
#> outer mgc:  0.02182604 
#> outer mgc:  0.02103326 
#> outer mgc:  0.02101124 
#> outer mgc:  0.020828 
#> outer mgc:  0.02080656 
#> outer mgc:  0.02015474 
#> outer mgc:  0.02013409 
#> outer mgc:  0.01918406 
#> outer mgc:  0.01916441 
#> outer mgc:  0.01724616 
#> outer mgc:  0.01722754 
#> outer mgc:  0.01613423 
#> outer mgc:  0.01611744 
#> outer mgc:  0.01507217 
#> outer mgc:  0.01505673 
#> outer mgc:  0.01431762 
#> outer mgc:  0.01430289 
#> outer mgc:  0.01372884 
#> outer mgc:  0.01371407 
#> outer mgc:  0.01405542 
#> outer mgc:  0.01404072 
#> outer mgc:  0.01454303 
#> outer mgc:  0.01452784 
#> outer mgc:  0.01530579 
#> outer mgc:  0.01528968 
#> outer mgc:  0.01606874 
#> outer mgc:  0.01605145 
#> outer mgc:  0.01686147 
#> outer mgc:  0.01684333 
#> outer mgc:  0.01710547 
#> outer mgc:  0.01708688 
#> outer mgc:  0.01700595 
#> outer mgc:  0.01698723 
#> outer mgc:  0.01647916 
#> outer mgc:  0.01646061 
#> outer mgc:  0.01565848 
#> outer mgc:  0.01564079 
#> outer mgc:  0.01462556 
#> outer mgc:  0.01460893 
#> outer mgc:  0.01370259 
#> outer mgc:  0.01368679 
#> outer mgc:  0.01283882 
#> outer mgc:  0.01282361 
#> outer mgc:  0.01236047 
#> outer mgc:  0.01234579 
#> outer mgc:  0.01190606 
#> outer mgc:  0.01189174 
#> outer mgc:  0.01141919 
#> outer mgc:  0.01140522 
#> outer mgc:  0.01082422 
#> outer mgc:  0.01081057 
#> outer mgc:  0.01058847 
#> outer mgc:  0.010575 
#> outer mgc:  0.01043105 
#> outer mgc:  0.01041747 
#> outer mgc:  0.01014746 
#> outer mgc:  0.01013383 
#> outer mgc:  0.009685618 
#> outer mgc:  0.009672127 
#> outer mgc:  0.009312784 
#> outer mgc:  0.009299263 
#> outer mgc:  0.008912726 
#> outer mgc:  0.008899084 
#> outer mgc:  0.008184648 
#> outer mgc:  0.008171483 
#> outer mgc:  0.007197504 
#> outer mgc:  0.007185436 
#> outer mgc:  0.006214566 
#> outer mgc:  0.00620384 
#> outer mgc:  0.005658703 
#> outer mgc:  0.00565157 
#> outer mgc:  0.005565039 
#> outer mgc:  0.005558498 
#> outer mgc:  0.005464202 
#> outer mgc:  0.005458053 
#> outer mgc:  0.005459641 
#> outer mgc:  0.005453687 
#> outer mgc:  0.005512257 
#> outer mgc:  0.005506202 
#> outer mgc:  0.005951965 
#> outer mgc:  0.005948283 
#> outer mgc:  0.006568977 
#> outer mgc:  0.006564679 
#> outer mgc:  0.007521996 
#> outer mgc:  0.007517236 
#> outer mgc:  0.008950055 
#> outer mgc:  0.008944883 
#> outer mgc:  0.01043485 
#> outer mgc:  0.01042924 
#> outer mgc:  0.01186668 
#> outer mgc:  0.01186071 
#> outer mgc:  0.01261979 
#> outer mgc:  0.01261368 
#> outer mgc:  0.01239226 
#> outer mgc:  0.012386 
#> outer mgc:  0.01199413 
#> outer mgc:  0.01198746 
#> outer mgc:  0.01185217 
#> outer mgc:  0.01184475 
#> outer mgc:  0.01225857 
#> outer mgc:  0.01225065 
#> outer mgc:  0.01298552 
#> outer mgc:  0.012977 
#> outer mgc:  0.01400607 
#> outer mgc:  0.01399681 
#> outer mgc:  0.01512083 
#> outer mgc:  0.01511059 
#> outer mgc:  0.0161317 
#> outer mgc:  0.01612097 
#> outer mgc:  0.01823648 
#> outer mgc:  0.01823396 
#> outer mgc:  0.01972271 
#> outer mgc:  0.01971924 
#> outer mgc:  0.02004934 
#> outer mgc:  0.020044 
#> outer mgc:  0.02000971 
#> outer mgc:  0.02000286 
#> outer mgc:  0.0198004 
#> outer mgc:  0.01979208 
#> outer mgc:  0.01998454 
#> outer mgc:  0.01997522 
#> outer mgc:  0.02086931 
#> outer mgc:  0.02085933 
#> outer mgc:  0.0219632 
#> outer mgc:  0.02195245 
#> outer mgc:  0.02321409 
#> outer mgc:  0.02320266 
#> outer mgc:  0.02416807 
#> outer mgc:  0.02415585 
#> outer mgc:  0.02420947 
#> outer mgc:  0.02419606 
#> outer mgc:  0.02415429 
#> outer mgc:  0.02414035 
#> outer mgc:  0.02429449 
#> outer mgc:  0.02428002 
#> outer mgc:  0.02527904 
#> outer mgc:  0.02526371 
#> outer mgc:  0.02726969 
#> outer mgc:  0.02725253 
#> outer mgc:  0.03057398 
#> outer mgc:  0.03055479 
#> outer mgc:  0.03640508 
#> outer mgc:  0.03638394 
#> outer mgc:  0.04475474 
#> outer mgc:  0.04473133 
#> outer mgc:  0.05563076 
#> outer mgc:  0.05560552 
#> outer mgc:  0.06843682 
#> outer mgc:  0.06841202 
#> outer mgc:  0.07521148 
#> outer mgc:  0.07518597 
#> outer mgc:  0.07407738 
#> outer mgc:  0.07404968 
#> outer mgc:  0.0674127 
#> outer mgc:  0.06738145 
#> outer mgc:  0.06054421 
#> outer mgc:  0.06051058 
#> outer mgc:  0.05713168 
#> outer mgc:  0.05709774 
#> outer mgc:  0.0577185 
#> outer mgc:  0.0576851 
#> outer mgc:  0.0607472 
#> outer mgc:  0.0607114 
#> outer mgc:  0.06511649 
#> outer mgc:  0.06507893 
#> outer mgc:  0.07434123 
#> outer mgc:  0.07430351 
#> outer mgc:  0.08572162 
#> outer mgc:  0.08568225 
#> outer mgc:  0.09931663 
#> outer mgc:  0.09927607 
#> outer mgc:  0.1118563 
#> outer mgc:  0.1118172 
#> outer mgc:  0.1171744 
#> outer mgc:  0.1171334 
#> outer mgc:  0.1119395 
#> outer mgc:  0.1118911 
#> outer mgc:  0.1083053 
#> outer mgc:  0.1082523 
#> outer mgc:  0.1085015 
#> outer mgc:  0.1084476 
#> outer mgc:  0.1131364 
#> outer mgc:  0.1130823 
#> outer mgc:  0.1168027 
#> outer mgc:  0.1167441 
#> outer mgc:  0.12824 
#> outer mgc:  0.1281782 
#> outer mgc:  0.1448289 
#> outer mgc:  0.1447682 
#> outer mgc:  0.1488663 
#> outer mgc:  0.148807 
#> outer mgc:  0.1331433 
#> outer mgc:  0.1330816 
#> outer mgc:  0.1160581 
#> outer mgc:  0.115995 
#> outer mgc:  0.1032122 
#> outer mgc:  0.1031486 
#> outer mgc:  0.09707163 
#> outer mgc:  0.09700854 
#> outer mgc:  0.1017144 
#> outer mgc:  0.1016516 
#> outer mgc:  0.1115651 
#> outer mgc:  0.111499 
#> outer mgc:  0.1248447 
#> outer mgc:  0.1247761 
#> outer mgc:  0.1382938 
#> outer mgc:  0.1382258 
#> outer mgc:  0.1434027 
#> outer mgc:  0.1433369 
#> outer mgc:  0.1357267 
#> outer mgc:  0.1356579 
#> outer mgc:  0.1304035 
#> outer mgc:  0.1303344 
#> outer mgc:  0.1326064 
#> outer mgc:  0.1325377 
#> outer mgc:  0.1399934 
#> outer mgc:  0.1399237 
#> outer mgc:  0.1478892 
#> outer mgc:  0.1478136 
#> outer mgc:  0.1607946 
#> outer mgc:  0.160719 
#> outer mgc:  0.1747331 
#> outer mgc:  0.1746599 
#> outer mgc:  0.1797385 
#> outer mgc:  0.1796667 
#> outer mgc:  0.1718467 
#> outer mgc:  0.1717704 
#> outer mgc:  0.1622535 
#> outer mgc:  0.1621772 
#> outer mgc:  0.1580853 
#> outer mgc:  0.1580109 
#> outer mgc:  0.1579289 
#> outer mgc:  0.1578556 
#> outer mgc:  0.1591844 
#> outer mgc:  0.1591053 
#> outer mgc:  0.1640026 
#> outer mgc:  0.1639217 
#> outer mgc:  0.1801885 
#> outer mgc:  0.180115 
#> outer mgc:  0.196168 
#> outer mgc:  0.1961033 
#> outer mgc:  0.1963476 
#> outer mgc:  0.1962794 
#> outer mgc:  0.1890751 
#> outer mgc:  0.1890018 
#> outer mgc:  0.1847872 
#> outer mgc:  0.1847152 
#> outer mgc:  0.1826456 
#> outer mgc:  0.1825769 
#> outer mgc:  0.1778515 
#> outer mgc:  0.1777826 
#> outer mgc:  0.1690758 
#> outer mgc:  0.1690048 
#> outer mgc:  0.1631956 
#> outer mgc:  0.1631231 
#> outer mgc:  0.163937 
#> outer mgc:  0.1638659 
#> outer mgc:  0.1739192 
#> outer mgc:  0.1738504 
#> outer mgc:  0.1800206 
#> outer mgc:  0.1799495 
#> outer mgc:  0.1942705 
#> outer mgc:  0.1942062 
#> outer mgc:  0.2003827 
#> outer mgc:  0.2003342 
#> outer mgc:  0.1823035 
#> outer mgc:  0.1822476 
#> outer mgc:  0.1519063 
#> outer mgc:  0.151837 
#> outer mgc:  0.1329732 
#> outer mgc:  0.1329013 
#> outer mgc:  0.1273073 
#> outer mgc:  0.1272425 
#> outer mgc:  0.1270056 
#> outer mgc:  0.1269448 
#> outer mgc:  0.1181345 
#> outer mgc:  0.1180691 
#> outer mgc:  0.1215816 
#> outer mgc:  0.1215165 
#> outer mgc:  0.1341553 
#> outer mgc:  0.1340888 
#> outer mgc:  0.1763906 
#> outer mgc:  0.1763399 
#> outer mgc:  0.2123527 
#> outer mgc:  0.2123155 
#> outer mgc:  0.2587935 
#> outer mgc:  0.2588034 
#> outer mgc:  0.1943547 
#> outer mgc:  0.1943242 
#> outer mgc:  0.1489883 
#> outer mgc:  0.148937 
#> outer mgc:  0.1166996 
#> outer mgc:  0.116642 
#> outer mgc:  0.1013811 
#> outer mgc:  0.1013279 
#> outer mgc:  0.09537271 
#> outer mgc:  0.09532768 
#> outer mgc:  0.09008129 
#> outer mgc:  0.09003705 
#> outer mgc:  0.08350529 
#> outer mgc:  0.08346109 
#> outer mgc:  0.08242265 
#> outer mgc:  0.08238205 
#> outer mgc:  0.08302936 
#> outer mgc:  0.08299201 
#> outer mgc:  0.08797039 
#> outer mgc:  0.08794074 
#> outer mgc:  0.08446968 
#> outer mgc:  0.08443955 
#> outer mgc:  0.08047273 
#> outer mgc:  0.08044473 
#> outer mgc:  0.07381704 
#> outer mgc:  0.07379129 
#> outer mgc:  0.06453203 
#> outer mgc:  0.06450665 
#> outer mgc:  0.04994719 
#> outer mgc:  0.0499207 
#> outer mgc:  0.03847386 
#> outer mgc:  0.03845168 
#> outer mgc:  0.02705258 
#> outer mgc:  0.02703592 
#> outer mgc:  0.01685738 
#> outer mgc:  0.01684572 
#> outer mgc:  0.009547628 
#> outer mgc:  0.009540424 
#> outer mgc:  0.00457244 
#> outer mgc:  0.004568829 
#> outer mgc:  0.002775735 
#> outer mgc:  0.002775718 
#> outer mgc:  0.002774406 
#> outer mgc:  0.002774408 
#> outer mgc:  0.002776984 
#> outer mgc:  0.002776985 
#> outer mgc:  0.002777671 
#> outer mgc:  0.002777671 
#> outer mgc:  0.002777778 
#> outer mgc:  0.002777778
```

### Data fits

Inspect predicted vs observed length compositions:

``` r
if(data$lf_switch == 0){
  cat("Length composition likelihood is switched off (lf_switch = 0), so no predicted length compositions to show.\n")
} else {
  lf_rep <- obj$report()
  lf_pred <- lf_rep$lf_pred
  nbins <- data$lf_maxbin[1] - data$lf_minbin[1] + 1

  # Name the list elements and columns
  names(lf_pred) <- data$lf_fishery_f
  for (k in seq_along(lf_pred)) colnames(lf_pred[[k]]) <- data$len_mid

  # Helper: convert a flat OBS vector back to list-of-matrices (proportions)
  flat_to_list <- function(v, lf_n_f, nbins) {
    out <- vector("list", length(lf_n_f))
    offset <- 0L
    for (k in seq_along(lf_n_f)) {
      m <- matrix(v[(offset + 1):(offset + lf_n_f[k] * nbins)], nrow = lf_n_f[k], ncol = nbins, byrow = TRUE)
      # Normalise each row to proportions
      rs <- rowSums(m)
      rs[rs == 0] <- 1
      out[[k]] <- m / rs
      offset <- offset + lf_n_f[k] * nbins
    }
    out
  }

  # Split year index by fishery (must match length of lf_year)
  lf_year_list <- split(data$lf_year, data$lf_fishery)

  # ---- Predicted: long data frame ----
  df_pred <- map_dfr(seq_along(lf_pred), function(i) {
    m <- lf_pred[[i]]
    yrs <- lf_year_list[[i]]
    as.data.frame.table(m, responseName = "pred", stringsAsFactors = FALSE) %>%
      rename(id = Var1, length = Var2) %>%
      mutate(fishery = names(lf_pred)[i],
            id = as.integer(factor(id, levels = unique(id))),
            year = yrs[id],
            length = as.numeric(length))
  })

  # ---- Observed: long data frame ----
  df_obs <- data.frame(data$lf_obs_in)
  names(df_obs) <- data$len_mid
  df_obs <- df_obs %>%
    mutate(year = data$lf_year, fishery = as.character(data$lf_fishery)) %>%
    group_by(fishery) %>%
    mutate(id = row_number()) %>%
    ungroup() %>%
    pivot_longer(cols = -c(year, fishery, id), names_to = "length", values_to = "obs") %>%
    mutate(length = as.numeric(length))

  # ---- Simulations: generate n_sim draws and convert to long format ----
  n_sim <- 5
  df_sim_all <- map_dfr(seq_len(n_sim), function(s) {
    # Determine which OBS vector was active based on lf_switch
    sim <- obj$simulate()
    if (data$lf_switch == 1) {
      sim_vec <- sim$lf_obs_flat
    } else if (data$lf_switch == 2) {
      sim_vec <- sim$lf_obs_prop
    } else {
      sim_vec <- sim$lf_obs_ints
    }
    sim_list <- flat_to_list(sim_vec, data$lf_n_f, nbins)
    names(sim_list) <- data$lf_fishery_f
    for (k in seq_along(sim_list)) colnames(sim_list[[k]]) <- data$len_mid
    map_dfr(seq_along(sim_list), function(i) {
      m <- sim_list[[i]]
      yrs <- lf_year_list[[i]]
      as.data.frame.table(m, responseName = "sim", stringsAsFactors = FALSE) %>%
        rename(id = Var1, length = Var2) %>%
        mutate(fishery = names(sim_list)[i],
              id = as.integer(factor(id, levels = unique(id))),
              year = yrs[id],
              length = as.numeric(length),
              sim_id = s)
    })
  })

  # ---- Join observed + predicted ----
  df <- left_join(df_obs, df_pred, by = c("id", "fishery", "length")) %>%
    select(-year.y) %>%
    rename(year = year.x)

  # ---- Plot fishery 8 ----
  yrs_plot <- 146:188
  ggplot(df %>% filter(fishery == "8", year %in% yrs_plot), aes(x = length)) +
    geom_col(aes(y = obs), fill = "grey70", width = 2) +
    geom_line(data = df_sim_all %>% filter(fishery == "8", year %in% yrs_plot),
              aes(y = sim, group = sim_id), colour = "steelblue", alpha = 0.3, linewidth = 0.4) +
    geom_line(aes(y = pred), colour = "red3", linewidth = 0.7) +
    facet_wrap(~ year, scales = "free_y", ncol = 6) +
    labs(x = "Length (cm)", y = "Proportion",
        title = "Length composition: Fishery 8",
        subtitle = "Grey bars = observed, red = predicted, blue = simulated") +
    theme(strip.text = element_text(size = 7),
          axis.text = element_text(size = 6))

  # ---- Plot fishery 9 ----
  ggplot(df %>% filter(fishery == "9", year %in% yrs_plot), aes(x = length)) +
    geom_col(aes(y = obs), fill = "grey70", width = 2) +
    geom_line(data = df_sim_all %>% filter(fishery == "9", year %in% yrs_plot),
              aes(y = sim, group = sim_id), colour = "steelblue", alpha = 0.3, linewidth = 0.4) +
    geom_line(aes(y = pred), colour = "red3", linewidth = 0.7) +
    facet_wrap(~ year, scales = "free_y", ncol = 6) +
    labs(x = "Length (cm)", y = "Proportion",
        title = "Length composition: Fishery 9",
        subtitle = "Grey bars = observed, red = predicted, blue = simulated") +
    theme(strip.text = element_text(size = 7),
          axis.text = element_text(size = 6))
}
#> Length composition likelihood is switched off (lf_switch = 0), so no predicted length compositions to show.
```

Inspect predicted vs observed weight compositions:

``` r
if(data$wf_switch == 0){
  cat("Weight composition likelihood is switched off (wf_switch = 0), so no predicted weight compositions to show.\n")
} else {
  wf_rep  <- obj$report()
  wf_pred <- wf_rep$wf_pred
  wt_mid_trim <- data$wt_mid[data$wf_minbin[data$wf_fishery_f[1]]:
                              data$wf_maxbin[data$wf_fishery_f[1]]]

  names(wf_pred) <- data$wf_fishery_f
  for (k in seq_along(wf_pred)) colnames(wf_pred[[k]]) <- wt_mid_trim

  wf_year_list <- split(data$wf_year, data$wf_fishery)

  # ---- Predicted: long data frame ----
  df_wf_pred <- map_dfr(seq_along(wf_pred), function(i) {
    m   <- wf_pred[[i]]
    yrs <- wf_year_list[[i]]
    as.data.frame.table(m, responseName = "pred", stringsAsFactors = FALSE) %>%
      rename(id = Var1, weight = Var2) %>%
      mutate(fishery = names(wf_pred)[i],
            id      = as.integer(factor(id, levels = unique(id))),
            year    = yrs[id],
            weight  = as.numeric(weight))
  })

  # ---- Observed: long data frame (first WF fishery only for this plot) ----
  df_wf_obs <- data.frame(data$wf_obs_in[
    seq_len(data$wf_n_f[1]), , drop = FALSE])
  names(df_wf_obs) <- data$wt_mid
  df_wf_obs <- df_wf_obs %>%
    mutate(year    = wf_year_list[[1]],
          fishery = as.character(data$wf_fishery_f[1])) %>%
    mutate(id = row_number()) %>%
    pivot_longer(cols = -c(year, fishery, id),
                names_to = "weight", values_to = "obs") %>%
    mutate(weight = as.numeric(weight)) %>%
    filter(weight >= wt_mid_trim[1], weight <= wt_mid_trim[length(wt_mid_trim)])

  df_wf <- left_join(df_wf_obs, df_wf_pred, by = c("id", "fishery", "weight"))

  # ---- Plot first fishery with WF data ----
  f_plot <- as.character(data$wf_fishery_f[1])
  max_yrs_wf_plot <- 24  # show at most 24 time steps for readability
  yrs_wf_plot <- sort(unique(df_wf$year))[seq_len(min(max_yrs_wf_plot, length(unique(df_wf$year))))]
  ggplot(df_wf %>% filter(fishery == f_plot, year %in% yrs_wf_plot),
        aes(x = weight)) +
    geom_col(aes(y = obs), fill = "grey70") +
    geom_line(aes(y = pred), colour = "red3", linewidth = 0.7) +
    facet_wrap(~ year, scales = "free_y", ncol = 6) +
    labs(x = "Weight (kg)", y = "Proportion",
        title = paste0("Weight composition: Fishery ", f_plot),
        subtitle = "Grey bars = observed, red = predicted") +
    theme(strip.text = element_text(size = 7),
          axis.text  = element_text(size = 6))
}
#> Weight composition likelihood is switched off (wf_switch = 0), so no predicted weight compositions to show.
```

Inspect fitted CPUE and catch:

``` r
if(data$cpue_switch == 0){
  cat("CPUE likelihood is switched off (cpue_switch = 0), so no predicted CPUE to show.\n")
} else {
  plot(data$cpue_data$value, col = 2, pch = 16,
      xlab = "Time step", ylab = "CPUE", main = "CPUE: observed vs predicted")
  lines(exp(obj$simulate()$cpue_log_obs), lwd = 2, col = "gray70")
  lines(obj$report()$cpue_pred, lwd = 2)
}
```

![](bet_files/figure-html/opt-cpue-1.png)

``` r
sum(obj$report()$catch_pred_ysf - data$catch_obs_ysf)
#> [1] 0.9998705
plot_catch(data = data, obj = obj)
#> [1] "The maximum catch difference was: 2.20916081161704e-06"
```

![](bet_files/figure-html/opt-catch-1.png)

``` r
plot_catch(data = data, obj = obj, plot_resid = TRUE)
#> [1] "The maximum catch difference was: 2.20916081161704e-06"
```

![](bet_files/figure-html/opt-catch-2.png)

Inspect estimated spawning biomass trajectory:

``` r
plot(obj$report()$spawning_biomass_y, type = "l",
     xlab = "Time step", ylab = "Spawning biomass (mt)",
     main = "Estimated spawning biomass trajectory",
     ylim = c(0, max(obj$report()$spawning_biomass_y) * 1.2))
```

![](bet_files/figure-html/opt-spbio-1.png)

### Selectivity

Visualise the selectivity curves by fleet (at length):

``` r
par_sel_est <- obj$env$parList()$par_sel

fleet_names <- c("F01_LL.NORTH", "F02_LL.US", "F03_LL.OFFSH",
                 "F04_LL.EQUAT", "F05_LL.WEST", "F06_LL.SOUTH", "F07_LL.AUS",
                 "F08_PS.ASSOC", "F09_PS.UNASS", "F10_DOM.MISC",
                 "F11_DOM.HL", "F12_JP.PS.N", "F13_JP.PL", "F14_EQ.PL",
                 "S01_INDEX")

sel_len_df <- do.call(rbind, lapply(seq_len(data$n_fishery), function(f) {
  par_f <- par_sel_est[f, ]
  sel_l <- if (data$sel_type_f[f] == 1L) {
    sel_logistic(len_mid, par_f)
  } else {
    sel_double_normal(len_mid, par_f)
  }
  data.frame(
    fishery    = f,
    fleet_name = fleet_names[f],
    length     = len_mid,
    selectivity = as.numeric(sel_l)
  )
}))

ggplot(sel_len_df, aes(x = length, y = selectivity)) +
  geom_line() +
  facet_wrap(~fleet_name, ncol = 4) +
  labs(x = "Length (cm)", y = "Selectivity") +
  ylim(0, 1)
```

![Selectivity-at-Length by
Fleet](bet_files/figure-html/plot-selectivity-length-1.png)

Selectivity-at-Length by Fleet

Visualise the selectivity curves by fleet (at age):

``` r
rep <- obj$report()
sel_fya <- rep$sel_fya

sel_df <- expand.grid(fishery = 1:data$n_fishery, age = 1:data$n_age)
sel_df$selectivity <- sapply(1:nrow(sel_df), function(i) {
  sel_fya[sel_df$fishery[i], 1, sel_df$age[i]]
})
sel_df$real_age <- sel_df$age / 4
sel_df$fleet_name <- fleet_names[sel_df$fishery]

ggplot(sel_df, aes(x = real_age, y = selectivity)) +
  geom_line() +
  facet_wrap(~fleet_name, ncol = 4) +
  labs(x = "Age (years)", y = "Selectivity") +
  ylim(0, 1)
```

![Selectivity-at-Age by
Fleet](bet_files/figure-html/plot-selectivity-1.png)

Selectivity-at-Age by Fleet

## Simulation

The CPUE series is set up using RTMB’s `OBS` mechanism inside the model,
which allows simulation via `obj$simulate()`. For example:

``` r
if(data$cpue_switch == 0){
  cat("CPUE likelihood is switched off (cpue_switch = 0), so no simulated CPUE to show.\n")
} else {
  plot(log(data$cpue_data$value), col = 2, pch = 16,
      xlab = "Time step", ylab = "log(CPUE)", main = "Observed vs simulated CPUE")
  lines(log(obj$report()$cpue_pred), lwd = 2)
  for (i in 1:5) lines(obj$simulate()$cpue_log_obs, col = "gray70")
}
```

![](bet_files/figure-html/sim-1-1.png)

## One step ahead (OSA) residuals

OSA residuals are a replacement for Pearson residuals:

``` r
if(data$cpue_switch == 0){
  cat("CPUE likelihood is switched off (cpue_switch = 0), so no OSA residuals to show.\n")
} else {
  osa_cpue <- oneStepPredict(obj = obj, observation.name = "cpue_log_obs",
                            method = "oneStepGeneric", trace = FALSE)
  qqnorm(osa_cpue$res)
  abline(0, 1)
  plot(osa_cpue$res, xlab = "Index", ylab = "OSA residual")
  abline(h = c(-2, 0, 2), lty = c(2, 1, 2))
}
```

![](bet_files/figure-html/osa1-1.png)![](bet_files/figure-html/osa1-2.png)
