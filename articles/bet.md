# The BET model

## Introduction

The `opal` package is an R package that contains example fisheries data
from western and central pacific bigeye tuna (BET). This page provides
examples using the `opal` BET model.

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
map_sel[8, 1] <- 1
map_sel[9, 1] <- 2
map_sel[8, 3] <- 3
map_sel[9, 3] <- 4
map_sel[8, 4] <- 5
map_sel[9, 4] <- 6
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
data$lf_switch <- 1 # multinomial likelihood on flat counts (default)
# data$lf_switch <- 2 # fails while fitting, need to sort out log_lf_tau for this, use simulate to tune and find
# data$lf_switch <- 3 # fails - wants integers - think this should be an issue to RTMBdist guys

# Optionally, we could also specify random effects here (e.g. `random = "rdev_y"`), but we'll start with a simpler fixed-effects model to check everything is working first.
obj <- MakeADFun(func = cmb(opal_model, data), parameters = parameters, map = map)
unique(names(obj$par))
#> [1] "log_B0"     "log_cpue_q" "par_sel"    "rdev_y"
obj$fn()
#> [1] 134424.1
obj$gr()
#> outer mgc:  133368
#>          [,1]     [,2]    [,3]     [,4]     [,5]      [,6]      [,7]     [,8]
#> [1,] 1.283559 847.8735 -133368 34580.67 10039.04 -1615.304 -2525.923 5595.121
#>           [,9]     [,10]     [,11]     [,12]     [,13]     [,14]     [,15]
#> [1,] -14.42373 -17.79091 -11.96625 -12.08163 -23.62672 -17.25554 -17.79412
#>        [,16]     [,17]     [,18]     [,19]     [,20]     [,21]     [,22]
#> [1,] -22.425 -21.47014 -9.492821 -9.487043 -34.25532 -13.94378 -24.41817
#>          [,23]     [,24]     [,25]     [,26]     [,27]     [,28]     [,29]
#> [1,] -21.90769 -13.48098 -7.581472 -20.25591 -26.55707 -11.59815 -17.71141
#>          [,30]     [,31]     [,32]    [,33]    [,34]     [,35]     [,36]
#> [1,] -14.77482 -11.33642 -8.327971 -8.96007 -11.6725 -14.40481 -18.79485
#>          [,37]     [,38]     [,39]    [,40]     [,41]     [,42]    [,43]
#> [1,] -10.89253 -12.52936 -13.77128 -16.7568 -7.612169 -19.80656 -13.2092
#>          [,44]     [,45]     [,46]     [,47]     [,48]    [,49]     [,50]
#> [1,] -10.29417 -14.54021 -12.02666 -28.96282 -10.74705 -10.2815 -14.14097
#>          [,51]    [,52]     [,53]     [,54]     [,55]     [,56]     [,57]
#> [1,] -22.33503 -9.90664 -8.281683 -11.64374 -9.064876 -10.29405 -9.664478
#>          [,58]     [,59]     [,60]     [,61]     [,62]     [,63]     [,64]
#> [1,] -13.90944 -8.889367 -10.53526 -21.16105 -6.615882 -7.256252 -8.705208
#>          [,65]    [,66]     [,67]     [,68]     [,69]     [,70]     [,71]
#> [1,] -6.171801 -5.83139 -6.230671 -4.644741 -12.92809 -4.819266 -4.237605
#>         [,72]     [,73]     [,74]     [,75]     [,76]     [,77]     [,78]
#> [1,] -4.83079 -4.725531 -3.042983 -2.802723 -2.735913 -3.041889 -2.956688
#>          [,79]     [,80]     [,81]     [,82]    [,83]    [,84]     [,85]
#> [1,] -2.832056 -2.950518 -3.325519 -3.373152 -2.71565 -1.90421 -1.532612
#>          [,86]     [,87]      [,88]     [,89]     [,90]     [,91]     [,92]
#> [1,] -1.110487 -2.979341 -0.7495767 -3.512314 -1.775151 -2.770009 -2.894131
#>          [,93]     [,94]     [,95]     [,96]     [,97]     [,98]    [,99]
#> [1,] -3.200336 -3.530947 -3.947154 -3.963475 -3.887601 -3.864565 -3.65377
#>        [,100]    [,101]    [,102]    [,103]    [,104]    [,105]    [,106]
#> [1,] -3.86449 -3.413029 -2.858822 -2.516414 -2.421942 -1.873611 0.1267772
#>         [,107]   [,108]  [,109]   [,110]   [,111]   [,112]  [,113]   [,114]
#> [1,] 0.9110396 4.349335 3.56786 -1.37208 3.879227 1.018149 1.27915 -1.85131
#>        [,115]  [,116]   [,117]   [,118]     [,119]   [,120]   [,121]   [,122]
#> [1,] 2.258086 1.64341 3.415289 3.686082 -0.3844454 10.15892 4.477783 8.043484
#>        [,123]   [,124]   [,125]   [,126]   [,127]   [,128]   [,129]   [,130]
#> [1,] 9.868662 2.966459 4.498426 10.61197 6.183317 5.398775 20.12911 1.603732
#>         [,131]    [,132]   [,133]   [,134]   [,135]   [,136]     [,137]
#> [1,] 0.1650621 0.9320455 10.95296 14.25796 5.440251 19.92351 -0.3016632
#>        [,138]   [,139]   [,140]   [,141]   [,142]   [,143]   [,144]   [,145]
#> [1,] 1.697023 4.131396 4.966179 4.839128 4.653109 18.31732 7.044734 3.955484
#>        [,146]  [,147]   [,148]   [,149]   [,150]    [,151]    [,152]   [,153]
#> [1,] 3.024506 92.8257 35.58872 45.51684 73.66507 -72.26733 -184.5494 65.96103
#>         [,154]    [,155]   [,156]   [,157]  [,158]   [,159]    [,160]    [,161]
#> [1,] -318.1379 -179.2867 255.5759 235.1899 187.914 -42.7416 -88.10128 -17.57306
#>         [,162]    [,163]    [,164]    [,165]   [,166]   [,167]     [,168]
#> [1,] -111.0605 -133.7156 -269.3846 -109.8594 201.4434 30.68756 -0.4588547
#>        [,169]   [,170]    [,171]   [,172]   [,173]    [,174]    [,175]   [,176]
#> [1,] 161.8331 103.0581 -16.94057 29.59203 312.2648 -125.2455 -15.82691 275.7488
#>         [,177]    [,178]    [,179]    [,180]    [,181]   [,182]   [,183]
#> [1,] -85.94097 -173.5342 -610.3124 -128.6917 -513.4837 25.29142 746.1919
#>        [,184]    [,185]    [,186]    [,187]   [,188]   [,189]   [,190]
#> [1,] 283.4917 -606.6894 -1198.241 -344.7203 1597.985 453.6914 616.1131
#>         [,191]    [,192]    [,193]    [,194]   [,195]   [,196]   [,197]
#> [1,] -291.8059 -470.7017 -2306.502 -349.7535 537.7777 1572.384 2951.304
#>        [,198]    [,199]    [,200]    [,201]    [,202]    [,203]    [,204]
#> [1,] 48.57432 -500.9493 -355.4578 -420.0505 -656.7727 -953.5728 -845.6572
#>        [,205]   [,206]  [,207]   [,208]   [,209]    [,210]    [,211]    [,212]
#> [1,] 406.8608 1047.019 178.954 498.7094 155.1375 -271.9153 -213.7467 -41.94074
#>        [,213]   [,214]    [,215]    [,216]   [,217]    [,218]    [,219]
#> [1,] 344.8872 355.1801 -360.7799 -287.8721 173.4785 -129.4135 -46.37659
#>        [,220]   [,221]   [,222]    [,223]   [,224]    [,225]   [,226]   [,227]
#> [1,] 204.7333 788.5733 -8.16132 -272.3347 19.99657 -71.28219 541.5011 197.5311
#>        [,228]    [,229]    [,230]   [,231]  [,232]   [,233]   [,234]    [,235]
#> [1,] 308.8046 -553.5524 -870.1115 -241.801 351.364 629.2007 413.6063 -79.50617
#>         [,236]    [,237]    [,238]    [,239]    [,240]    [,241]   [,242]
#> [1,] -594.2893 -298.4981 -1081.499 -415.2285 -130.2232 -60.01887 3185.805
#>        [,243]   [,244]   [,245]   [,246]    [,247]    [,248]   [,249]   [,250]
#> [1,] 408.0679 1235.054 142.8773 29.82557 -1253.426 -671.8458 812.6589 69.88638
#>         [,251]    [,252]   [,253]    [,254]    [,255]   [,256]    [,257]
#> [1,] -600.1491 -141.6265 373.5507 -780.2418 -706.1567 233.5268 -158.4664
#>         [,258]    [,259]  [,260]    [,261]    [,262]   [,263]  [,264]   [,265]
#> [1,] -406.3222 -515.1598 161.159 -725.4307 -536.3243 869.2207 386.896 93.48939
#>        [,266]   [,267]   [,268]    [,269]   [,270]   [,271]   [,272]   [,273]
#> [1,] 692.0882 595.3564 1460.613 -1.630662 -933.118 13.88269 253.6477 -613.092
#>         [,274]    [,275]    [,276]
#> [1,] -1011.023 -42.43524 0.6256006
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
#> [1] "The maximum catch difference was: 1.10399014374707e-08"
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
```

Compare initial and estimated parameter values:

``` r
get_par_table(obj, parameters, map, lower = Lwr, upper = Upr, grad_tol = 1e-2, digits = 2L)
#> outer mgc:  0.005814213
#>          par   init     est  lwr  upr       gr gr_chk bd_chk
#> 1     log_B0 20.000 14.0000  0.0 22.0 -4.3e-05     OK     OK
#> 2 log_cpue_q  0.000  0.0079 -2.3  2.3 -5.4e-08     OK     OK
#> 3   par_sel8 -0.980 -0.9700 -Inf  Inf -3.4e-03     OK     OK
#> 4   par_sel9 -0.780 -1.0000 -Inf  Inf -5.8e-03     OK     OK
#> 5  par_sel38  0.044 -9.4000 -Inf  Inf  8.7e-06     OK     OK
#> 6  par_sel39  1.600 -9.3000 -Inf  Inf  1.6e-05     OK     OK
#> 7  par_sel53  0.450 -2.3000 -Inf  Inf -4.1e-06     OK     OK
#> 8  par_sel54  3.000  1.7000 -Inf  Inf -8.4e-06     OK     OK
```

### Diagnostics

Check that all parameters are estimable using the `check_estimability`
function:

``` r
check_estimability(obj = obj)
#> outer mgc:  0.005814213 
#> outer mgc:  347.42 
#> outer mgc:  351.1373 
#> outer mgc:  7.425739 
#> outer mgc:  7.425739 
#> outer mgc:  3890.009 
#> outer mgc:  8195.263 
#> outer mgc:  465.0801 
#> outer mgc:  992.1283 
#> outer mgc:  0.005814301 
#> outer mgc:  0.005814126 
#> outer mgc:  0.005918232 
#> outer mgc:  0.005711919 
#> outer mgc:  58.15584 
#> outer mgc:  58.12709 
#> outer mgc:  21.19431 
#> outer mgc:  21.20637 
#> outer mgc:  0.3042916 
#> outer mgc:  0.3039646 
#> outer mgc:  0.01688228 
#> outer mgc:  0.01686521 
#> outer mgc:  0.01673787 
#> outer mgc:  0.01672093 
#> outer mgc:  0.01686755 
#> outer mgc:  0.01685047 
#> outer mgc:  0.01659516 
#> outer mgc:  0.01657793 
#> outer mgc:  0.01697671 
#> outer mgc:  0.01695946 
#> outer mgc:  0.0168465 
#> outer mgc:  0.01682941 
#> outer mgc:  0.01652358 
#> outer mgc:  0.01650683 
#> outer mgc:  0.01566825 
#> outer mgc:  0.01565201 
#> outer mgc:  0.01544898 
#> outer mgc:  0.0154333 
#> outer mgc:  0.01504045 
#> outer mgc:  0.01502522 
#> outer mgc:  0.271319 
#> outer mgc:  0.2710141 
#> outer mgc:  0.01437019 
#> outer mgc:  0.01435543 
#> outer mgc:  0.01477366 
#> outer mgc:  0.0147586 
#> outer mgc:  0.01492105 
#> outer mgc:  0.01490586 
#> outer mgc:  0.01501281 
#> outer mgc:  0.0149975 
#> outer mgc:  0.0141765 
#> outer mgc:  0.01416136 
#> outer mgc:  0.01457218 
#> outer mgc:  0.01455717 
#> outer mgc:  0.01406493 
#> outer mgc:  0.01405044 
#> outer mgc:  0.01329359 
#> outer mgc:  0.01327984 
#> outer mgc:  0.01163088 
#> outer mgc:  0.01161817 
#> outer mgc:  0.01141218 
#> outer mgc:  0.01140019 
#> outer mgc:  0.01065471 
#> outer mgc:  0.01064348 
#> outer mgc:  0.009949391 
#> outer mgc:  0.009938813 
#> outer mgc:  0.008817761 
#> outer mgc:  0.008807948 
#> outer mgc:  0.008512605 
#> outer mgc:  0.008503329 
#> outer mgc:  0.008523387 
#> outer mgc:  0.008428557 
#> outer mgc:  0.008684793 
#> outer mgc:  0.008589807 
#> outer mgc:  0.008689872 
#> outer mgc:  0.008594781 
#> outer mgc:  0.0089683 
#> outer mgc:  0.008873037 
#> outer mgc:  0.009113287 
#> outer mgc:  0.009017896 
#> outer mgc:  0.009370842 
#> outer mgc:  0.009275203 
#> outer mgc:  0.009725815 
#> outer mgc:  0.009629775 
#> outer mgc:  0.01031767 
#> outer mgc:  0.01022108 
#> outer mgc:  0.01087413 
#> outer mgc:  0.01077705 
#> outer mgc:  0.01112326 
#> outer mgc:  0.01102594 
#> outer mgc:  0.01089429 
#> outer mgc:  0.01079712 
#> outer mgc:  0.01052345 
#> outer mgc:  0.01042666 
#> outer mgc:  0.0101078 
#> outer mgc:  0.01001146 
#> outer mgc:  0.009755902 
#> outer mgc:  0.009659914 
#> outer mgc:  0.009568277 
#> outer mgc:  0.009472438 
#> outer mgc:  0.009589216 
#> outer mgc:  0.009493369 
#> outer mgc:  0.009743978 
#> outer mgc:  0.009648008 
#> outer mgc:  0.01009576 
#> outer mgc:  0.00999945 
#> outer mgc:  0.01063821 
#> outer mgc:  0.01054131 
#> outer mgc:  0.01143027 
#> outer mgc:  0.01133263 
#> outer mgc:  0.01188737 
#> outer mgc:  0.01178934 
#> outer mgc:  0.0118354 
#> outer mgc:  0.01173742 
#> outer mgc:  0.01159499 
#> outer mgc:  0.0114972 
#> outer mgc:  0.01123659 
#> outer mgc:  0.01113916 
#> outer mgc:  0.01066512 
#> outer mgc:  0.01056827 
#> outer mgc:  0.01016346 
#> outer mgc:  0.01006709 
#> outer mgc:  0.009771673 
#> outer mgc:  0.009675632 
#> outer mgc:  0.009451658 
#> outer mgc:  0.009355942 
#> outer mgc:  0.009101834 
#> outer mgc:  0.009006492 
#> outer mgc:  0.008855073 
#> outer mgc:  0.00875997 
#> outer mgc:  0.008804128 
#> outer mgc:  0.008709024 
#> outer mgc:  0.009033596 
#> outer mgc:  0.008938283 
#> outer mgc:  0.009428661 
#> outer mgc:  0.009332995 
#> outer mgc:  0.009937349 
#> outer mgc:  0.009841185 
#> outer mgc:  0.01016883 
#> outer mgc:  0.01007234 
#> outer mgc:  0.009903808 
#> outer mgc:  0.009807598 
#> outer mgc:  0.009019507 
#> outer mgc:  0.008924225 
#> outer mgc:  0.008224952 
#> outer mgc:  0.007970253 
#> outer mgc:  0.007791051 
#> outer mgc:  0.007217033 
#> outer mgc:  0.007592337 
#> outer mgc:  0.006852852 
#> outer mgc:  0.007567097 
#> outer mgc:  0.006773861 
#> outer mgc:  0.007685113 
#> outer mgc:  0.006920103 
#> outer mgc:  0.00800093 
#> outer mgc:  0.007382354 
#> outer mgc:  0.008560641 
#> outer mgc:  0.008226598 
#> outer mgc:  0.009261893 
#> outer mgc:  0.009166291 
#> outer mgc:  0.01004301 
#> outer mgc:  0.009946636 
#> outer mgc:  0.0105068 
#> outer mgc:  0.01040991 
#> outer mgc:  0.01018391 
#> outer mgc:  0.01008735 
#> outer mgc:  0.009198618 
#> outer mgc:  0.009008019 
#> outer mgc:  0.00851311 
#> outer mgc:  0.007917411 
#> outer mgc:  0.007978154 
#> outer mgc:  0.007070931 
#> outer mgc:  0.007613693 
#> outer mgc:  0.006555783 
#> outer mgc:  0.007372657 
#> outer mgc:  0.006313704 
#> outer mgc:  0.00726984 
#> outer mgc:  0.006255238 
#> outer mgc:  0.007331487 
#> outer mgc:  0.006450664 
#> outer mgc:  0.00754939 
#> outer mgc:  0.006863183 
#> outer mgc:  0.007867802 
#> outer mgc:  0.007360846 
#> outer mgc:  0.008280704 
#> outer mgc:  0.007927066 
#> outer mgc:  0.00871985 
#> outer mgc:  0.008455316 
#> outer mgc:  0.009132536 
#> outer mgc:  0.008871843 
#> outer mgc:  0.009370432 
#> outer mgc:  0.009005246 
#> outer mgc:  0.009414892 
#> outer mgc:  0.008976727 
#> outer mgc:  0.009287637 
#> outer mgc:  0.00875525 
#> outer mgc:  0.009012937 
#> outer mgc:  0.008310484 
#> outer mgc:  0.008713204 
#> outer mgc:  0.007830172 
#> outer mgc:  0.008577917 
#> outer mgc:  0.007592788 
#> outer mgc:  0.008668877 
#> outer mgc:  0.007684435 
#> outer mgc:  0.008944078 
#> outer mgc:  0.008044959 
#> outer mgc:  0.009152839 
#> outer mgc:  0.008619397 
#> outer mgc:  0.009027656 
#> outer mgc:  0.008895536 
#> outer mgc:  0.008561575 
#> outer mgc:  0.008554818 
#> outer mgc:  0.008256897 
#> outer mgc:  0.008250267 
#> outer mgc:  0.008033942 
#> outer mgc:  0.008027415 
#> outer mgc:  0.00805381 
#> outer mgc:  0.008047196 
#> outer mgc:  0.008202076 
#> outer mgc:  0.008195175 
#> outer mgc:  0.008433295 
#> outer mgc:  0.008426284 
#> outer mgc:  0.008300739 
#> outer mgc:  0.008293856 
#> outer mgc:  0.007847038 
#> outer mgc:  0.00784052 
#> outer mgc:  0.007042925 
#> outer mgc:  0.007036829 
#> outer mgc:  0.006423638 
#> outer mgc:  0.00641819 
#> outer mgc:  0.005985482 
#> outer mgc:  0.00586921 
#> outer mgc:  0.005969096 
#> outer mgc:  0.005659481 
#> outer mgc:  0.005956586 
#> outer mgc:  0.005671981 
#> outer mgc:  0.005950864 
#> outer mgc:  0.005677695 
#> outer mgc:  0.005946057 
#> outer mgc:  0.005785156 
#> outer mgc:  0.006121307 
#> outer mgc:  0.006116402 
#> outer mgc:  0.006447632 
#> outer mgc:  0.00644256 
#> outer mgc:  0.006666683 
#> outer mgc:  0.006661648 
#> outer mgc:  0.006547545 
#> outer mgc:  0.00654277 
#> outer mgc:  0.006238794 
#> outer mgc:  0.006234416 
#> outer mgc:  0.005833388 
#> outer mgc:  0.005828263 
#> outer mgc:  0.005799077 
#> outer mgc:  0.005829345 
#> outer mgc:  0.005765781 
#> outer mgc:  0.005862611 
#> outer mgc:  0.005733158 
#> outer mgc:  0.005895205 
#> outer mgc:  0.005701619 
#> outer mgc:  0.005926713 
#> outer mgc:  0.005666037 
#> outer mgc:  0.005962263 
#> outer mgc:  0.005628721 
#> outer mgc:  0.005999547 
#> outer mgc:  0.005588318 
#> outer mgc:  0.006039917 
#> outer mgc:  0.005548072 
#> outer mgc:  0.006080122 
#> outer mgc:  0.005495784 
#> outer mgc:  0.006132368 
#> outer mgc:  0.005435181 
#> outer mgc:  0.006192926 
#> outer mgc:  0.005371342 
#> outer mgc:  0.006256709 
#> outer mgc:  0.005302249 
#> outer mgc:  0.006325736 
#> outer mgc:  0.005212001 
#> outer mgc:  0.006415908 
#> outer mgc:  0.005093913 
#> outer mgc:  0.006533899 
#> outer mgc:  0.005090693 
#> outer mgc:  0.006691797 
#> outer mgc:  0.005328581 
#> outer mgc:  0.006911909 
#> outer mgc:  0.005354578 
#> outer mgc:  0.007265944 
#> outer mgc:  0.005072143 
#> outer mgc:  0.00785322 
#> outer mgc:  0.004544398 
#> outer mgc:  0.008800212 
#> outer mgc:  0.004274964 
#> outer mgc:  0.01102501 
#> outer mgc:  0.00837789 
#> outer mgc:  0.01517419 
#> outer mgc:  0.0148112 
#> outer mgc:  0.02160205 
#> outer mgc:  0.02432195 
#> outer mgc:  0.03110304 
#> outer mgc:  0.03591293 
#> outer mgc:  0.04267642 
#> outer mgc:  0.05634112 
#> outer mgc:  0.06308545 
#> outer mgc:  0.7418375 
#> outer mgc:  0.7481334 
#> outer mgc:  2.155934 
#> outer mgc:  2.16165 
#> outer mgc:  0.07896789 
#> outer mgc:  0.08569386 
#> outer mgc:  0.5787951 
#> outer mgc:  0.5846762 
#> outer mgc:  0.8247599 
#> outer mgc:  0.831068 
#> outer mgc:  0.1043998 
#> outer mgc:  0.1111215 
#> outer mgc:  0.1210241 
#> outer mgc:  0.1141382 
#> outer mgc:  1.481842 
#> outer mgc:  1.47538 
#> outer mgc:  0.8193408 
#> outer mgc:  0.8121152 
#> outer mgc:  0.3935518 
#> outer mgc:  0.4005771 
#> outer mgc:  0.06584296 
#> outer mgc:  0.07258385 
#> outer mgc:  1.111809 
#> outer mgc:  1.104889 
#> outer mgc:  0.2947721 
#> outer mgc:  0.3012738 
#> outer mgc:  6.879335 
#> outer mgc:  6.885098 
#> outer mgc:  0.1283628 
#> outer mgc:  0.13504 
#> outer mgc:  1.184257 
#> outer mgc:  1.177687 
#> outer mgc:  1.515492 
#> outer mgc:  1.507782 
#> outer mgc:  0.2788924 
#> outer mgc:  0.2855042 
#> outer mgc:  0.6796601 
#> outer mgc:  0.6860741 
#> outer mgc:  0.5504502 
#> outer mgc:  0.5568282 
#> outer mgc:  0.7450391 
#> outer mgc:  0.7517251 
#> outer mgc:  0.3958923 
#> outer mgc:  0.4022327 
#> outer mgc:  0.06147695 
#> outer mgc:  0.06825507 
#> outer mgc:  0.4966519 
#> outer mgc:  0.5041882 
#> outer mgc:  0.7254516 
#> outer mgc:  0.7181549 
#> outer mgc:  0.07289171 
#> outer mgc:  0.06592496 
#> outer mgc:  0.6026561 
#> outer mgc:  0.6086213 
#> outer mgc:  0.5543888 
#> outer mgc:  0.5607264 
#> outer mgc:  1.300134 
#> outer mgc:  1.305874 
#> outer mgc:  0.1926756 
#> outer mgc:  0.1993436 
#> outer mgc:  5.574363 
#> outer mgc:  5.57848 
#> outer mgc:  0.0709951 
#> outer mgc:  0.07772749 
#> outer mgc:  0.2460567 
#> outer mgc:  0.2456817 
#> outer mgc:  1.117081 
#> outer mgc:  1.122467 
#> outer mgc:  0.5795667 
#> outer mgc:  0.5788785 
#> outer mgc:  1.205114 
#> outer mgc:  1.204542 
#> outer mgc:  0.4461172 
#> outer mgc:  0.452669 
#> outer mgc:  0.9735659 
#> outer mgc:  0.9809175 
#> outer mgc:  0.01533163 
#> outer mgc:  0.01538453 
#> outer mgc:  2.779095 
#> outer mgc:  2.770233 
#> outer mgc:  3.013237 
#> outer mgc:  3.005858 
#> outer mgc:  7.590464 
#> outer mgc:  7.592583 
#> outer mgc:  10.9211 
#> outer mgc:  10.92415 
#> outer mgc:  3.537771 
#> outer mgc:  3.54227 
#> outer mgc:  0.304755 
#> outer mgc:  0.311397 
#> outer mgc:  0.1706194 
#> outer mgc:  0.1635631 
#> outer mgc:  0.8772995 
#> outer mgc:  0.8700635 
#> outer mgc:  1.522929 
#> outer mgc:  1.521785 
#> outer mgc:  2.873678 
#> outer mgc:  2.881836 
#> outer mgc:  1.06335 
#> outer mgc:  1.06258 
#> outer mgc:  1.233552 
#> outer mgc:  1.225529 
#> outer mgc:  8.522805 
#> outer mgc:  8.524982 
#> outer mgc:  6.760632 
#> outer mgc:  6.763145 
#> outer mgc:  1.842182 
#> outer mgc:  1.847702 
#> outer mgc:  4.006611 
#> outer mgc:  3.997721 
#> outer mgc:  4.248198 
#> outer mgc:  4.23627 
#> outer mgc:  1.268859 
#> outer mgc:  1.261518 
#> outer mgc:  0.1612047 
#> outer mgc:  0.1680634 
#> outer mgc:  0.2803039 
#> outer mgc:  0.2869764 
#> outer mgc:  4.211275 
#> outer mgc:  4.216344 
#> outer mgc:  3.41939 
#> outer mgc:  3.424782 
#> outer mgc:  3.232833 
#> outer mgc:  3.237916 
#> outer mgc:  1.486452 
#> outer mgc:  1.491889 
#> outer mgc:  0.03107697 
#> outer mgc:  0.03785861 
#> outer mgc:  3.450172 
#> outer mgc:  3.45573 
#> outer mgc:  2.65533 
#> outer mgc:  2.662088 
#> outer mgc:  1.17144 
#> outer mgc:  1.177339 
#> outer mgc:  3.432255 
#> outer mgc:  3.437952 
#> outer mgc:  3.82763 
#> outer mgc:  3.833155 
#> outer mgc:  0.7330496 
#> outer mgc:  0.7393711 
#> outer mgc:  0.2860441 
#> outer mgc:  0.2790369 
#> outer mgc:  0.534149 
#> outer mgc:  0.533641 
#> outer mgc:  5.023069 
#> outer mgc:  5.027728 
#> outer mgc:  3.23317 
#> outer mgc:  3.239002 
#> outer mgc:  5.141974 
#> outer mgc:  5.148815 
#> outer mgc:  0.6807245 
#> outer mgc:  0.6868668 
#> outer mgc:  3.708177 
#> outer mgc:  3.712172 
#> outer mgc:  1.022393 
#> outer mgc:  1.028466 
#> outer mgc:  5.778456 
#> outer mgc:  5.784651 
#> outer mgc:  4.67385 
#> outer mgc:  4.679382 
#> outer mgc:  2.909998 
#> outer mgc:  2.915709 
#> outer mgc:  0.3652791 
#> outer mgc:  0.3647707 
#> outer mgc:  0.8994371 
#> outer mgc:  0.8927196 
#> outer mgc:  1.174924 
#> outer mgc:  1.182068 
#> outer mgc:  2.058319 
#> outer mgc:  2.064897 
#> outer mgc:  5.807428 
#> outer mgc:  5.810696 
#> outer mgc:  3.950383 
#> outer mgc:  3.954096 
#> outer mgc:  12.19628 
#> outer mgc:  12.19668 
#> outer mgc:  7.180289 
#> outer mgc:  7.182716 
#> outer mgc:  1.688183 
#> outer mgc:  1.693474 
#> outer mgc:  10.89622 
#> outer mgc:  10.89814 
#> outer mgc:  6.896194 
#> outer mgc:  6.903373 
#> outer mgc:  6.130411 
#> outer mgc:  6.138583 
#> outer mgc:  4.239735 
#> outer mgc:  4.246409 
#> outer mgc:  8.230351 
#> outer mgc:  8.234589 
#> outer mgc:  4.239054 
#> outer mgc:  4.252703 
#> outer mgc:  5.644993 
#> outer mgc:  5.658023 
#> outer mgc:  5.7308 
#> outer mgc:  5.738049 
#> outer mgc:  5.531328 
#> outer mgc:  5.535201 
#> outer mgc:  1.420927 
#> outer mgc:  1.419385 
#> outer mgc:  1.021527 
#> outer mgc:  1.020361 
#> outer mgc:  3.468327 
#> outer mgc:  3.47458 
#> outer mgc:  1.757969 
#> outer mgc:  1.763484 
#> outer mgc:  1.904011 
#> outer mgc:  1.903157 
#> outer mgc:  1.730432 
#> outer mgc:  1.736458 
#> outer mgc:  1.752873 
#> outer mgc:  1.759543 
#> outer mgc:  2.815207 
#> outer mgc:  2.819144 
#> outer mgc:  0.8573765 
#> outer mgc:  0.8683017 
#> outer mgc:  3.501611 
#> outer mgc:  3.505564 
#> outer mgc:  2.522105 
#> outer mgc:  2.527065 
#> outer mgc:  9.112851 
#> outer mgc:  9.114084 
#> outer mgc:  6.64043 
#> outer mgc:  6.644873 
#> outer mgc:  12.88707 
#> outer mgc:  12.89751 
#> outer mgc:  5.109293 
#> outer mgc:  5.11293 
#> outer mgc:  16.9188 
#> outer mgc:  16.93082 
#> outer mgc:  10.50079 
#> outer mgc:  10.5095 
#> outer mgc:  5.66418 
#> outer mgc:  5.670516 
#> outer mgc:  10.29922 
#> outer mgc:  10.30999 
#> outer mgc:  20.4458 
#> outer mgc:  20.47105 
#> outer mgc:  32.5316 
#> outer mgc:  32.64331 
#> outer mgc:  24.96667 
#> outer mgc:  25.04121 
#> outer mgc:  3.503004 
#> outer mgc:  3.507844 
#> outer mgc:  1.510851 
#> outer mgc:  1.522215 
#> outer mgc:  13.79604 
#> outer mgc:  13.78759 
#> outer mgc:  0.08780279 
#> outer mgc:  0.08091153 
#> outer mgc:  0.005814213 
#> outer mgc:  0.005814213
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
#> outer mgc:  0.005814213 
#> outer mgc:  347.42 
#> outer mgc:  351.1373 
#> outer mgc:  7.425739 
#> outer mgc:  7.425739 
#> outer mgc:  3890.009 
#> outer mgc:  8195.263 
#> outer mgc:  465.0801 
#> outer mgc:  992.1283 
#> outer mgc:  0.005814301 
#> outer mgc:  0.005814126 
#> outer mgc:  0.005918232 
#> outer mgc:  0.005711919 
#> outer mgc:  58.15584 
#> outer mgc:  58.12709 
#> outer mgc:  21.19431 
#> outer mgc:  21.20637 
#> outer mgc:  0.3042916 
#> outer mgc:  0.3039646 
#> outer mgc:  0.01688228 
#> outer mgc:  0.01686521 
#> outer mgc:  0.01673787 
#> outer mgc:  0.01672093 
#> outer mgc:  0.01686755 
#> outer mgc:  0.01685047 
#> outer mgc:  0.01659516 
#> outer mgc:  0.01657793 
#> outer mgc:  0.01697671 
#> outer mgc:  0.01695946 
#> outer mgc:  0.0168465 
#> outer mgc:  0.01682941 
#> outer mgc:  0.01652358 
#> outer mgc:  0.01650683 
#> outer mgc:  0.01566825 
#> outer mgc:  0.01565201 
#> outer mgc:  0.01544898 
#> outer mgc:  0.0154333 
#> outer mgc:  0.01504045 
#> outer mgc:  0.01502522 
#> outer mgc:  0.271319 
#> outer mgc:  0.2710141 
#> outer mgc:  0.01437019 
#> outer mgc:  0.01435543 
#> outer mgc:  0.01477366 
#> outer mgc:  0.0147586 
#> outer mgc:  0.01492105 
#> outer mgc:  0.01490586 
#> outer mgc:  0.01501281 
#> outer mgc:  0.0149975 
#> outer mgc:  0.0141765 
#> outer mgc:  0.01416136 
#> outer mgc:  0.01457218 
#> outer mgc:  0.01455717 
#> outer mgc:  0.01406493 
#> outer mgc:  0.01405044 
#> outer mgc:  0.01329359 
#> outer mgc:  0.01327984 
#> outer mgc:  0.01163088 
#> outer mgc:  0.01161817 
#> outer mgc:  0.01141218 
#> outer mgc:  0.01140019 
#> outer mgc:  0.01065471 
#> outer mgc:  0.01064348 
#> outer mgc:  0.009949391 
#> outer mgc:  0.009938813 
#> outer mgc:  0.008817761 
#> outer mgc:  0.008807948 
#> outer mgc:  0.008512605 
#> outer mgc:  0.008503329 
#> outer mgc:  0.008523387 
#> outer mgc:  0.008428557 
#> outer mgc:  0.008684793 
#> outer mgc:  0.008589807 
#> outer mgc:  0.008689872 
#> outer mgc:  0.008594781 
#> outer mgc:  0.0089683 
#> outer mgc:  0.008873037 
#> outer mgc:  0.009113287 
#> outer mgc:  0.009017896 
#> outer mgc:  0.009370842 
#> outer mgc:  0.009275203 
#> outer mgc:  0.009725815 
#> outer mgc:  0.009629775 
#> outer mgc:  0.01031767 
#> outer mgc:  0.01022108 
#> outer mgc:  0.01087413 
#> outer mgc:  0.01077705 
#> outer mgc:  0.01112326 
#> outer mgc:  0.01102594 
#> outer mgc:  0.01089429 
#> outer mgc:  0.01079712 
#> outer mgc:  0.01052345 
#> outer mgc:  0.01042666 
#> outer mgc:  0.0101078 
#> outer mgc:  0.01001146 
#> outer mgc:  0.009755902 
#> outer mgc:  0.009659914 
#> outer mgc:  0.009568277 
#> outer mgc:  0.009472438 
#> outer mgc:  0.009589216 
#> outer mgc:  0.009493369 
#> outer mgc:  0.009743978 
#> outer mgc:  0.009648008 
#> outer mgc:  0.01009576 
#> outer mgc:  0.00999945 
#> outer mgc:  0.01063821 
#> outer mgc:  0.01054131 
#> outer mgc:  0.01143027 
#> outer mgc:  0.01133263 
#> outer mgc:  0.01188737 
#> outer mgc:  0.01178934 
#> outer mgc:  0.0118354 
#> outer mgc:  0.01173742 
#> outer mgc:  0.01159499 
#> outer mgc:  0.0114972 
#> outer mgc:  0.01123659 
#> outer mgc:  0.01113916 
#> outer mgc:  0.01066512 
#> outer mgc:  0.01056827 
#> outer mgc:  0.01016346 
#> outer mgc:  0.01006709 
#> outer mgc:  0.009771673 
#> outer mgc:  0.009675632 
#> outer mgc:  0.009451658 
#> outer mgc:  0.009355942 
#> outer mgc:  0.009101834 
#> outer mgc:  0.009006492 
#> outer mgc:  0.008855073 
#> outer mgc:  0.00875997 
#> outer mgc:  0.008804128 
#> outer mgc:  0.008709024 
#> outer mgc:  0.009033596 
#> outer mgc:  0.008938283 
#> outer mgc:  0.009428661 
#> outer mgc:  0.009332995 
#> outer mgc:  0.009937349 
#> outer mgc:  0.009841185 
#> outer mgc:  0.01016883 
#> outer mgc:  0.01007234 
#> outer mgc:  0.009903808 
#> outer mgc:  0.009807598 
#> outer mgc:  0.009019507 
#> outer mgc:  0.008924225 
#> outer mgc:  0.008224952 
#> outer mgc:  0.007970253 
#> outer mgc:  0.007791051 
#> outer mgc:  0.007217033 
#> outer mgc:  0.007592337 
#> outer mgc:  0.006852852 
#> outer mgc:  0.007567097 
#> outer mgc:  0.006773861 
#> outer mgc:  0.007685113 
#> outer mgc:  0.006920103 
#> outer mgc:  0.00800093 
#> outer mgc:  0.007382354 
#> outer mgc:  0.008560641 
#> outer mgc:  0.008226598 
#> outer mgc:  0.009261893 
#> outer mgc:  0.009166291 
#> outer mgc:  0.01004301 
#> outer mgc:  0.009946636 
#> outer mgc:  0.0105068 
#> outer mgc:  0.01040991 
#> outer mgc:  0.01018391 
#> outer mgc:  0.01008735 
#> outer mgc:  0.009198618 
#> outer mgc:  0.009008019 
#> outer mgc:  0.00851311 
#> outer mgc:  0.007917411 
#> outer mgc:  0.007978154 
#> outer mgc:  0.007070931 
#> outer mgc:  0.007613693 
#> outer mgc:  0.006555783 
#> outer mgc:  0.007372657 
#> outer mgc:  0.006313704 
#> outer mgc:  0.00726984 
#> outer mgc:  0.006255238 
#> outer mgc:  0.007331487 
#> outer mgc:  0.006450664 
#> outer mgc:  0.00754939 
#> outer mgc:  0.006863183 
#> outer mgc:  0.007867802 
#> outer mgc:  0.007360846 
#> outer mgc:  0.008280704 
#> outer mgc:  0.007927066 
#> outer mgc:  0.00871985 
#> outer mgc:  0.008455316 
#> outer mgc:  0.009132536 
#> outer mgc:  0.008871843 
#> outer mgc:  0.009370432 
#> outer mgc:  0.009005246 
#> outer mgc:  0.009414892 
#> outer mgc:  0.008976727 
#> outer mgc:  0.009287637 
#> outer mgc:  0.00875525 
#> outer mgc:  0.009012937 
#> outer mgc:  0.008310484 
#> outer mgc:  0.008713204 
#> outer mgc:  0.007830172 
#> outer mgc:  0.008577917 
#> outer mgc:  0.007592788 
#> outer mgc:  0.008668877 
#> outer mgc:  0.007684435 
#> outer mgc:  0.008944078 
#> outer mgc:  0.008044959 
#> outer mgc:  0.009152839 
#> outer mgc:  0.008619397 
#> outer mgc:  0.009027656 
#> outer mgc:  0.008895536 
#> outer mgc:  0.008561575 
#> outer mgc:  0.008554818 
#> outer mgc:  0.008256897 
#> outer mgc:  0.008250267 
#> outer mgc:  0.008033942 
#> outer mgc:  0.008027415 
#> outer mgc:  0.00805381 
#> outer mgc:  0.008047196 
#> outer mgc:  0.008202076 
#> outer mgc:  0.008195175 
#> outer mgc:  0.008433295 
#> outer mgc:  0.008426284 
#> outer mgc:  0.008300739 
#> outer mgc:  0.008293856 
#> outer mgc:  0.007847038 
#> outer mgc:  0.00784052 
#> outer mgc:  0.007042925 
#> outer mgc:  0.007036829 
#> outer mgc:  0.006423638 
#> outer mgc:  0.00641819 
#> outer mgc:  0.005985482 
#> outer mgc:  0.00586921 
#> outer mgc:  0.005969096 
#> outer mgc:  0.005659481 
#> outer mgc:  0.005956586 
#> outer mgc:  0.005671981 
#> outer mgc:  0.005950864 
#> outer mgc:  0.005677695 
#> outer mgc:  0.005946057 
#> outer mgc:  0.005785156 
#> outer mgc:  0.006121307 
#> outer mgc:  0.006116402 
#> outer mgc:  0.006447632 
#> outer mgc:  0.00644256 
#> outer mgc:  0.006666683 
#> outer mgc:  0.006661648 
#> outer mgc:  0.006547545 
#> outer mgc:  0.00654277 
#> outer mgc:  0.006238794 
#> outer mgc:  0.006234416 
#> outer mgc:  0.005833388 
#> outer mgc:  0.005828263 
#> outer mgc:  0.005799077 
#> outer mgc:  0.005829345 
#> outer mgc:  0.005765781 
#> outer mgc:  0.005862611 
#> outer mgc:  0.005733158 
#> outer mgc:  0.005895205 
#> outer mgc:  0.005701619 
#> outer mgc:  0.005926713 
#> outer mgc:  0.005666037 
#> outer mgc:  0.005962263 
#> outer mgc:  0.005628721 
#> outer mgc:  0.005999547 
#> outer mgc:  0.005588318 
#> outer mgc:  0.006039917 
#> outer mgc:  0.005548072 
#> outer mgc:  0.006080122 
#> outer mgc:  0.005495784 
#> outer mgc:  0.006132368 
#> outer mgc:  0.005435181 
#> outer mgc:  0.006192926 
#> outer mgc:  0.005371342 
#> outer mgc:  0.006256709 
#> outer mgc:  0.005302249 
#> outer mgc:  0.006325736 
#> outer mgc:  0.005212001 
#> outer mgc:  0.006415908 
#> outer mgc:  0.005093913 
#> outer mgc:  0.006533899 
#> outer mgc:  0.005090693 
#> outer mgc:  0.006691797 
#> outer mgc:  0.005328581 
#> outer mgc:  0.006911909 
#> outer mgc:  0.005354578 
#> outer mgc:  0.007265944 
#> outer mgc:  0.005072143 
#> outer mgc:  0.00785322 
#> outer mgc:  0.004544398 
#> outer mgc:  0.008800212 
#> outer mgc:  0.004274964 
#> outer mgc:  0.01102501 
#> outer mgc:  0.00837789 
#> outer mgc:  0.01517419 
#> outer mgc:  0.0148112 
#> outer mgc:  0.02160205 
#> outer mgc:  0.02432195 
#> outer mgc:  0.03110304 
#> outer mgc:  0.03591293 
#> outer mgc:  0.04267642 
#> outer mgc:  0.05634112 
#> outer mgc:  0.06308545 
#> outer mgc:  0.7418375 
#> outer mgc:  0.7481334 
#> outer mgc:  2.155934 
#> outer mgc:  2.16165 
#> outer mgc:  0.07896789 
#> outer mgc:  0.08569386 
#> outer mgc:  0.5787951 
#> outer mgc:  0.5846762 
#> outer mgc:  0.8247599 
#> outer mgc:  0.831068 
#> outer mgc:  0.1043998 
#> outer mgc:  0.1111215 
#> outer mgc:  0.1210241 
#> outer mgc:  0.1141382 
#> outer mgc:  1.481842 
#> outer mgc:  1.47538 
#> outer mgc:  0.8193408 
#> outer mgc:  0.8121152 
#> outer mgc:  0.3935518 
#> outer mgc:  0.4005771 
#> outer mgc:  0.06584296 
#> outer mgc:  0.07258385 
#> outer mgc:  1.111809 
#> outer mgc:  1.104889 
#> outer mgc:  0.2947721 
#> outer mgc:  0.3012738 
#> outer mgc:  6.879335 
#> outer mgc:  6.885098 
#> outer mgc:  0.1283628 
#> outer mgc:  0.13504 
#> outer mgc:  1.184257 
#> outer mgc:  1.177687 
#> outer mgc:  1.515492 
#> outer mgc:  1.507782 
#> outer mgc:  0.2788924 
#> outer mgc:  0.2855042 
#> outer mgc:  0.6796601 
#> outer mgc:  0.6860741 
#> outer mgc:  0.5504502 
#> outer mgc:  0.5568282 
#> outer mgc:  0.7450391 
#> outer mgc:  0.7517251 
#> outer mgc:  0.3958923 
#> outer mgc:  0.4022327 
#> outer mgc:  0.06147695 
#> outer mgc:  0.06825507 
#> outer mgc:  0.4966519 
#> outer mgc:  0.5041882 
#> outer mgc:  0.7254516 
#> outer mgc:  0.7181549 
#> outer mgc:  0.07289171 
#> outer mgc:  0.06592496 
#> outer mgc:  0.6026561 
#> outer mgc:  0.6086213 
#> outer mgc:  0.5543888 
#> outer mgc:  0.5607264 
#> outer mgc:  1.300134 
#> outer mgc:  1.305874 
#> outer mgc:  0.1926756 
#> outer mgc:  0.1993436 
#> outer mgc:  5.574363 
#> outer mgc:  5.57848 
#> outer mgc:  0.0709951 
#> outer mgc:  0.07772749 
#> outer mgc:  0.2460567 
#> outer mgc:  0.2456817 
#> outer mgc:  1.117081 
#> outer mgc:  1.122467 
#> outer mgc:  0.5795667 
#> outer mgc:  0.5788785 
#> outer mgc:  1.205114 
#> outer mgc:  1.204542 
#> outer mgc:  0.4461172 
#> outer mgc:  0.452669 
#> outer mgc:  0.9735659 
#> outer mgc:  0.9809175 
#> outer mgc:  0.01533163 
#> outer mgc:  0.01538453 
#> outer mgc:  2.779095 
#> outer mgc:  2.770233 
#> outer mgc:  3.013237 
#> outer mgc:  3.005858 
#> outer mgc:  7.590464 
#> outer mgc:  7.592583 
#> outer mgc:  10.9211 
#> outer mgc:  10.92415 
#> outer mgc:  3.537771 
#> outer mgc:  3.54227 
#> outer mgc:  0.304755 
#> outer mgc:  0.311397 
#> outer mgc:  0.1706194 
#> outer mgc:  0.1635631 
#> outer mgc:  0.8772995 
#> outer mgc:  0.8700635 
#> outer mgc:  1.522929 
#> outer mgc:  1.521785 
#> outer mgc:  2.873678 
#> outer mgc:  2.881836 
#> outer mgc:  1.06335 
#> outer mgc:  1.06258 
#> outer mgc:  1.233552 
#> outer mgc:  1.225529 
#> outer mgc:  8.522805 
#> outer mgc:  8.524982 
#> outer mgc:  6.760632 
#> outer mgc:  6.763145 
#> outer mgc:  1.842182 
#> outer mgc:  1.847702 
#> outer mgc:  4.006611 
#> outer mgc:  3.997721 
#> outer mgc:  4.248198 
#> outer mgc:  4.23627 
#> outer mgc:  1.268859 
#> outer mgc:  1.261518 
#> outer mgc:  0.1612047 
#> outer mgc:  0.1680634 
#> outer mgc:  0.2803039 
#> outer mgc:  0.2869764 
#> outer mgc:  4.211275 
#> outer mgc:  4.216344 
#> outer mgc:  3.41939 
#> outer mgc:  3.424782 
#> outer mgc:  3.232833 
#> outer mgc:  3.237916 
#> outer mgc:  1.486452 
#> outer mgc:  1.491889 
#> outer mgc:  0.03107697 
#> outer mgc:  0.03785861 
#> outer mgc:  3.450172 
#> outer mgc:  3.45573 
#> outer mgc:  2.65533 
#> outer mgc:  2.662088 
#> outer mgc:  1.17144 
#> outer mgc:  1.177339 
#> outer mgc:  3.432255 
#> outer mgc:  3.437952 
#> outer mgc:  3.82763 
#> outer mgc:  3.833155 
#> outer mgc:  0.7330496 
#> outer mgc:  0.7393711 
#> outer mgc:  0.2860441 
#> outer mgc:  0.2790369 
#> outer mgc:  0.534149 
#> outer mgc:  0.533641 
#> outer mgc:  5.023069 
#> outer mgc:  5.027728 
#> outer mgc:  3.23317 
#> outer mgc:  3.239002 
#> outer mgc:  5.141974 
#> outer mgc:  5.148815 
#> outer mgc:  0.6807245 
#> outer mgc:  0.6868668 
#> outer mgc:  3.708177 
#> outer mgc:  3.712172 
#> outer mgc:  1.022393 
#> outer mgc:  1.028466 
#> outer mgc:  5.778456 
#> outer mgc:  5.784651 
#> outer mgc:  4.67385 
#> outer mgc:  4.679382 
#> outer mgc:  2.909998 
#> outer mgc:  2.915709 
#> outer mgc:  0.3652791 
#> outer mgc:  0.3647707 
#> outer mgc:  0.8994371 
#> outer mgc:  0.8927196 
#> outer mgc:  1.174924 
#> outer mgc:  1.182068 
#> outer mgc:  2.058319 
#> outer mgc:  2.064897 
#> outer mgc:  5.807428 
#> outer mgc:  5.810696 
#> outer mgc:  3.950383 
#> outer mgc:  3.954096 
#> outer mgc:  12.19628 
#> outer mgc:  12.19668 
#> outer mgc:  7.180289 
#> outer mgc:  7.182716 
#> outer mgc:  1.688183 
#> outer mgc:  1.693474 
#> outer mgc:  10.89622 
#> outer mgc:  10.89814 
#> outer mgc:  6.896194 
#> outer mgc:  6.903373 
#> outer mgc:  6.130411 
#> outer mgc:  6.138583 
#> outer mgc:  4.239735 
#> outer mgc:  4.246409 
#> outer mgc:  8.230351 
#> outer mgc:  8.234589 
#> outer mgc:  4.239054 
#> outer mgc:  4.252703 
#> outer mgc:  5.644993 
#> outer mgc:  5.658023 
#> outer mgc:  5.7308 
#> outer mgc:  5.738049 
#> outer mgc:  5.531328 
#> outer mgc:  5.535201 
#> outer mgc:  1.420927 
#> outer mgc:  1.419385 
#> outer mgc:  1.021527 
#> outer mgc:  1.020361 
#> outer mgc:  3.468327 
#> outer mgc:  3.47458 
#> outer mgc:  1.757969 
#> outer mgc:  1.763484 
#> outer mgc:  1.904011 
#> outer mgc:  1.903157 
#> outer mgc:  1.730432 
#> outer mgc:  1.736458 
#> outer mgc:  1.752873 
#> outer mgc:  1.759543 
#> outer mgc:  2.815207 
#> outer mgc:  2.819144 
#> outer mgc:  0.8573765 
#> outer mgc:  0.8683017 
#> outer mgc:  3.501611 
#> outer mgc:  3.505564 
#> outer mgc:  2.522105 
#> outer mgc:  2.527065 
#> outer mgc:  9.112851 
#> outer mgc:  9.114084 
#> outer mgc:  6.64043 
#> outer mgc:  6.644873 
#> outer mgc:  12.88707 
#> outer mgc:  12.89751 
#> outer mgc:  5.109293 
#> outer mgc:  5.11293 
#> outer mgc:  16.9188 
#> outer mgc:  16.93082 
#> outer mgc:  10.50079 
#> outer mgc:  10.5095 
#> outer mgc:  5.66418 
#> outer mgc:  5.670516 
#> outer mgc:  10.29922 
#> outer mgc:  10.30999 
#> outer mgc:  20.4458 
#> outer mgc:  20.47105 
#> outer mgc:  32.5316 
#> outer mgc:  32.64331 
#> outer mgc:  24.96667 
#> outer mgc:  25.04121 
#> outer mgc:  3.503004 
#> outer mgc:  3.507844 
#> outer mgc:  1.510851 
#> outer mgc:  1.522215 
#> outer mgc:  13.79604 
#> outer mgc:  13.78759 
#> outer mgc:  0.08780279 
#> outer mgc:  0.08091153 
#> outer mgc:  0.005814213 
#> outer mgc:  0.005814213
```

### Data fits

Inspect predicted vs observed length compositions:

``` r
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
```

![](bet_files/figure-html/lf-diagnostics-1.png)

``` r

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
```

![](bet_files/figure-html/lf-diagnostics-2.png)

Inspect fitted CPUE and catch:

``` r
plot(data$cpue_data$value, col = 2, pch = 16,
     xlab = "Time step", ylab = "CPUE", main = "CPUE: observed vs predicted")
lines(exp(obj$simulate()$cpue_log_obs), lwd = 2, col = "gray70")
lines(obj$report()$cpue_pred, lwd = 2)
```

![](bet_files/figure-html/opt-cpue-1.png)

``` r
sum(obj$report()$catch_pred_ysf - data$catch_obs_ysf)
#> [1] 0.9998036
plot_catch(data = data, obj = obj)
#> [1] "The maximum catch difference was: 7.91827187640592e-06"
```

![](bet_files/figure-html/opt-catch-1.png)

``` r
plot_catch(data = data, obj = obj, plot_resid = TRUE)
#> [1] "The maximum catch difference was: 7.91827187640592e-06"
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
plot(log(data$cpue_data$value), col = 2, pch = 16,
     xlab = "Time step", ylab = "log(CPUE)", main = "Observed vs simulated CPUE")
lines(log(obj$report()$cpue_pred), lwd = 2)
for (i in 1:5) lines(obj$simulate()$cpue_log_obs, col = "gray70")
```

![](bet_files/figure-html/sim-1-1.png)

## One step ahead (OSA) residuals

OSA residuals are a replacement for Pearson residuals:

``` r
osa_cpue <- oneStepPredict(obj = obj, observation.name = "cpue_log_obs",
                           method = "oneStepGeneric", trace = FALSE)
qqnorm(osa_cpue$res)
abline(0, 1)
```

![](bet_files/figure-html/osa1-1.png)

``` r
plot(osa_cpue$res, xlab = "Index", ylab = "OSA residual")
abline(h = c(-2, 0, 2), lty = c(2, 1, 2))
```

![](bet_files/figure-html/osa1-2.png)
