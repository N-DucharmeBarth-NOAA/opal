# The BET model

## Introduction

The `opal` package is an R package that contains example fisheries data
from western and central pacific bigeye tuna (BET). This page provides
examples using the `opal` BET model.

## Load inputs

Load the `opal` package and the `RTMB` dependency. The `ggplot2` package
is used for plotting.

``` r
library(ggplot2)
library(dplyr)
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

## Model setup

### Parameters

Define the initial parameter values. Growth parameters (`log_L1`,
`log_L2`, `log_k`, `log_CV1`, `log_CV2`) and selectivity (`par_sel`) are
initialised at reasonable starting values. Load these from the bundled
`wcpo_bet_parameters` data object:

``` r
data(wcpo_bet_parameters)

parameters <- list(
  log_B0 = as.numeric(wcpo_bet_parameters$log_B0),
  log_h = as.numeric(wcpo_bet_parameters$log_h),
  log_sigma_r = as.numeric(wcpo_bet_parameters$log_sigma_r),
  log_cpue_q = as.numeric(wcpo_bet_parameters$log_cpue_q),
  cpue_creep = as.numeric(wcpo_bet_parameters$cpue_creep),
  log_cpue_tau = log(0.1),
  log_cpue_omega = as.numeric(wcpo_bet_parameters$log_cpue_omega),
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
#> [1] 92336.23
```

### Parameter map

Use RTMB’s `map` option to turn parameters on/off. Parameters mapped to
`factor(NA)` are fixed at their initial values:

``` r
map_sel <- matrix(NA, nrow(parameters$par_sel), ncol(parameters$par_sel))
map_rdev <- rep(NA, length(parameters$rdev_y))

map <- list(
  # log_B0 = factor(NA),
  log_h = factor(NA),
  log_sigma_r = factor(NA),
  # log_cpue_q = factor(NA),
  cpue_creep = factor(NA),
  log_cpue_tau = factor(NA),
  log_cpue_omega = factor(NA),
  # rdev_y = as.factor(map_rdev),
  par_sel = as.factor(map_sel),
  log_L1  = factor(NA),
  log_L2  = factor(NA),
  log_k   = factor(NA),
  log_CV1 = factor(NA),
  log_CV2 = factor(NA)
)
```

### Build the AD object

Using the `data`, the `parameters`, the parameter `map`, and the model
(`bet_model`), the AD object is created using RTMB’s `MakeADFun`
function:

``` r
obj <- MakeADFun(func = cmb(bet_model, data), parameters = parameters, map = map)
unique(names(obj$par))
#> [1] "log_B0"     "log_cpue_q" "rdev_y"
obj$fn()
#> [1] 522.3479
obj$gr()
#> outer mgc:  655.9751
#>         [,1]     [,2]      [,3]      [,4]      [,5]      [,6]      [,7]
#> [1,] 63.7911 655.9751 -11.25735 -13.78893 -9.534916 -9.631811 -18.07066
#>           [,8]      [,9]     [,10]     [,11]     [,12]     [,13]     [,14]
#> [1,] -13.36007 -13.75146 -17.17998 -16.41135 -7.876104 -7.872097 -26.18511
#>          [,15]    [,16]     [,17]     [,18]     [,19]     [,20]     [,21]
#> [1,] -10.84389 -18.6981 -16.80455 -10.55175 -6.476894 -15.40873 -20.03406
#>          [,22]     [,23]    [,24]     [,25]     [,26]     [,27]     [,28]
#> [1,] -9.011048 -13.20024 -11.0775 -8.638446 -6.687169 -6.926078 -8.702665
#>          [,29]    [,30]     [,31]     [,32]     [,33]    [,34]     [,35]
#> [1,] -10.48124 -13.4115 -8.016917 -9.036088 -9.778807 -11.6629 -6.161929
#>          [,36]     [,37]     [,38]     [,39]    [,40]     [,41]     [,42]
#> [1,] -13.52363 -9.260501 -7.509813 -10.00469 -8.45336 -19.22816 -7.701762
#>          [,43]     [,44]     [,45]     [,46]     [,47]    [,48]     [,49]
#> [1,] -7.429175 -9.678863 -14.82738 -7.265268 -6.477925 -8.21443 -6.811462
#>          [,50]     [,51]     [,52]     [,53]    [,54]     [,55]     [,56]
#> [1,] -7.404575 -6.978893 -9.358331 -6.371974 -7.12056 -13.37546 -4.910538
#>          [,57]     [,58]     [,59]     [,60]     [,61]     [,62]     [,63]
#> [1,] -4.930101 -5.393831 -4.104039 -3.807634 -3.745984 -3.539606 -6.624154
#>          [,64]     [,65]     [,66]     [,67]     [,68]     [,69]     [,70]
#> [1,] -2.853447 -2.540404 -2.171635 -4.343456 -1.870132 -1.542732 -1.569754
#>          [,71]     [,72]     [,73]     [,74]      [,75]     [,76]      [,77]
#> [1,] -2.222998 -1.181415 -1.564512 -1.139099 -0.9027357 -2.822676 -0.2199903
#>           [,78]      [,79]      [,80]     [,81]     [,82]     [,83]     [,84]
#> [1,] -0.7827129 -0.4187449 0.08235421 -2.610672 0.5154272 -3.202673 -0.774114
#>           [,85]     [,86]     [,87]     [,88]     [,89]     [,90]     [,91]
#> [1,] -0.2704935 -1.176879 -1.938762 -2.372391 -3.239055 -3.172927 -2.835908
#>          [,92]    [,93]     [,94]     [,95]     [,96]    [,97]     [,98]
#> [1,] -2.590655 -2.76712 -2.048967 -1.766134 -1.700995 -1.57042 -1.788408
#>          [,99]    [,100]   [,101]   [,102]   [,103]    [,104]   [,105]   [,106]
#> [1,] -1.304719 0.9793638 1.503228 5.140861 3.969402 -1.295671 4.010146 1.022691
#>        [,107]    [,108]   [,109]   [,110]   [,111]   [,112]    [,113]   [,114]
#> [1,] 1.194639 -1.928964 2.002905 1.361108 2.983877 3.221787 -0.616015 9.173514
#>        [,115]  [,116]  [,117]  [,118]   [,119]   [,120]   [,121]   [,122]
#> [1,] 3.944085 7.23007 8.88986 2.53452 3.939087 9.481959 5.418169 4.662833
#>        [,123]   [,124]     [,125]    [,126]   [,127]   [,128]   [,129]   [,130]
#> [1,] 17.65694 1.133542 -0.2114242 0.4775341 9.305087 12.03776 4.399057 16.26963
#>          [,131]    [,132]   [,133]   [,134]   [,135]  [,136]   [,137]   [,138]
#> [1,] -0.8088866 0.8202585 2.592461 2.865422 2.304695 1.64971 7.931158 1.537033
#>          [,139]    [,140]   [,141]   [,142]   [,143]   [,144]    [,145]
#> [1,] -0.7883139 -2.036391 14.25958 2.483436 1.597492 5.022001 0.9194292
#>        [,146]   [,147]   [,148]  [,149]    [,150]     [,151]   [,152]   [,153]
#> [1,] 1.835415 2.135943 5.435077 1.32418 0.9917521 -0.1023073 6.148105 2.261003
#>        [,154]   [,155]   [,156]   [,157]   [,158]   [,159]   [,160]   [,161]
#> [1,] 7.653568 3.285446 3.917342 2.580898 4.303675 6.324247 3.126986 3.249181
#>        [,162]   [,163]   [,164]  [,165]     [,166]   [,167]   [,168]    [,169]
#> [1,] 1.144573 3.126676 1.879378 2.06351 0.01057533 1.542719 6.622598 -2.189544
#>        [,170]   [,171]    [,172]   [,173]   [,174]   [,175]  [,176]   [,177]
#> [1,] 6.831933 2.994982 0.4446768 4.612368 3.734476 3.382808 3.25774 3.857272
#>        [,178]   [,179]   [,180]     [,181]   [,182]   [,183]   [,184]   [,185]
#> [1,] 3.151026 2.749792 1.125806 -0.1384281 16.75225 5.810109 9.222575 8.499391
#>       [,186]   [,187]   [,188]  [,189]   [,190]   [,191]   [,192]   [,193]
#> [1,] 1.08887 20.18758 10.15954 4.44349 4.160749 9.601755 0.233639 7.065445
#>         [,194]   [,195]   [,196]   [,197]   [,198]   [,199]   [,200]   [,201]
#> [1,] -1.510034 3.997019 2.107863 1.768163 5.095416 3.106601 9.354154 9.335153
#>       [,202]   [,203]   [,204]   [,205]   [,206]   [,207]   [,208]    [,209]
#> [1,] 15.7247 12.10172 10.06144 11.57801 3.541204 5.788505 3.524187 0.1072272
#>        [,210]   [,211]   [,212]   [,213]   [,214]   [,215]   [,216]   [,217]
#> [1,] 6.600613 22.64246 25.05418 24.42968 9.316963 19.87895 4.475798 14.11842
#>        [,218]   [,219]   [,220]   [,221]   [,222]   [,223]  [,224]   [,225]
#> [1,] 17.98253 11.83972 11.66908 3.079964 13.18136 1.736702 8.64408 7.047542
#>        [,226]   [,227]   [,228]   [,229]   [,230]   [,231]   [,232]   [,233]
#> [1,] 8.674707 11.76245 11.10025 10.19942 14.92978 4.723927 10.76472 6.921362
#>          [,234]   [,235]   [,236]   [,237]   [,238]    [,239]   [,240]   [,241]
#> [1,] -0.9621026 8.099005 24.56407 6.287994 7.257496 0.8778224 4.335976 5.891532
#>        [,242]   [,243]   [,244]   [,245]   [,246]  [,247]   [,248]   [,249]
#> [1,] 2.616097 15.02018 7.147951 2.201741 1.054701 3.59571 6.544637 3.999719
#>        [,250]   [,251]     [,252]  [,253]   [,254]   [,255]   [,256]   [,257]
#> [1,] 4.120162 4.614858 -0.7223892 6.84508 4.845268 4.360661 1.624197 10.24088
#>        [,258]   [,259]   [,260]   [,261]  [,262]    [,263]   [,264]   [,265]
#> [1,] 2.141329 4.395949 3.760758 2.030657 4.47666 0.5880699 1.322109 2.399602
#>         [,266]     [,267]    [,268]    [,269]    [,270]
#> [1,] -1.102332 -0.7025284 -2.285256 -2.101558 0.6256006
```

Inspect initial model outputs:

``` r
plot(obj$report()$spawning_biomass_y, type = "l",
     xlab = "Time step", ylab = "Spawning biomass",
     main = "Initial spawning biomass trajectory")
```

![](bet_files/figure-html/init-checks-1.png)

``` r

plot_catch(data = data, obj = obj)
#> [1] "The maximum catch difference was: 6.55189978715498e-07"
```

![](bet_files/figure-html/init-checks-2.png)

### Parameter bounds

``` r
Lwr <- rep(-Inf, length(obj$par))
Upr <- rep(Inf, length(obj$par))
Lwr[grep("log_B0", names(obj$par))] <- log(1)
Upr[grep("log_B0", names(obj$par))] <- log(exp(20))
Lwr[grep("log_cpue_q", names(obj$par))] <- log(0.1)
Upr[grep("log_cpue_q", names(obj$par))] <- log(10)
Lwr[grep("rdev_y", names(obj$par))] <- rep(-5, length(parameters$rdev_y))
Upr[grep("rdev_y", names(obj$par))] <- rep(5, length(parameters$rdev_y))
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
opt <- nlminb(start = opt$par, objective = obj$fn, gradient = obj$gr,
              hessian = obj$he, lower = Lwr, upper = Upr, control = control)
max(obj$gr())
```

Compare initial and estimated parameter values:

``` r
data.frame(init = unlist(parameters), 
           value = unlist(obj$env$parList(obj$env$last.par.best))) %>%
  mutate(par = rownames(.)) %>%
  left_join(bounds) %>%
  head()
#>          init       value          par     lower     upper
#> 1 12.00000000 10.21002973       log_B0  0.000000 20.000000
#> 2 -0.05129329 -0.05129329        log_h        NA        NA
#> 3 -0.51082562 -0.51082562  log_sigma_r        NA        NA
#> 4  0.00000000 -0.02397139   log_cpue_q -2.302585  2.302585
#> 5  0.00000000  0.00000000   cpue_creep        NA        NA
#> 6 -2.30258509 -2.30258509 log_cpue_tau        NA        NA
```

Inspect fitted CPUE and catch:

``` r
plot(data$cpue_data$value, col = 2, pch = 16,
     xlab = "Time step", ylab = "CPUE", main = "CPUE: observed vs predicted")
lines(exp(obj$simulate()$cpue_log_obs), lwd = 2, col = "gray70")
lines(obj$report()$cpue_pred, lwd = 2)
```

![](bet_files/figure-html/opt-catch-1.png)

``` r
sum(obj$report()$catch_pred_ysf - data$catch_obs_ysf)
#> [1] 0.9998737
plot_catch(data = data, obj = obj)
#> [1] "The maximum catch difference was: 2.13476050703321e-06"
```

![](bet_files/figure-html/opt-cpue-1.png)

``` r
plot_catch(data = data, obj = obj, plot_resid = TRUE)
#> [1] "The maximum catch difference was: 2.13476050703321e-06"
```

![](bet_files/figure-html/opt-cpue-2.png)

### Selectivity

Visualise the selectivity curves by fleet:

``` r
rep <- obj$report()
sel_fya <- rep$sel_fya

fleet_names <- c("F01_LL.NORTH", "F02_LL.US", "F03_LL.OFFSH",
                 "F04_LL.EQUAT", "F05_LL.WEST", "F06_LL.SOUTH", "F07_LL.AUS",
                 "F08_PS.ASSOC", "F09_PS.UNASS", "F10_DOM.MISC",
                 "F11_DOM.HL", "F12_JP.PS.N", "F13_JP.PL", "F14_EQ.PL",
                 "S01_INDEX")

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

### Diagnostics

Check that all parameters are estimable using the `check_estimability`
function:

``` r
check_estimability(obj = obj)
#> outer mgc:  3.595992e-12 
#> outer mgc:  7.794049 
#> outer mgc:  5.160031 
#> outer mgc:  5.741698 
#> outer mgc:  5.741698 
#> outer mgc:  0.02762746 
#> outer mgc:  0.02759883 
#> outer mgc:  0.02810483 
#> outer mgc:  0.02807635 
#> outer mgc:  0.02860791 
#> outer mgc:  0.02857893 
#> outer mgc:  0.02975023 
#> outer mgc:  0.02972006 
#> outer mgc:  0.0298251 
#> outer mgc:  0.02979406 
#> outer mgc:  0.03130944 
#> outer mgc:  0.03127751 
#> outer mgc:  0.03150733 
#> outer mgc:  0.03147524 
#> outer mgc:  0.03169537 
#> outer mgc:  0.03166309 
#> outer mgc:  0.03080802 
#> outer mgc:  0.03077602 
#> outer mgc:  0.03162572 
#> outer mgc:  0.03159342 
#> outer mgc:  0.03156911 
#> outer mgc:  0.03153688 
#> outer mgc:  0.03167046 
#> outer mgc:  0.0316381 
#> outer mgc:  0.02959097 
#> outer mgc:  0.02956069 
#> outer mgc:  0.02993642 
#> outer mgc:  0.02990675 
#> outer mgc:  0.02818155 
#> outer mgc:  0.02815369 
#> outer mgc:  0.02665255 
#> outer mgc:  0.02662621 
#> outer mgc:  0.02368975 
#> outer mgc:  0.02366498 
#> outer mgc:  0.02353168 
#> outer mgc:  0.02350827 
#> outer mgc:  0.0223005 
#> outer mgc:  0.02227844 
#> outer mgc:  0.0213247 
#> outer mgc:  0.02130362 
#> outer mgc:  0.01901342 
#> outer mgc:  0.01899275 
#> outer mgc:  0.01981084 
#> outer mgc:  0.01979107 
#> outer mgc:  0.01974936 
#> outer mgc:  0.01972983 
#> outer mgc:  0.01989437 
#> outer mgc:  0.01987472 
#> outer mgc:  0.01880266 
#> outer mgc:  0.01878244 
#> outer mgc:  0.02003778 
#> outer mgc:  0.02001785 
#> outer mgc:  0.02029417 
#> outer mgc:  0.02027415 
#> outer mgc:  0.02050945 
#> outer mgc:  0.02048924 
#> outer mgc:  0.02018664 
#> outer mgc:  0.02016605 
#> outer mgc:  0.02081754 
#> outer mgc:  0.02079699 
#> outer mgc:  0.02108022 
#> outer mgc:  0.02105949 
#> outer mgc:  0.02147688 
#> outer mgc:  0.02145577 
#> outer mgc:  0.02175205 
#> outer mgc:  0.02173028 
#> outer mgc:  0.0224593 
#> outer mgc:  0.02243698 
#> outer mgc:  0.02319538 
#> outer mgc:  0.02317264 
#> outer mgc:  0.02355973 
#> outer mgc:  0.02353666 
#> outer mgc:  0.02323155 
#> outer mgc:  0.02320829 
#> outer mgc:  0.02304847 
#> outer mgc:  0.02302554 
#> outer mgc:  0.02287536 
#> outer mgc:  0.02285296 
#> outer mgc:  0.02257094 
#> outer mgc:  0.02254884 
#> outer mgc:  0.02229777 
#> outer mgc:  0.02227571 
#> outer mgc:  0.022396 
#> outer mgc:  0.02237392 
#> outer mgc:  0.02272323 
#> outer mgc:  0.022701 
#> outer mgc:  0.02324648 
#> outer mgc:  0.02322376 
#> outer mgc:  0.02364574 
#> outer mgc:  0.02362212 
#> outer mgc:  0.02446716 
#> outer mgc:  0.02444292 
#> outer mgc:  0.0250582 
#> outer mgc:  0.02503374 
#> outer mgc:  0.02508387 
#> outer mgc:  0.0250594 
#> outer mgc:  0.02470856 
#> outer mgc:  0.02468414 
#> outer mgc:  0.02437486 
#> outer mgc:  0.02435085 
#> outer mgc:  0.0238723 
#> outer mgc:  0.02384899 
#> outer mgc:  0.02328065 
#> outer mgc:  0.02325792 
#> outer mgc:  0.02255523 
#> outer mgc:  0.02253283 
#> outer mgc:  0.02207247 
#> outer mgc:  0.02205064 
#> outer mgc:  0.02162473 
#> outer mgc:  0.02160358 
#> outer mgc:  0.02123299 
#> outer mgc:  0.02121223 
#> outer mgc:  0.02095015 
#> outer mgc:  0.02092932 
#> outer mgc:  0.02118985 
#> outer mgc:  0.02116888 
#> outer mgc:  0.0216899 
#> outer mgc:  0.02166868 
#> outer mgc:  0.02223035 
#> outer mgc:  0.02220863 
#> outer mgc:  0.02201166 
#> outer mgc:  0.02198951 
#> outer mgc:  0.02178478 
#> outer mgc:  0.02176289 
#> outer mgc:  0.0210625 
#> outer mgc:  0.02104187 
#> outer mgc:  0.01987674 
#> outer mgc:  0.01985726 
#> outer mgc:  0.01850472 
#> outer mgc:  0.01848585 
#> outer mgc:  0.01796921 
#> outer mgc:  0.01795104 
#> outer mgc:  0.01804781 
#> outer mgc:  0.01803005 
#> outer mgc:  0.01835737 
#> outer mgc:  0.0183391 
#> outer mgc:  0.01906321 
#> outer mgc:  0.01904376 
#> outer mgc:  0.02046756 
#> outer mgc:  0.02044693 
#> outer mgc:  0.02185774 
#> outer mgc:  0.02183609 
#> outer mgc:  0.02266597 
#> outer mgc:  0.0226435 
#> outer mgc:  0.02288703 
#> outer mgc:  0.02286393 
#> outer mgc:  0.02259358 
#> outer mgc:  0.02257088 
#> outer mgc:  0.0215267 
#> outer mgc:  0.02150531 
#> outer mgc:  0.02006179 
#> outer mgc:  0.0200419 
#> outer mgc:  0.01856533 
#> outer mgc:  0.01854662 
#> outer mgc:  0.01750568 
#> outer mgc:  0.01748808 
#> outer mgc:  0.01678885 
#> outer mgc:  0.01677213 
#> outer mgc:  0.01643471 
#> outer mgc:  0.0164184 
#> outer mgc:  0.01656892 
#> outer mgc:  0.01655235 
#> outer mgc:  0.0171564 
#> outer mgc:  0.01713933 
#> outer mgc:  0.01794004 
#> outer mgc:  0.01792222 
#> outer mgc:  0.01895983 
#> outer mgc:  0.01894102 
#> outer mgc:  0.01989237 
#> outer mgc:  0.01987247 
#> outer mgc:  0.02080414 
#> outer mgc:  0.02078344 
#> outer mgc:  0.02136109 
#> outer mgc:  0.02133992 
#> outer mgc:  0.02152558 
#> outer mgc:  0.02150422 
#> outer mgc:  0.02134894 
#> outer mgc:  0.0213276 
#> outer mgc:  0.02092865 
#> outer mgc:  0.02090777 
#> outer mgc:  0.0203813 
#> outer mgc:  0.02036101 
#> outer mgc:  0.02005317 
#> outer mgc:  0.02003314 
#> outer mgc:  0.02006957 
#> outer mgc:  0.02004934 
#> outer mgc:  0.02041343 
#> outer mgc:  0.02039289 
#> outer mgc:  0.02059462 
#> outer mgc:  0.02057392 
#> outer mgc:  0.02023784 
#> outer mgc:  0.02021746 
#> outer mgc:  0.01872124 
#> outer mgc:  0.01870174 
#> outer mgc:  0.01760136 
#> outer mgc:  0.01758336 
#> outer mgc:  0.01660747 
#> outer mgc:  0.01659067 
#> outer mgc:  0.01591364 
#> outer mgc:  0.01589749 
#> outer mgc:  0.01527405 
#> outer mgc:  0.01525808 
#> outer mgc:  0.01506308 
#> outer mgc:  0.0150476 
#> outer mgc:  0.01458621 
#> outer mgc:  0.0145713 
#> outer mgc:  0.01398981 
#> outer mgc:  0.01397552 
#> outer mgc:  0.01272776 
#> outer mgc:  0.01271407 
#> outer mgc:  0.01207866 
#> outer mgc:  0.01206614 
#> outer mgc:  0.01142566 
#> outer mgc:  0.01141399 
#> outer mgc:  0.01093307 
#> outer mgc:  0.01092187 
#> outer mgc:  0.01049712 
#> outer mgc:  0.01048591 
#> outer mgc:  0.01069531 
#> outer mgc:  0.0106842 
#> outer mgc:  0.01096984 
#> outer mgc:  0.01095848 
#> outer mgc:  0.01140623 
#> outer mgc:  0.01139437 
#> outer mgc:  0.01180962 
#> outer mgc:  0.01179712 
#> outer mgc:  0.01223443 
#> outer mgc:  0.01222156 
#> outer mgc:  0.01229857 
#> outer mgc:  0.01228559 
#> outer mgc:  0.01215303 
#> outer mgc:  0.01214015 
#> outer mgc:  0.01175118 
#> outer mgc:  0.01173856 
#> outer mgc:  0.0112088 
#> outer mgc:  0.01119687 
#> outer mgc:  0.01053371 
#> outer mgc:  0.01052257 
#> outer mgc:  0.009904207 
#> outer mgc:  0.009893767 
#> outer mgc:  0.009278032 
#> outer mgc:  0.009268171 
#> outer mgc:  0.008905848 
#> outer mgc:  0.008896667 
#> outer mgc:  0.008536971 
#> outer mgc:  0.008528407 
#> outer mgc:  0.008133423 
#> outer mgc:  0.008125561 
#> outer mgc:  0.007645132 
#> outer mgc:  0.007637991 
#> outer mgc:  0.007398273 
#> outer mgc:  0.00739207 
#> outer mgc:  0.00718916 
#> outer mgc:  0.007184012 
#> outer mgc:  0.006880353 
#> outer mgc:  0.006876504 
#> outer mgc:  0.006446132 
#> outer mgc:  0.006443818 
#> outer mgc:  0.006069381 
#> outer mgc:  0.006069062 
#> outer mgc:  0.005681134 
#> outer mgc:  0.005683312 
#> outer mgc:  0.00537566 
#> outer mgc:  0.005368684 
#> outer mgc:  0.005271237 
#> outer mgc:  0.005264455 
#> outer mgc:  0.005200924 
#> outer mgc:  0.005194585 
#> outer mgc:  0.00513219 
#> outer mgc:  0.005126388 
#> outer mgc:  0.005104001 
#> outer mgc:  0.005098565 
#> outer mgc:  0.005053469 
#> outer mgc:  0.005048286 
#> outer mgc:  0.005061828 
#> outer mgc:  0.005056758 
#> outer mgc:  0.005100011 
#> outer mgc:  0.005094835 
#> outer mgc:  0.005126851 
#> outer mgc:  0.005121398 
#> outer mgc:  0.00542231 
#> outer mgc:  0.005418701 
#> outer mgc:  0.006066755 
#> outer mgc:  0.006062858 
#> outer mgc:  0.007021041 
#> outer mgc:  0.007016893 
#> outer mgc:  0.007958359 
#> outer mgc:  0.007953936 
#> outer mgc:  0.008826869 
#> outer mgc:  0.00882219 
#> outer mgc:  0.009355745 
#> outer mgc:  0.009350998 
#> outer mgc:  0.009340247 
#> outer mgc:  0.009335444 
#> outer mgc:  0.009133438 
#> outer mgc:  0.009128366 
#> outer mgc:  0.01026722 
#> outer mgc:  0.009681111 
#> outer mgc:  0.01421542 
#> outer mgc:  0.01285767 
#> outer mgc:  0.01829322 
#> outer mgc:  0.01592722 
#> outer mgc:  0.02526456 
#> outer mgc:  0.02046567 
#> outer mgc:  0.2212568 
#> outer mgc:  0.04729101 
#> outer mgc:  0.01177351 
#> outer mgc:  0.01176565 
#> outer mgc:  0.01246214 
#> outer mgc:  0.01246054 
#> outer mgc:  0.01350969 
#> outer mgc:  0.01350746 
#> outer mgc:  0.0138395 
#> outer mgc:  0.01383603 
#> outer mgc:  0.0139337 
#> outer mgc:  0.01392919 
#> outer mgc:  0.01389454 
#> outer mgc:  0.013889 
#> outer mgc:  0.01408321 
#> outer mgc:  0.01407693 
#> outer mgc:  0.01469526 
#> outer mgc:  0.01468854 
#> outer mgc:  0.01539377 
#> outer mgc:  0.0153865 
#> outer mgc:  0.01617976 
#> outer mgc:  0.01617202 
#> outer mgc:  0.01681233 
#> outer mgc:  0.01680404 
#> outer mgc:  0.01692655 
#> outer mgc:  0.01691741 
#> outer mgc:  0.01704681 
#> outer mgc:  0.01703727 
#> outer mgc:  0.01730443 
#> outer mgc:  0.01729447 
#> outer mgc:  0.01810717 
#> outer mgc:  0.01809657 
#> outer mgc:  0.01955576 
#> outer mgc:  0.01954386 
#> outer mgc:  0.02186468 
#> outer mgc:  0.02185137 
#> outer mgc:  0.02586182 
#> outer mgc:  0.0258472 
#> outer mgc:  0.03145777 
#> outer mgc:  0.03144165 
#> outer mgc:  0.03836433 
#> outer mgc:  0.03834686 
#> outer mgc:  0.04593506 
#> outer mgc:  0.04591743 
#> outer mgc:  0.05005644 
#> outer mgc:  0.05003804 
#> outer mgc:  0.05037716 
#> outer mgc:  0.05035758 
#> outer mgc:  0.0471638 
#> outer mgc:  0.047142 
#> outer mgc:  0.0432265 
#> outer mgc:  0.04320293 
#> outer mgc:  0.04123472 
#> outer mgc:  0.04121071 
#> outer mgc:  0.0418061 
#> outer mgc:  0.04178236 
#> outer mgc:  0.04395702 
#> outer mgc:  0.04393163 
#> outer mgc:  0.04690856 
#> outer mgc:  0.04688191 
#> outer mgc:  0.05303744 
#> outer mgc:  0.05301068 
#> outer mgc:  0.06032949 
#> outer mgc:  0.0603016 
#> outer mgc:  0.06869544 
#> outer mgc:  0.06866645 
#> outer mgc:  0.07634493 
#> outer mgc:  0.07631659 
#> outer mgc:  0.07977887 
#> outer mgc:  0.07974901 
#> outer mgc:  0.07701175 
#> outer mgc:  0.0769772 
#> outer mgc:  0.07545341 
#> outer mgc:  0.07541578 
#> outer mgc:  0.07624481 
#> outer mgc:  0.07620649 
#> outer mgc:  0.07957008 
#> outer mgc:  0.07953152 
#> outer mgc:  0.08147101 
#> outer mgc:  0.08142916 
#> outer mgc:  0.08765489 
#> outer mgc:  0.08761025 
#> outer mgc:  0.09635715 
#> outer mgc:  0.09631201 
#> outer mgc:  0.09878914 
#> outer mgc:  0.09874375 
#> outer mgc:  0.09066831 
#> outer mgc:  0.09062091 
#> outer mgc:  0.08078171 
#> outer mgc:  0.08073153 
#> outer mgc:  0.07272687 
#> outer mgc:  0.0726729 
#> outer mgc:  0.06872288 
#> outer mgc:  0.0686655 
#> outer mgc:  0.0717681 
#> outer mgc:  0.07171009 
#> outer mgc:  0.07800383 
#> outer mgc:  0.07778134 
#> outer mgc:  0.08534276 
#> outer mgc:  0.08529113 
#> outer mgc:  0.0933315 
#> outer mgc:  0.09327859 
#> outer mgc:  0.09744778 
#> outer mgc:  0.09739474 
#> outer mgc:  0.09404602 
#> outer mgc:  0.09399055 
#> outer mgc:  0.09176721 
#> outer mgc:  0.09171078 
#> outer mgc:  0.09397378 
#> outer mgc:  0.0939161 
#> outer mgc:  0.0989484 
#> outer mgc:  0.09888764 
#> outer mgc:  0.1034547 
#> outer mgc:  0.1033871 
#> outer mgc:  0.1109489 
#> outer mgc:  0.1108774 
#> outer mgc:  0.1193166 
#> outer mgc:  0.1192411 
#> outer mgc:  0.1230887 
#> outer mgc:  0.1230098 
#> outer mgc:  0.119385 
#> outer mgc:  0.1193019 
#> outer mgc:  0.1146076 
#> outer mgc:  0.1145237 
#> outer mgc:  0.1129869 
#> outer mgc:  0.1129015 
#> outer mgc:  0.1132372 
#> outer mgc:  0.1131482 
#> outer mgc:  0.1137906 
#> outer mgc:  0.1136929 
#> outer mgc:  0.1164561 
#> outer mgc:  0.1163502 
#> outer mgc:  0.1266481 
#> outer mgc:  0.1265316 
#> outer mgc:  0.1366754 
#> outer mgc:  0.1365445 
#> outer mgc:  0.1368572 
#> outer mgc:  0.1367135 
#> outer mgc:  0.1327684 
#> outer mgc:  0.1326165 
#> outer mgc:  0.1309959 
#> outer mgc:  0.1308357 
#> outer mgc:  0.130547 
#> outer mgc:  0.130376 
#> outer mgc:  0.128059 
#> outer mgc:  0.1278786 
#> outer mgc:  0.1226086 
#> outer mgc:  0.1224252 
#> outer mgc:  0.1189223 
#> outer mgc:  0.1187341 
#> outer mgc:  0.1195131 
#> outer mgc:  0.1193115 
#> outer mgc:  0.1256355 
#> outer mgc:  0.125404 
#> outer mgc:  0.1284233 
#> outer mgc:  0.1281621 
#> outer mgc:  0.1369664 
#> outer mgc:  0.1366387 
#> outer mgc:  0.1422792 
#> outer mgc:  0.1418418 
#> outer mgc:  0.1322344 
#> outer mgc:  0.1316985 
#> outer mgc:  0.1127646 
#> outer mgc:  0.1121742 
#> outer mgc:  0.1013517 
#> outer mgc:  0.1005376 
#> outer mgc:  0.1011769 
#> outer mgc:  0.09943067 
#> outer mgc:  0.1073404 
#> outer mgc:  0.1031177 
#> outer mgc:  0.1072716 
#> outer mgc:  0.09971096 
#> outer mgc:  0.1187644 
#> outer mgc:  0.1054645 
#> outer mgc:  0.1186366 
#> outer mgc:  0.1112838 
#> outer mgc:  0.1923087 
#> outer mgc:  0.1495574 
#> outer mgc:  0.1459135 
#> outer mgc:  0.1458323 
#> outer mgc:  0.1633655 
#> outer mgc:  0.1632838 
#> outer mgc:  0.1357662 
#> outer mgc:  0.135682 
#> outer mgc:  0.1068614 
#> outer mgc:  0.1067662 
#> outer mgc:  0.08457142 
#> outer mgc:  0.08446327 
#> outer mgc:  0.07409151 
#> outer mgc:  0.07394482 
#> outer mgc:  0.07021158 
#> outer mgc:  0.06996034 
#> outer mgc:  0.06765164 
#> outer mgc:  0.0668997 
#> outer mgc:  0.09055419 
#> outer mgc:  0.06667119 
#> outer mgc:  0.05877811 
#> outer mgc:  0.05874964 
#> outer mgc:  0.05878693 
#> outer mgc:  0.05876087 
#> outer mgc:  0.06142128 
#> outer mgc:  0.06140022 
#> outer mgc:  0.05860144 
#> outer mgc:  0.05858009 
#> outer mgc:  0.05553062 
#> outer mgc:  0.05551064 
#> outer mgc:  0.05078459 
#> outer mgc:  0.05076612 
#> outer mgc:  0.04431982 
#> outer mgc:  0.04430154 
#> outer mgc:  0.03453287 
#> outer mgc:  0.03451406 
#> outer mgc:  0.02690462 
#> outer mgc:  0.02688886 
#> outer mgc:  0.01921414 
#> outer mgc:  0.01920223 
#> outer mgc:  0.0120897 
#> outer mgc:  0.0120813 
#> outer mgc:  0.006862722 
#> outer mgc:  0.006857488 
#> outer mgc:  0.003281701 
#> outer mgc:  0.003279032 
#> outer mgc:  0.00277574 
#> outer mgc:  0.002775722 
#> outer mgc:  0.00277571 
#> outer mgc:  0.002775533 
#> outer mgc:  0.002777148 
#> outer mgc:  0.002777148 
#> outer mgc:  0.002777693 
#> outer mgc:  0.002777693 
#> outer mgc:  0.002777778 
#> outer mgc:  0.002777778 
#>          Param           MLE Param_check
#> 1       log_B0 10.2100297275         Bad
#> 2   log_cpue_q -0.0239713899          OK
#> 3       rdev_y  0.2947498364          OK
#> 4       rdev_y  0.3201184623          OK
#> 5       rdev_y  0.3578416401          OK
#> 6       rdev_y  0.4143144684          OK
#> 7       rdev_y  0.4544750213          OK
#> 8       rdev_y  0.5046777777          OK
#> 9       rdev_y  0.5220028584          OK
#> 10      rdev_y  0.5362673019          OK
#> 11      rdev_y  0.5330519605          OK
#> 12      rdev_y  0.5568928630          OK
#> 13      rdev_y  0.5625147473          OK
#> 14      rdev_y  0.5737690376          OK
#> 15      rdev_y  0.5320143738          OK
#> 16      rdev_y  0.5156607488          OK
#> 17      rdev_y  0.4446630072          OK
#> 18      rdev_y  0.3743849619          OK
#> 19      rdev_y  0.2912691979          OK
#> 20      rdev_y  0.2380222993          OK
#> 21      rdev_y  0.1657125410          OK
#> 22      rdev_y  0.0981384678          OK
#> 23      rdev_y  0.0439054924          OK
#> 24      rdev_y  0.0059109607          OK
#> 25      rdev_y -0.0227714805          OK
#> 26      rdev_y -0.0397972198          OK
#> 27      rdev_y -0.0436558471          OK
#> 28      rdev_y -0.0524344796          OK
#> 29      rdev_y -0.0617543723          OK
#> 30      rdev_y -0.0702135281          OK
#> 31      rdev_y -0.0746205840          OK
#> 32      rdev_y -0.0846504720          OK
#> 33      rdev_y -0.0939668879          OK
#> 34      rdev_y -0.0968323600          OK
#> 35      rdev_y -0.0887344942          OK
#> 36      rdev_y -0.0777983413          OK
#> 37      rdev_y -0.0713132863          OK
#> 38      rdev_y -0.0722072156          OK
#> 39      rdev_y -0.0807972084          OK
#> 40      rdev_y -0.1007261500          OK
#> 41      rdev_y -0.1276200776          OK
#> 42      rdev_y -0.1501109259          OK
#> 43      rdev_y -0.1609280699          OK
#> 44      rdev_y -0.1637825747          OK
#> 45      rdev_y -0.1593386783          OK
#> 46      rdev_y -0.1421057810          OK
#> 47      rdev_y -0.1092630523          OK
#> 48      rdev_y -0.0773149733          OK
#> 49      rdev_y -0.0618156367          OK
#> 50      rdev_y -0.0614520041          OK
#> 51      rdev_y -0.0661456239          OK
#> 52      rdev_y -0.0810709085          OK
#> 53      rdev_y -0.1099407183          OK
#> 54      rdev_y -0.1406268867          OK
#> 55      rdev_y -0.1646051173          OK
#> 56      rdev_y -0.1923388822          OK
#> 57      rdev_y -0.2247101179          OK
#> 58      rdev_y -0.2473230650          OK
#> 59      rdev_y -0.2495942915          OK
#> 60      rdev_y -0.2409672769          OK
#> 61      rdev_y -0.2265188121          OK
#> 62      rdev_y -0.2071824133          OK
#> 63      rdev_y -0.1892899540          OK
#> 64      rdev_y -0.2033233011          OK
#> 65      rdev_y -0.2549404412          OK
#> 66      rdev_y -0.3167030564          OK
#> 67      rdev_y -0.3572783667          OK
#> 68      rdev_y -0.3884406693          OK
#> 69      rdev_y -0.4022852533          OK
#> 70      rdev_y -0.3767500422          OK
#> 71      rdev_y -0.3171384519          OK
#> 72      rdev_y -0.2492586999          OK
#> 73      rdev_y -0.1897409170          OK
#> 74      rdev_y -0.1443976928          OK
#> 75      rdev_y -0.1192926518          OK
#> 76      rdev_y -0.1326909848          OK
#> 77      rdev_y -0.1890943941          OK
#> 78      rdev_y -0.2647439593          OK
#> 79      rdev_y -0.3335119676          OK
#> 80      rdev_y -0.3953225846          OK
#> 81      rdev_y -0.4431878277          OK
#> 82      rdev_y -0.4637295186          OK
#> 83      rdev_y -0.4519706577          OK
#> 84      rdev_y -0.4137612750          OK
#> 85      rdev_y -0.3578025725          OK
#> 86      rdev_y -0.2947964885          OK
#> 87      rdev_y -0.2307968935          OK
#> 88      rdev_y -0.1786194938          OK
#> 89      rdev_y -0.1432145227          OK
#> 90      rdev_y -0.1231435310          OK
#> 91      rdev_y -0.1178181837          OK
#> 92      rdev_y -0.1291317547          OK
#> 93      rdev_y -0.1467952243          OK
#> 94      rdev_y -0.1507936116          OK
#> 95      rdev_y -0.1324817583          OK
#> 96      rdev_y -0.1056201001          OK
#> 97      rdev_y -0.0932372781          OK
#> 98      rdev_y -0.1146049607          OK
#> 99      rdev_y -0.1630080079          OK
#> 100     rdev_y -0.2325781322          OK
#> 101     rdev_y -0.2958706224          OK
#> 102     rdev_y -0.3333109849          OK
#> 103     rdev_y -0.3418751918          OK
#> 104     rdev_y -0.3548084898          OK
#> 105     rdev_y -0.3796651839          OK
#> 106     rdev_y -0.4256852773          OK
#> 107     rdev_y -0.4639516522          OK
#> 108     rdev_y -0.5273689104          OK
#> 109     rdev_y -0.5794862602          OK
#> 110     rdev_y -0.6059451305          OK
#> 111     rdev_y -0.5934382640          OK
#> 112     rdev_y -0.5723416240          OK
#> 113     rdev_y -0.5215316304          OK
#> 114     rdev_y -0.4524088908          OK
#> 115     rdev_y -0.3735767016          OK
#> 116     rdev_y -0.3053035518          OK
#> 117     rdev_y -0.2509932676          OK
#> 118     rdev_y -0.2127047772          OK
#> 119     rdev_y -0.1929754916          OK
#> 120     rdev_y -0.1980195345          OK
#> 121     rdev_y -0.2162600708          OK
#> 122     rdev_y -0.2317209919          OK
#> 123     rdev_y -0.2353123882          OK
#> 124     rdev_y -0.2359811140          OK
#> 125     rdev_y -0.2288465907          OK
#> 126     rdev_y -0.2117704652          OK
#> 127     rdev_y -0.1848643318          OK
#> 128     rdev_y -0.1577508868          OK
#> 129     rdev_y -0.1262491979          OK
#> 130     rdev_y -0.0899024317          OK
#> 131     rdev_y -0.0552489841          OK
#> 132     rdev_y -0.0314308729          OK
#> 133     rdev_y -0.0221930940          OK
#> 134     rdev_y -0.0270671534          OK
#> 135     rdev_y -0.0491269720          OK
#> 136     rdev_y -0.0944076895          OK
#> 137     rdev_y -0.1570052676          OK
#> 138     rdev_y -0.2199588735          OK
#> 139     rdev_y -0.2614350963          OK
#> 140     rdev_y -0.2848940787          OK
#> 141     rdev_y -0.2818196024          OK
#> 142     rdev_y -0.2448927838          OK
#> 143     rdev_y -0.1773947977          OK
#> 144     rdev_y -0.0883506607          OK
#> 145     rdev_y  0.0218553801          OK
#> 146     rdev_y  0.1401948066          OK
#> 147     rdev_y  0.2418182203          OK
#> 148     rdev_y  0.2843607415          OK
#> 149     rdev_y  0.2540182462          OK
#> 150     rdev_y  0.1953973652          OK
#> 151     rdev_y  0.1421351621          OK
#> 152     rdev_y  0.1011538448          OK
#> 153     rdev_y  0.0769109134          OK
#> 154     rdev_y  0.0788468606          OK
#> 155     rdev_y  0.0932329395         Bad
#> 156     rdev_y  0.1274245552          OK
#> 157     rdev_y  0.1292698662          OK
#> 158     rdev_y  0.0987140181          OK
#> 159     rdev_y  0.0507453861          OK
#> 160     rdev_y -0.0016010831          OK
#> 161     rdev_y -0.0505489833          OK
#> 162     rdev_y -0.0905019396          OK
#> 163     rdev_y -0.1161412619          OK
#> 164     rdev_y -0.1236066010          OK
#> 165     rdev_y -0.1241863586          OK
#> 166     rdev_y -0.1306922614          OK
#> 167     rdev_y -0.1478554486          OK
#> 168     rdev_y -0.1725727067          OK
#> 169     rdev_y -0.1860123794          OK
#> 170     rdev_y -0.1702016137          OK
#> 171     rdev_y -0.1171244579          OK
#> 172     rdev_y -0.0337497165          OK
#> 173     rdev_y  0.0743054567          OK
#> 174     rdev_y  0.2005095937          OK
#> 175     rdev_y  0.3290222103          OK
#> 176     rdev_y  0.4304774780          OK
#> 177     rdev_y  0.4386561213          OK
#> 178     rdev_y  0.3498965136          OK
#> 179     rdev_y  0.2044409692          OK
#> 180     rdev_y  0.0624361094          OK
#> 181     rdev_y -0.0487432361          OK
#> 182     rdev_y -0.1264612417          OK
#> 183     rdev_y -0.1566533271          OK
#> 184     rdev_y -0.1424753452          OK
#> 185     rdev_y -0.0988464953          OK
#> 186     rdev_y -0.0268115896          OK
#> 187     rdev_y  0.0545599880          OK
#> 188     rdev_y  0.1226822748          OK
#> 189     rdev_y  0.1555096426          OK
#> 190     rdev_y  0.1503343308          OK
#> 191     rdev_y  0.1347939845          OK
#> 192     rdev_y  0.1274477917          OK
#> 193     rdev_y  0.1485774236          OK
#> 194     rdev_y  0.1969013723          OK
#> 195     rdev_y  0.2762714551          OK
#> 196     rdev_y  0.3541332871          OK
#> 197     rdev_y  0.3532847446          OK
#> 198     rdev_y  0.2578386453          OK
#> 199     rdev_y  0.1459437021          OK
#> 200     rdev_y  0.0604232669          OK
#> 201     rdev_y  0.0116409713          OK
#> 202     rdev_y  0.0087581670          OK
#> 203     rdev_y  0.0554393207          OK
#> 204     rdev_y  0.1224288786          OK
#> 205     rdev_y  0.1693602028          OK
#> 206     rdev_y  0.1653113791          OK
#> 207     rdev_y  0.1248022116          OK
#> 208     rdev_y  0.0921890955          OK
#> 209     rdev_y  0.0888530259          OK
#> 210     rdev_y  0.1198499683          OK
#> 211     rdev_y  0.1675840030          OK
#> 212     rdev_y  0.2181949812          OK
#> 213     rdev_y  0.2494627002          OK
#> 214     rdev_y  0.2394653538          OK
#> 215     rdev_y  0.1885066658          OK
#> 216     rdev_y  0.1202262207          OK
#> 217     rdev_y  0.0606122713          OK
#> 218     rdev_y  0.0308641935          OK
#> 219     rdev_y  0.0313241993          OK
#> 220     rdev_y  0.0482726567          OK
#> 221     rdev_y  0.0756249375          OK
#> 222     rdev_y  0.0994605342          OK
#> 223     rdev_y  0.1023162767          OK
#> 224     rdev_y  0.0882338948          OK
#> 225     rdev_y  0.0655959350          OK
#> 226     rdev_y  0.0411852700          OK
#> 227     rdev_y  0.0178518204          OK
#> 228     rdev_y -0.0025397172          OK
#> 229     rdev_y -0.0099955036          OK
#> 230     rdev_y  0.0091170405          OK
#> 231     rdev_y  0.0626604010          OK
#> 232     rdev_y  0.1260148143          OK
#> 233     rdev_y  0.1768306645          OK
#> 234     rdev_y  0.1733281370          OK
#> 235     rdev_y  0.1212395401          OK
#> 236     rdev_y  0.0603891868          OK
#> 237     rdev_y  0.0062833438          OK
#> 238     rdev_y -0.0427288894          OK
#> 239     rdev_y -0.0582748703         Bad
#> 240     rdev_y -0.0290384383         Bad
#> 241     rdev_y  0.0369063832         Bad
#> 242     rdev_y  0.1697896814         Bad
#> 243     rdev_y  0.4075496896         Bad
#> 244     rdev_y  0.6217180546         Bad
#> 245     rdev_y  0.7500795231         Bad
#> 246     rdev_y  0.5574558419          OK
#> 247     rdev_y  0.3399615055          OK
#> 248     rdev_y  0.1685600656          OK
#> 249     rdev_y  0.0363854697          OK
#> 250     rdev_y -0.0697387647          OK
#> 251     rdev_y -0.1256426662          OK
#> 252     rdev_y -0.1429415131         Bad
#> 253     rdev_y -0.1455397310          OK
#> 254     rdev_y -0.1257326100          OK
#> 255     rdev_y -0.0875285279          OK
#> 256     rdev_y -0.0421088668          OK
#> 257     rdev_y -0.0016265477          OK
#> 258     rdev_y  0.0355035373          OK
#> 259     rdev_y  0.0619023919          OK
#> 260     rdev_y  0.0657579189          OK
#> 261     rdev_y  0.0568316530          OK
#> 262     rdev_y  0.0299943613          OK
#> 263     rdev_y  0.0116226421          OK
#> 264     rdev_y  0.0066989396          OK
#> 265     rdev_y  0.0058030548          OK
#> 266     rdev_y  0.0039326308          OK
#> 267     rdev_y  0.0013547517          OK
#> 268     rdev_y  0.0003079437          OK
#> 269     rdev_y  0.0000306292          OK
#> 270     rdev_y  0.0000000000          OK
```

Calculate standard deviations of all model parameters:

``` r
Report <- sdreport(obj)
#> outer mgc:  3.595992e-12 
#> outer mgc:  7.794049 
#> outer mgc:  5.160031 
#> outer mgc:  5.741698 
#> outer mgc:  5.741698 
#> outer mgc:  0.02762746 
#> outer mgc:  0.02759883 
#> outer mgc:  0.02810483 
#> outer mgc:  0.02807635 
#> outer mgc:  0.02860791 
#> outer mgc:  0.02857893 
#> outer mgc:  0.02975023 
#> outer mgc:  0.02972006 
#> outer mgc:  0.0298251 
#> outer mgc:  0.02979406 
#> outer mgc:  0.03130944 
#> outer mgc:  0.03127751 
#> outer mgc:  0.03150733 
#> outer mgc:  0.03147524 
#> outer mgc:  0.03169537 
#> outer mgc:  0.03166309 
#> outer mgc:  0.03080802 
#> outer mgc:  0.03077602 
#> outer mgc:  0.03162572 
#> outer mgc:  0.03159342 
#> outer mgc:  0.03156911 
#> outer mgc:  0.03153688 
#> outer mgc:  0.03167046 
#> outer mgc:  0.0316381 
#> outer mgc:  0.02959097 
#> outer mgc:  0.02956069 
#> outer mgc:  0.02993642 
#> outer mgc:  0.02990675 
#> outer mgc:  0.02818155 
#> outer mgc:  0.02815369 
#> outer mgc:  0.02665255 
#> outer mgc:  0.02662621 
#> outer mgc:  0.02368975 
#> outer mgc:  0.02366498 
#> outer mgc:  0.02353168 
#> outer mgc:  0.02350827 
#> outer mgc:  0.0223005 
#> outer mgc:  0.02227844 
#> outer mgc:  0.0213247 
#> outer mgc:  0.02130362 
#> outer mgc:  0.01901342 
#> outer mgc:  0.01899275 
#> outer mgc:  0.01981084 
#> outer mgc:  0.01979107 
#> outer mgc:  0.01974936 
#> outer mgc:  0.01972983 
#> outer mgc:  0.01989437 
#> outer mgc:  0.01987472 
#> outer mgc:  0.01880266 
#> outer mgc:  0.01878244 
#> outer mgc:  0.02003778 
#> outer mgc:  0.02001785 
#> outer mgc:  0.02029417 
#> outer mgc:  0.02027415 
#> outer mgc:  0.02050945 
#> outer mgc:  0.02048924 
#> outer mgc:  0.02018664 
#> outer mgc:  0.02016605 
#> outer mgc:  0.02081754 
#> outer mgc:  0.02079699 
#> outer mgc:  0.02108022 
#> outer mgc:  0.02105949 
#> outer mgc:  0.02147688 
#> outer mgc:  0.02145577 
#> outer mgc:  0.02175205 
#> outer mgc:  0.02173028 
#> outer mgc:  0.0224593 
#> outer mgc:  0.02243698 
#> outer mgc:  0.02319538 
#> outer mgc:  0.02317264 
#> outer mgc:  0.02355973 
#> outer mgc:  0.02353666 
#> outer mgc:  0.02323155 
#> outer mgc:  0.02320829 
#> outer mgc:  0.02304847 
#> outer mgc:  0.02302554 
#> outer mgc:  0.02287536 
#> outer mgc:  0.02285296 
#> outer mgc:  0.02257094 
#> outer mgc:  0.02254884 
#> outer mgc:  0.02229777 
#> outer mgc:  0.02227571 
#> outer mgc:  0.022396 
#> outer mgc:  0.02237392 
#> outer mgc:  0.02272323 
#> outer mgc:  0.022701 
#> outer mgc:  0.02324648 
#> outer mgc:  0.02322376 
#> outer mgc:  0.02364574 
#> outer mgc:  0.02362212 
#> outer mgc:  0.02446716 
#> outer mgc:  0.02444292 
#> outer mgc:  0.0250582 
#> outer mgc:  0.02503374 
#> outer mgc:  0.02508387 
#> outer mgc:  0.0250594 
#> outer mgc:  0.02470856 
#> outer mgc:  0.02468414 
#> outer mgc:  0.02437486 
#> outer mgc:  0.02435085 
#> outer mgc:  0.0238723 
#> outer mgc:  0.02384899 
#> outer mgc:  0.02328065 
#> outer mgc:  0.02325792 
#> outer mgc:  0.02255523 
#> outer mgc:  0.02253283 
#> outer mgc:  0.02207247 
#> outer mgc:  0.02205064 
#> outer mgc:  0.02162473 
#> outer mgc:  0.02160358 
#> outer mgc:  0.02123299 
#> outer mgc:  0.02121223 
#> outer mgc:  0.02095015 
#> outer mgc:  0.02092932 
#> outer mgc:  0.02118985 
#> outer mgc:  0.02116888 
#> outer mgc:  0.0216899 
#> outer mgc:  0.02166868 
#> outer mgc:  0.02223035 
#> outer mgc:  0.02220863 
#> outer mgc:  0.02201166 
#> outer mgc:  0.02198951 
#> outer mgc:  0.02178478 
#> outer mgc:  0.02176289 
#> outer mgc:  0.0210625 
#> outer mgc:  0.02104187 
#> outer mgc:  0.01987674 
#> outer mgc:  0.01985726 
#> outer mgc:  0.01850472 
#> outer mgc:  0.01848585 
#> outer mgc:  0.01796921 
#> outer mgc:  0.01795104 
#> outer mgc:  0.01804781 
#> outer mgc:  0.01803005 
#> outer mgc:  0.01835737 
#> outer mgc:  0.0183391 
#> outer mgc:  0.01906321 
#> outer mgc:  0.01904376 
#> outer mgc:  0.02046756 
#> outer mgc:  0.02044693 
#> outer mgc:  0.02185774 
#> outer mgc:  0.02183609 
#> outer mgc:  0.02266597 
#> outer mgc:  0.0226435 
#> outer mgc:  0.02288703 
#> outer mgc:  0.02286393 
#> outer mgc:  0.02259358 
#> outer mgc:  0.02257088 
#> outer mgc:  0.0215267 
#> outer mgc:  0.02150531 
#> outer mgc:  0.02006179 
#> outer mgc:  0.0200419 
#> outer mgc:  0.01856533 
#> outer mgc:  0.01854662 
#> outer mgc:  0.01750568 
#> outer mgc:  0.01748808 
#> outer mgc:  0.01678885 
#> outer mgc:  0.01677213 
#> outer mgc:  0.01643471 
#> outer mgc:  0.0164184 
#> outer mgc:  0.01656892 
#> outer mgc:  0.01655235 
#> outer mgc:  0.0171564 
#> outer mgc:  0.01713933 
#> outer mgc:  0.01794004 
#> outer mgc:  0.01792222 
#> outer mgc:  0.01895983 
#> outer mgc:  0.01894102 
#> outer mgc:  0.01989237 
#> outer mgc:  0.01987247 
#> outer mgc:  0.02080414 
#> outer mgc:  0.02078344 
#> outer mgc:  0.02136109 
#> outer mgc:  0.02133992 
#> outer mgc:  0.02152558 
#> outer mgc:  0.02150422 
#> outer mgc:  0.02134894 
#> outer mgc:  0.0213276 
#> outer mgc:  0.02092865 
#> outer mgc:  0.02090777 
#> outer mgc:  0.0203813 
#> outer mgc:  0.02036101 
#> outer mgc:  0.02005317 
#> outer mgc:  0.02003314 
#> outer mgc:  0.02006957 
#> outer mgc:  0.02004934 
#> outer mgc:  0.02041343 
#> outer mgc:  0.02039289 
#> outer mgc:  0.02059462 
#> outer mgc:  0.02057392 
#> outer mgc:  0.02023784 
#> outer mgc:  0.02021746 
#> outer mgc:  0.01872124 
#> outer mgc:  0.01870174 
#> outer mgc:  0.01760136 
#> outer mgc:  0.01758336 
#> outer mgc:  0.01660747 
#> outer mgc:  0.01659067 
#> outer mgc:  0.01591364 
#> outer mgc:  0.01589749 
#> outer mgc:  0.01527405 
#> outer mgc:  0.01525808 
#> outer mgc:  0.01506308 
#> outer mgc:  0.0150476 
#> outer mgc:  0.01458621 
#> outer mgc:  0.0145713 
#> outer mgc:  0.01398981 
#> outer mgc:  0.01397552 
#> outer mgc:  0.01272776 
#> outer mgc:  0.01271407 
#> outer mgc:  0.01207866 
#> outer mgc:  0.01206614 
#> outer mgc:  0.01142566 
#> outer mgc:  0.01141399 
#> outer mgc:  0.01093307 
#> outer mgc:  0.01092187 
#> outer mgc:  0.01049712 
#> outer mgc:  0.01048591 
#> outer mgc:  0.01069531 
#> outer mgc:  0.0106842 
#> outer mgc:  0.01096984 
#> outer mgc:  0.01095848 
#> outer mgc:  0.01140623 
#> outer mgc:  0.01139437 
#> outer mgc:  0.01180962 
#> outer mgc:  0.01179712 
#> outer mgc:  0.01223443 
#> outer mgc:  0.01222156 
#> outer mgc:  0.01229857 
#> outer mgc:  0.01228559 
#> outer mgc:  0.01215303 
#> outer mgc:  0.01214015 
#> outer mgc:  0.01175118 
#> outer mgc:  0.01173856 
#> outer mgc:  0.0112088 
#> outer mgc:  0.01119687 
#> outer mgc:  0.01053371 
#> outer mgc:  0.01052257 
#> outer mgc:  0.009904207 
#> outer mgc:  0.009893767 
#> outer mgc:  0.009278032 
#> outer mgc:  0.009268171 
#> outer mgc:  0.008905848 
#> outer mgc:  0.008896667 
#> outer mgc:  0.008536971 
#> outer mgc:  0.008528407 
#> outer mgc:  0.008133423 
#> outer mgc:  0.008125561 
#> outer mgc:  0.007645132 
#> outer mgc:  0.007637991 
#> outer mgc:  0.007398273 
#> outer mgc:  0.00739207 
#> outer mgc:  0.00718916 
#> outer mgc:  0.007184012 
#> outer mgc:  0.006880353 
#> outer mgc:  0.006876504 
#> outer mgc:  0.006446132 
#> outer mgc:  0.006443818 
#> outer mgc:  0.006069381 
#> outer mgc:  0.006069062 
#> outer mgc:  0.005681134 
#> outer mgc:  0.005683312 
#> outer mgc:  0.00537566 
#> outer mgc:  0.005368684 
#> outer mgc:  0.005271237 
#> outer mgc:  0.005264455 
#> outer mgc:  0.005200924 
#> outer mgc:  0.005194585 
#> outer mgc:  0.00513219 
#> outer mgc:  0.005126388 
#> outer mgc:  0.005104001 
#> outer mgc:  0.005098565 
#> outer mgc:  0.005053469 
#> outer mgc:  0.005048286 
#> outer mgc:  0.005061828 
#> outer mgc:  0.005056758 
#> outer mgc:  0.005100011 
#> outer mgc:  0.005094835 
#> outer mgc:  0.005126851 
#> outer mgc:  0.005121398 
#> outer mgc:  0.00542231 
#> outer mgc:  0.005418701 
#> outer mgc:  0.006066755 
#> outer mgc:  0.006062858 
#> outer mgc:  0.007021041 
#> outer mgc:  0.007016893 
#> outer mgc:  0.007958359 
#> outer mgc:  0.007953936 
#> outer mgc:  0.008826869 
#> outer mgc:  0.00882219 
#> outer mgc:  0.009355745 
#> outer mgc:  0.009350998 
#> outer mgc:  0.009340247 
#> outer mgc:  0.009335444 
#> outer mgc:  0.009133438 
#> outer mgc:  0.009128366 
#> outer mgc:  0.01026722 
#> outer mgc:  0.009681111 
#> outer mgc:  0.01421542 
#> outer mgc:  0.01285767 
#> outer mgc:  0.01829322 
#> outer mgc:  0.01592722 
#> outer mgc:  0.02526456 
#> outer mgc:  0.02046567 
#> outer mgc:  0.2212568 
#> outer mgc:  0.04729101 
#> outer mgc:  0.01177351 
#> outer mgc:  0.01176565 
#> outer mgc:  0.01246214 
#> outer mgc:  0.01246054 
#> outer mgc:  0.01350969 
#> outer mgc:  0.01350746 
#> outer mgc:  0.0138395 
#> outer mgc:  0.01383603 
#> outer mgc:  0.0139337 
#> outer mgc:  0.01392919 
#> outer mgc:  0.01389454 
#> outer mgc:  0.013889 
#> outer mgc:  0.01408321 
#> outer mgc:  0.01407693 
#> outer mgc:  0.01469526 
#> outer mgc:  0.01468854 
#> outer mgc:  0.01539377 
#> outer mgc:  0.0153865 
#> outer mgc:  0.01617976 
#> outer mgc:  0.01617202 
#> outer mgc:  0.01681233 
#> outer mgc:  0.01680404 
#> outer mgc:  0.01692655 
#> outer mgc:  0.01691741 
#> outer mgc:  0.01704681 
#> outer mgc:  0.01703727 
#> outer mgc:  0.01730443 
#> outer mgc:  0.01729447 
#> outer mgc:  0.01810717 
#> outer mgc:  0.01809657 
#> outer mgc:  0.01955576 
#> outer mgc:  0.01954386 
#> outer mgc:  0.02186468 
#> outer mgc:  0.02185137 
#> outer mgc:  0.02586182 
#> outer mgc:  0.0258472 
#> outer mgc:  0.03145777 
#> outer mgc:  0.03144165 
#> outer mgc:  0.03836433 
#> outer mgc:  0.03834686 
#> outer mgc:  0.04593506 
#> outer mgc:  0.04591743 
#> outer mgc:  0.05005644 
#> outer mgc:  0.05003804 
#> outer mgc:  0.05037716 
#> outer mgc:  0.05035758 
#> outer mgc:  0.0471638 
#> outer mgc:  0.047142 
#> outer mgc:  0.0432265 
#> outer mgc:  0.04320293 
#> outer mgc:  0.04123472 
#> outer mgc:  0.04121071 
#> outer mgc:  0.0418061 
#> outer mgc:  0.04178236 
#> outer mgc:  0.04395702 
#> outer mgc:  0.04393163 
#> outer mgc:  0.04690856 
#> outer mgc:  0.04688191 
#> outer mgc:  0.05303744 
#> outer mgc:  0.05301068 
#> outer mgc:  0.06032949 
#> outer mgc:  0.0603016 
#> outer mgc:  0.06869544 
#> outer mgc:  0.06866645 
#> outer mgc:  0.07634493 
#> outer mgc:  0.07631659 
#> outer mgc:  0.07977887 
#> outer mgc:  0.07974901 
#> outer mgc:  0.07701175 
#> outer mgc:  0.0769772 
#> outer mgc:  0.07545341 
#> outer mgc:  0.07541578 
#> outer mgc:  0.07624481 
#> outer mgc:  0.07620649 
#> outer mgc:  0.07957008 
#> outer mgc:  0.07953152 
#> outer mgc:  0.08147101 
#> outer mgc:  0.08142916 
#> outer mgc:  0.08765489 
#> outer mgc:  0.08761025 
#> outer mgc:  0.09635715 
#> outer mgc:  0.09631201 
#> outer mgc:  0.09878914 
#> outer mgc:  0.09874375 
#> outer mgc:  0.09066831 
#> outer mgc:  0.09062091 
#> outer mgc:  0.08078171 
#> outer mgc:  0.08073153 
#> outer mgc:  0.07272687 
#> outer mgc:  0.0726729 
#> outer mgc:  0.06872288 
#> outer mgc:  0.0686655 
#> outer mgc:  0.0717681 
#> outer mgc:  0.07171009 
#> outer mgc:  0.07800383 
#> outer mgc:  0.07778134 
#> outer mgc:  0.08534276 
#> outer mgc:  0.08529113 
#> outer mgc:  0.0933315 
#> outer mgc:  0.09327859 
#> outer mgc:  0.09744778 
#> outer mgc:  0.09739474 
#> outer mgc:  0.09404602 
#> outer mgc:  0.09399055 
#> outer mgc:  0.09176721 
#> outer mgc:  0.09171078 
#> outer mgc:  0.09397378 
#> outer mgc:  0.0939161 
#> outer mgc:  0.0989484 
#> outer mgc:  0.09888764 
#> outer mgc:  0.1034547 
#> outer mgc:  0.1033871 
#> outer mgc:  0.1109489 
#> outer mgc:  0.1108774 
#> outer mgc:  0.1193166 
#> outer mgc:  0.1192411 
#> outer mgc:  0.1230887 
#> outer mgc:  0.1230098 
#> outer mgc:  0.119385 
#> outer mgc:  0.1193019 
#> outer mgc:  0.1146076 
#> outer mgc:  0.1145237 
#> outer mgc:  0.1129869 
#> outer mgc:  0.1129015 
#> outer mgc:  0.1132372 
#> outer mgc:  0.1131482 
#> outer mgc:  0.1137906 
#> outer mgc:  0.1136929 
#> outer mgc:  0.1164561 
#> outer mgc:  0.1163502 
#> outer mgc:  0.1266481 
#> outer mgc:  0.1265316 
#> outer mgc:  0.1366754 
#> outer mgc:  0.1365445 
#> outer mgc:  0.1368572 
#> outer mgc:  0.1367135 
#> outer mgc:  0.1327684 
#> outer mgc:  0.1326165 
#> outer mgc:  0.1309959 
#> outer mgc:  0.1308357 
#> outer mgc:  0.130547 
#> outer mgc:  0.130376 
#> outer mgc:  0.128059 
#> outer mgc:  0.1278786 
#> outer mgc:  0.1226086 
#> outer mgc:  0.1224252 
#> outer mgc:  0.1189223 
#> outer mgc:  0.1187341 
#> outer mgc:  0.1195131 
#> outer mgc:  0.1193115 
#> outer mgc:  0.1256355 
#> outer mgc:  0.125404 
#> outer mgc:  0.1284233 
#> outer mgc:  0.1281621 
#> outer mgc:  0.1369664 
#> outer mgc:  0.1366387 
#> outer mgc:  0.1422792 
#> outer mgc:  0.1418418 
#> outer mgc:  0.1322344 
#> outer mgc:  0.1316985 
#> outer mgc:  0.1127646 
#> outer mgc:  0.1121742 
#> outer mgc:  0.1013517 
#> outer mgc:  0.1005376 
#> outer mgc:  0.1011769 
#> outer mgc:  0.09943067 
#> outer mgc:  0.1073404 
#> outer mgc:  0.1031177 
#> outer mgc:  0.1072716 
#> outer mgc:  0.09971096 
#> outer mgc:  0.1187644 
#> outer mgc:  0.1054645 
#> outer mgc:  0.1186366 
#> outer mgc:  0.1112838 
#> outer mgc:  0.1923087 
#> outer mgc:  0.1495574 
#> outer mgc:  0.1459135 
#> outer mgc:  0.1458323 
#> outer mgc:  0.1633655 
#> outer mgc:  0.1632838 
#> outer mgc:  0.1357662 
#> outer mgc:  0.135682 
#> outer mgc:  0.1068614 
#> outer mgc:  0.1067662 
#> outer mgc:  0.08457142 
#> outer mgc:  0.08446327 
#> outer mgc:  0.07409151 
#> outer mgc:  0.07394482 
#> outer mgc:  0.07021158 
#> outer mgc:  0.06996034 
#> outer mgc:  0.06765164 
#> outer mgc:  0.0668997 
#> outer mgc:  0.09055419 
#> outer mgc:  0.06667119 
#> outer mgc:  0.05877811 
#> outer mgc:  0.05874964 
#> outer mgc:  0.05878693 
#> outer mgc:  0.05876087 
#> outer mgc:  0.06142128 
#> outer mgc:  0.06140022 
#> outer mgc:  0.05860144 
#> outer mgc:  0.05858009 
#> outer mgc:  0.05553062 
#> outer mgc:  0.05551064 
#> outer mgc:  0.05078459 
#> outer mgc:  0.05076612 
#> outer mgc:  0.04431982 
#> outer mgc:  0.04430154 
#> outer mgc:  0.03453287 
#> outer mgc:  0.03451406 
#> outer mgc:  0.02690462 
#> outer mgc:  0.02688886 
#> outer mgc:  0.01921414 
#> outer mgc:  0.01920223 
#> outer mgc:  0.0120897 
#> outer mgc:  0.0120813 
#> outer mgc:  0.006862722 
#> outer mgc:  0.006857488 
#> outer mgc:  0.003281701 
#> outer mgc:  0.003279032 
#> outer mgc:  0.00277574 
#> outer mgc:  0.002775722 
#> outer mgc:  0.00277571 
#> outer mgc:  0.002775533 
#> outer mgc:  0.002777148 
#> outer mgc:  0.002777148 
#> outer mgc:  0.002777693 
#> outer mgc:  0.002777693 
#> outer mgc:  0.002777778 
#> outer mgc:  0.002777778
```

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

### One step ahead (OSA) residuals

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
