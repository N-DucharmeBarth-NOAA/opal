# Quickstart Guide

## The `opal` Package

`opal` is an open, modular R package for fisheries stock assessment
built on [RTMB](https://github.com/kaskr/RTMB) for automatic
differentiation and gradient-based optimization.

The core workflow involves:

1.  **Specifying a model** — Define population dynamics, selectivity,
    biological parameters, and data likelihoods
2.  **Running an optimization** — Fit the model to data and estimate key
    parameters

See
[`vignette("bet")`](https://n-ducharmebarth-noaa.github.io/opal/articles/bet.md)
for a detailed worked example.

``` r
library(opal)
```
