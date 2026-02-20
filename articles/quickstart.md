# Quickstart Guide

## The `opal` Package

`opal` is an R package for building fast, transparent, and reproducible
stock assessment models. It leverages
[RTMB](https://kaskr.github.io/adcomp/group__R__interface.html) for
automatic differentiation, enabling efficient optimization and
uncertainty quantification.

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
