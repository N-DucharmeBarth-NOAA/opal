
<!-- README.md is generated from README.Rmd. Please edit that file -->

# opal

<!-- badges: start -->

[![R-CMD-check](https://github.com/N-DucharmeBarth-NOAA/opal/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/N-DucharmeBarth-NOAA/opal/actions/workflows/R-CMD-check.yaml)
[![codecov-main](https://codecov.io/gh/N-DucharmeBarth-NOAA/opal/graph/badge.svg?token=6JY6W5MDDN)](https://codecov.io/gh/N-DucharmeBarth-NOAA/opal)
[![codecov-dev](https://codecov.io/gh/N-DucharmeBarth-NOAA/opal/branch/dev/graph/badge.svg?token=6JY6W5MDDN)](https://codecov.io/gh/N-DucharmeBarth-NOAA/opal)
<!-- badges: end -->

## Overview

**opal**, the **o**pen **p**opulation **a**ssessment **l**ibrary, is an
open-source, modular R package for fisheries stock assessment. It is
built on [RTMB](https://github.com/kaskr/RTMB), which provides automatic
differentiation of the objective function, gradient-based optimization,
and the Laplace approximation for random effects. The full model is
written in R, so it can be read, modified, and extended without working
in C++.

opal builds on the `sbt` package developed by
[Quantifish](https://www.quantifish.co.nz/) for the CCSBT southern
bluefin tuna assessment, and draws on design elements from SS3,
MULTIFAN-CL, WHAM, SPoRC, and CASAL2. It is being developed in support
of **WCPFC Project 123: Scoping the next generation of tuna stock
assessment software**.

> **Development status:** opal is under active development. Case studies
> are illustrative and are not intended to inform management advice.
> Interfaces may change between versions.

## Resources

| Resource | Contents |
|----|----|
| [Project website](https://connect.fisheries.noaa.gov/opal/) | Project background, development plan and governance, the WCPFC SC22 working paper and presentation, and the ’opakapaka case study |
| [Package documentation](https://n-ducharmebarth-noaa.github.io/opal/) | Function reference, vignettes, and changelog (pkgdown) |
| [opal-documentation](https://github.com/N-DucharmeBarth-NOAA/opal-documentation) | Quarto source for the project website and WCPFC documents |
| [Issues](https://github.com/N-DucharmeBarth-NOAA/opal/issues) | Bug reports, feature requests, and development discussion |

## Features

- **Population dynamics**: Age- and season-structured dynamics
  conditioned on observed catch (in weight or numbers). Seasonal harvest
  rates are solved from catch and vulnerable abundance, with a
  `posfun()` penalty that keeps total harvest below one. Beverton–Holt
  recruitment with log-normal deviates and SS3-style bias-adjustment
  ramps, fished or unfished initial equilibrium, and static and dynamic
  depletion.
- **Growth and biology**: von Bertalanffy growth (L1/L2/k) with
  variability in length at age, represented as a
  probability-of-length-at-age matrix. Maturity, weight, and fecundity
  can be supplied at age or at length; length-based inputs are converted
  on the AD tape so gradients propagate when growth is estimated.
- **Selectivity**: Logistic, double-normal, double-Richards, and
  length-based forms, with helpers to convert between SS3 and opal
  selectivity parameterizations.
- **Data and likelihoods**: Log-normal CPUE with multiple indices;
  length and weight compositions with multinomial, Dirichlet, or
  Dirichlet-multinomial likelihoods (weight compositions are predicted
  by rebinning from length); and user-specified priors.
- **Estimation and diagnostics**: Optimization with `nlminb()` and
  `get_bounds()`; estimability checks, correlated-parameter detection,
  and parameter tables (`check_estimability()`, `get_cor_pairs()`,
  `get_par_table()`). Observations are marked with RTMB’s `OBS()`, which
  supports simulation with `obj$simulate()` and one-step-ahead
  residuals.
- **Bayesian inference**: No-U-turn sampling through
  [`SparseNUTS`](https://github.com/noaa-afsc/SparseNUTS), using
  `opal_globals()`.
- **Projections**: Forward projections of dynamics, recruitment
  deviations, and selectivity (`project_dynamics()`,
  `project_rec_devs()`, `project_selectivity()`).
- **Portable fitted models**: An `opal_fit` object stores data, fitted
  parameters, optimizer output, MCMC draws, diagnostics, derived
  results, and provenance, and rebuilds the RTMB objective on demand.

### Case studies

- **’Opakapaka** (*Pristipomoides filamentosus*): a fixed-effects
  assessment benchmarked against its SS3 reference model. It is the
  basis for the Quickstart vignette, the bundled example fit, and the
  package’s numerical regression and simulation self-tests.
- **WCPO bigeye tuna** (*Thunnus obesus*): a high-dimensional quarterly
  assessment (15 fisheries, 40 age classes, length and weight
  compositions) compared against SS3 and MULTIFAN-CL. Development
  workflows are in `dev/vignettes/`.

Both case studies are described in the SC22 working paper on the
[project website](https://connect.fisheries.noaa.gov/opal/).

## Installation

Install the development version of opal from GitHub:

``` r
# install.packages("remotes")
remotes::install_github("N-DucharmeBarth-NOAA/opal")
```

opal depends on [`SparseNUTS`](https://github.com/noaa-afsc/SparseNUTS),
which is not on CRAN. If installation fails because `SparseNUTS` is
unavailable, install `StanEstimators` first (preferably in a fresh R
session):

``` r
install.packages(
  "StanEstimators",
  repos = c("https://andrjohns.r-universe.dev", "https://cloud.r-project.org")
)
```

Then install `SparseNUTS` and retry the opal installation:

``` r
remotes::install_github("noaa-afsc/SparseNUTS")
```

For development, the package environment is pinned with `renv`; run
`renv::restore()` from the repository root.

## Quick start

Fit the bundled ’opakapaka example and store the result as a portable
`opal_fit`. The full walkthrough is in `vignette("quickstart")`.

``` r
library(opal)
library(RTMB)

inputs <- opaka_quickstart_inputs()

obj <- MakeADFun(
  func = cmb(opal_model, inputs$data),
  parameters = inputs$parameters,
  map = inputs$map,
  silent = TRUE
)
bounds <- get_bounds(obj, inputs$parameters)
control <- list(eval.max = 10000, iter.max = 10000)

opt <- nlminb(obj$par, obj$fn, obj$gr,
              lower = bounds$lower, upper = bounds$upper, control = control)
opt <- nlminb(opt$par, obj$fn, obj$gr,
              lower = bounds$lower, upper = bounds$upper, control = control)

fit <- opal_fit(
  data = inputs$data,
  obj = obj,
  opt = opt,
  bounds = list(lower = bounds$lower, upper = bounds$upper),
  control = control,
  estimability = check_estimability(obj),
  metadata = list(stock = "Opakapaka")
)
summary(fit)
```

To review a fit without refitting, open the bundled example
(`vignette("opakapaka-fit-review")`):

``` r
fit <- read_opal_fit(
  system.file("extdata", "opaka_quickstart_fit.rds", package = "opal"),
  strict = TRUE,
  integrity = "portable"
)
```

The bundled fit may have been saved under a different R version.
Portable integrity permits that difference while still checking model
compatibility and requiring the rebuilt objective to match the saved
value.

## Saving fitted models

`opal_fit` does not serialize the transient RTMB objective.
`read_opal_fit()` validates the saved object against the current
scientific model version, and `opal_fit_object()` rebuilds the objective
from the stored data, parameters, map, and random-effect specification,
then checks that it reproduces the saved objective value.

``` r
save_opal_fit(fit, "fit.rds")
fit <- read_opal_fit("fit.rds", strict = TRUE)

obj <- opal_fit_object(fit)    # rebuilt objective, cached for the session
report <- opal_fit_report(fit) # REPORT() output at the fitted parameters

# Attach later results without changing the fitted state
fit <- update_opal_fit(
  fit,
  mcmc = mcmc_fit,
  derived = list(retrospective = retrospective_results)
)
```

Fits saved under a different R version can be read with
`integrity = "portable"`, which accepts the change provided the rebuilt
objective matches. `opal_as_tmbfit(fit)` converts stored MCMC output to
the `SparseNUTS` `tmbfit` structure for existing plotting and diagnostic
functions.

## Contributing

Contributions are welcome. Please read
[CONTRIBUTING.md](.github/CONTRIBUTING.md) before opening a pull
request. In brief:

- Open an [issue](https://github.com/N-DucharmeBarth-NOAA/opal/issues)
  to discuss substantial changes before starting work.
- Submit pull requests against the `dev` branch; `main` is the stable
  release branch.
- Any change that alters model numerics must follow the regeneration
  protocol described in
  [CONTRIBUTING.md](.github/CONTRIBUTING.md#numerical-changes): show the
  golden test failing, document the change in `NEWS.md`, bump
  `.opal_model_scientific_version`, and regenerate the test reference
  and bundled example fit.

## License

opal is released under the GNU General Public License, version 3 or
later. See [LICENSE](LICENSE).

## Citation

A formal package citation is not yet available. In the meantime, please
cite the WCPFC SC22 working paper describing opal, available from the
[project website](https://connect.fisheries.noaa.gov/opal/).

## Disclaimer

“The United States Department of Commerce (DOC) GitHub project code is
provided on an ‘as is’ basis and the user assumes responsibility for its
use. DOC has relinquished control of the information and no longer has
responsibility to protect the integrity, confidentiality, or
availability of the information. Any claims against the Department of
Commerce stemming from the use of its GitHub project will be governed by
all applicable Federal law. Any reference to specific commercial
products, processes, or services by service mark, trademark,
manufacturer, or otherwise, does not constitute or imply their
endorsement, recommendation or favoring by the Department of Commerce.
The Department of Commerce seal and logo, or the seal and logo of a DOC
bureau, shall not be used in any manner to imply endorsement of any
commercial product or activity by DOC or the United States Government.”

------------------------------------------------------------------------

<a href="https://www.fisheries.noaa.gov/"><img src="man/figures/noaa-fisheries-rgb-2line-horizontal-small.png" alt="NOAA Fisheries" height="75"/></a>

[U.S. Department of Commerce](https://www.commerce.gov/) \| [National
Oceanographic and Atmospheric Administration](https://www.noaa.gov) \|
[NOAA Fisheries](https://www.fisheries.noaa.gov/)
