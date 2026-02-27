
<!-- README.md is generated from README.Rmd. Please edit that file -->

# opal

<!-- badges: start -->

[![R-CMD-check](https://github.com/N-DucharmeBarth-NOAA/opal/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/N-DucharmeBarth-NOAA/opal/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## Overview

**opal** — the **o**pen **p**opulation **a**ssessment **l**ibrary — is
an open, modular R package for fisheries stock assessment built on
[RTMB](https://github.com/kaskr/RTMB) for automatic differentiation and
gradient-based optimization.

### Features

- **Age-structured population dynamics with length-based processes**:
  Forward simulation using the Baranov catch equation with flexible
  natural mortality, growth, and selectivity options
- **Modular design**: Population model components (growth, selectivity,
  recruitment, likelihoods) are implemented as standalone functions that
  can be combined as needed
- **Multiple likelihood options**: Supports log-normal CPUE, and
  multinomial, Dirichlet, and Dirichlet-multinomial length & weight
  composition likelihoods
- **Parameter estimation and diagnostics**: Gradient-based optimization
  via `nlminb`, OSA residuals, likelihood profiling, estimability
  checks, and parameter correlation summaries
- **Uncertainty quantification and projections (*planned*)**: MCMC
  sampling via `SparseNUTS`, and forward projections from estimated or
  sampled parameter sets
- **Open source**: Released under an open-source license with documented
  functions and reproducible workflows

### Purpose

opal is designed for fisheries stock assessment analysts and population
modelling researchers who need a transparent, reproducible framework for
quantitative population assessment.

## Installation

You can install the development version of opal from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("N-DucharmeBarth-NOAA/opal")
```

Ensure you have RTMB installed:

``` r
install.packages("RTMB")
```

## Quick Start

``` r
library(opal)

# Build a basic population model
# (Example code will depend on your API—update with actual function calls)
```

## Documentation

For detailed examples and function reference, see the
[vignettes](vignettes/) and inline function documentation (`?opal`).

## Contributing

Contributions are welcome. Please open an issue or submit a pull request
with enhancements or bug fixes.

### Branch Structure

opal uses a standard development workflow:

- **`main`**: Stable release branch. Code here represents tested,
  production-ready versions of the package.
- **`dev`**: Active development branch. New features, bug fixes, and
  improvements are developed and tested here before being merged into
  `main`.

**Workflow**: Submit pull requests against the `dev` branch. Once
features are tested and validated, they are incorporated into `main` for
release.

### Using GitHub Copilot

If you use GitHub Copilot to assist with development, please select
`dev` as the base branch from the dropdown menu when assigning work to
Copilot (the default is `main`). This ensures your feature branches and
pull requests follow our development workflow.

## License

See [LICENSE](LICENSE) for terms and conditions.

## Citation

If you use opal in your research, please cite it appropriately. Details
will be provided as the package matures.

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

<a href="https://www.fisheries.noaa.gov/"><img src="man/figures/noaa-fisheries-rgb-2line-horizontal-small.png" alt="NOAA Fisheries" height="75"/>

[U.S. Department of Commerce](https://www.commerce.gov/) \| [National
Oceanographic and Atmospheric Administration](https://www.noaa.gov) \|
[NOAA Fisheries](https://www.fisheries.noaa.gov/)
