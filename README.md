# wcpfc-bet-single-region
Exploration of alternative model structures using [2020 WCPFC bigeye tuna stock assessment](https://meetings.wcpfc.int/node/11693) as a proof-of-concept.

Any model runs contained within this repository represent preliminary model explorations and should not be used '*as is*' for the consideration of management advice. Only models that have been presented and accepted to the Western and Central Pacific Fisheries Commission ([WCPFC](https://www.wcpfc.int/)) Scientific Committee should be used as the basis for management advice.

## Background

### WCPO Bigeye Tuna Stock

Bigeye tuna (*Thunnus obesus*) in the Western and Central Pacific Ocean (WCPO) are distributed throughout tropical and subtropical waters west of 150°W. The stock is assessed regularly by the WCPFC using the integrated stock assessment software MULTIFAN-CL ([MFCL](https://github.com/PacificCommunity/multifan-cl/)). 

Key characteristics:
- **Distribution**: Throughout the WCPO, between 50°N and 40°S, west of 150°W
- **Fisheries**: Primarily caught by longline gear (targeting larger adults for high-value sashimi markets) and purse seine sets on fish aggregating devices/FADs (incidental catch of juveniles). Additional catches come from pole-and-line and various domestic fisheries in Southeast Asia
- **Biology**: Fast-growing with maximum fork length around ~155 cm; reproductively active from about 80 cm FL with nearly all fish >120 cm FL reproductively mature

### 2020 Diagnostic Case Model

The [2020 bigeye tuna stock assessment](https://meetings.wcpfc.int/node/11693) was conducted using MFCL and used a complex spatial structure with:
- **9 spatial regions** across the WCPO (with the northern boundary of regions 3 and 4 set at 10°N)
- **41 fisheries** defined by gear type, region, and fleet/flag
- **Quarterly time steps** from 1952-2018
- **Quarterly recruitment deviations** (2,403 parameters)
- **Seasonal movement dynamics** with 104 diffusion coefficients estimated between regions
- **Multiple data sources**: catch, effort, CPUE indices, length and weight composition data, and extensive tag-recapture data

## Simplified Single Region Model

### Motivation

In response to concerns about model complexity and sensitivity to tagging data treatment, a simplified single-region assessment model was developed in tandem with the 2020 diagnostic case. The goals were to:
1. Explore the implications of a drastically different set of assumptions about population dynamics
2. Test whether model complexity was leading to inappropriate spatial buffering
3. Provide a starting point for discussions on appropriate model structure
4. Evaluate the sensitivity of stock status estimates to spatial assumptions

### Approach

The initial simplified model implemented in MFCL used a **"fleets-as-areas"** approach where spatial information was retained through fleet definitions rather than explicit spatial regions. This assumes the stock is well-mixed with no capacity for spatial buffering or partitioning of abundance and fishing mortality.

**Key simplifications from the 2020 diagnostic case:**

1. **Spatial structure**: Collapsed from 9 explicit spatial regions to a single, well-mixed population

2. **Fishery aggregation**: Combined the 41 original fisheries into 15 fisheries based on gear type and broad geographic area:
   - Northern longline (Fisheries 1 & 2)
   - US longline (Fishery 3)
   - Offshore longline (Fisheries 5 & 6)
   - Equatorial longlines (Fisheries 4, 8 & 9)
   - Western longline (Fishery 7)
   - Southern longline (Fisheries 11, 12 & 29)
   - Australian longline (Fisheries 10 & 27)
   - Associated purse seines (Fisheries 13, 15, 25, 30, & 24)
   - Unassociated purse seines (14, 16, 26, & 31)
   - Domestic miscellaneous fisheries (Fisheries 17, 23, & 22)
   - Domestic handline (Fishery 18)
   - Northern Japanese purse-seine (Fishery 19)
   - Northern Japanese pole-and-line (Fishery 20)
   - Equatorial pole-and-line (Fisheries 21, 22, & 28)
   - Index fishery (Fisheries 33-41)

3. **Data aggregation**: Catch and size frequency data were aggregated from the appropriate original fisheries

4. **Movement**: Eliminated quarterly movement parameters entirely (removing 104 diffusion coefficients)

5. **CPUE index**: Calculated a single standardized CPUE index over the entire WCPO assessment region

6. **Seasonal dynamics**: Added a seasonal catchability component for the index fishery to compensate for lack of movement

7. **Selectivity**: 
   - Removed non-decreasing penalties from most longline fisheries (except index and domestic handline)
   - Fixed selectivity of first 3 quarterly ages to zero for some longline fisheries
   - Retained non-decreasing selectivity for index fishery (assumed to sample largest fish) and domestic handline

## Stock Synthesis (SS3) Model Development

Building on the initial simplified MFCL model, an equivalent model was developed using Stock Synthesis ([SS3](https://nmfs-ost.github.io/ss3-website/)). **SS3** is a widely used and established integrated stock assessment framework that has been applied by multiple tuna Regional Fisheries Management Organizations (tRFMOs), including the Inter-American Tropical Tuna Commission (IATTC) and Indian Ocean Tuna Commission (IOTC) for assessing tropical tuna stocks. While SS3 is approaching end-of-life as software development winds down, it serves as a valuable "bridge" model between the current MFCL framework and next-generation modeling approaches.

### SS3 Workflow

The SS3 implementation workflow consists of two main stages:

#### Stage 1: Baseline Model Creation ([01-make-baseline-ss3-model.r](code/ss3/01-make-baseline-ss3-model.r))
- Parses MFCL outputs (frequency data, parameters) using FLR4MFCL and frqit packages
- Reads template SS3 files from `00-swpo-mls-base-file/` (fixed reference baseline)
- Maps MFCL biological parameters (growth, maturity, length-weight) into SS3 parameter blocks
- Configures fleets, selectivity patterns, and data components to match MFCL structure
- Writes SS3 input files (`starter.ss`, `data.ss`, `control.ss`, `forecast.ss`) to `01-bet-base/`
- Automatically launches SS3 executable and runs the baseline model

#### Stage 2: Parameter Refinement ([02-fix-sel.r](code/ss3/02-fix-sel.r))
- **New workflow** for iteratively addressing parameters flagged at bounds during optimization
- Parses `Report.sso` to identify parameters hitting bounds (status HI/LO)
- Intelligently expands parameter bounds based on parameter type:
  - Selectivity curve parameters (ascend_se, descend_se): expanded bounds allow more shape flexibility
  - Width parameters: adjusted for improved model convergence
- Updates selectivity parameter matrices in control file
- Reconfigures likelihood weights (lambdas) for all data components (surveys, length composition, weight composition)
- Writes refined SS3 files to `02-fix-sel/` and re-runs model
- Provides diagnostic output for bound expansion strategy

### Fleet Structure & Data Integration

The SS3 model maintains the same 15-fishery aggregation scheme as the MFCL simplified model:
- **Fleet definitions** retain spatial/gear information through explicit fleet naming
- **Selectivity** is specified separately for each fleet (peak, steepness, ascending/descending width for double-normal; inflection and width for logistic)
- **Data components** (catch, length/weight compositions, CPUE indices) are mapped from aggregated MFCL data
- **Likelihood weights** (lambdas) can be adjusted to emphasize or de-emphasize particular data streams in successive model refinements

#### Catch Units

Catch data in the SS3 model are specified in mixed units depending on fleet:
- **Numbers (in 1000s of fish)** for longline fleets (F01–F07: Northern, US, Offshore, Equatorial, Western, Southern, Australian longlines) and the index fishery (S01)
- **Metric tons (mt)** for all other fleets (F08–F14: Purse seines, domestic miscellaneous, domestic handline, Japanese purse-seine, and pole-and-line fleets)

When converting catch data from MFCL (recorded in numbers), fleets specified in numbers are divided by 1,000 to match SS3's convention of recording catch in thousands of fish.

## Alternative Models Framework

Future development of an [RTMB](https://github.com/kaskr/RTMB)-based model is planned. RTMB is a modern statistical modeling language built on Template Model Builder ([TMB](https://github.com/kaskr/adcomp)) that provides a flexible framework for developing bespoke population dynamics models. As a next-generation modeling platform, RTMB enables the implementation of random effects structures and seamless integration with Bayesian modeling approaches, facilitating more robust model diagnostics and improved characterization of uncertainty compared to traditional deterministic optimization approaches.

## Extractor utilities

This repository includes modular extractor utilities and unit tests to standardize and compare MFCL and SS3 model outputs.

### Biomass Extractors

- `extract_ss3_biomass()` (located in `code/ss3/helper-fns/`) reads an SS3 `Report.sso` and extracts spawning biomass (SSB) and Bratio-derived depletion.
- `extract_mfcl_biomass()` (located in `code/mfcl/helper-fns/`) reads MFCL `.rep` output and extracts fished and unfished SSB to compute depletion.

Both extractors return a `data.table` with the same columns: `model, year, ts, season, ssb, ssb_se, depletion, depletion_se`. The functions accept `quarterly=TRUE/FALSE`; when `quarterly=FALSE` they aggregate quarterly values into annual summaries (season=1) so outputs remain directly comparable.

### Length Composition Extractors

- `extract_ss3_length_comp()` (located in `code/ss3/helper-fns/`) reads an SS3 `Report.sso` and extracts observed and expected length composition data fits.
- `extract_mfcl_length_comp()` (located in `code/mfcl/helper-fns/`) reads MFCL `length.fit` files and extracts length composition data.

Both extractors return a `data.table` with identical columns: `id, Fleet, Fleet_name, Used, Kind, Sex, Bin, Obs, Exp, Dev, effN, Nsamp_in, Nsamp_adj`. They aggregate data across time and optionally support bin harmonization using `rebin_composition()` to enable direct cross-model comparison even when bin structures differ.

## Unit testing

Unit tests are under `tests/testthat/` and validate structure, data types, value ranges, and compatibility between extractors. Run the test suite from an R session with either:

```r
devtools::test()
# or
testthat::test_dir("tests/testthat")
```

## License

The code contained in this repository is licensed under the GNU GENERAL PUBLIC LICENSE version 3 ([GPLv3](https://www.gnu.org/licenses/gpl-3.0.html)).

## Disclaimer

This repository is a scientific product and is not official communication of the National Oceanic and Atmospheric Administration, or the United States Department of Commerce. All NOAA GitHub project code is provided on an ‘as is’ basis and the user assumes responsibility for its use. Any claims against the Department of Commerce or Department of Commerce bureaus stemming from the use of this GitHub project will be governed by all applicable Federal law. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or favoring by the Department of Commerce. The Department of Commerce seal and logo, or the seal and logo of a DOC bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by DOC or the United States Government.
