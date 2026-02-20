# Catch data prepared for RTMB modeling

A dataset containing catch observations from MFCL (Multifan-CL)
assessments, standardized and structured for integration with RTMB (R
Template Model Builder) stock assessment models. The data represents
catch-at-time by fishery/fleet over a 67-year period (1952-2018) at
quarterly resolution.

## Usage

``` r
catch_data
```

## Format

A data frame (implied dimensions: 14,868 × 8)

## Source

Prepared from MFCL (Multifan-CL) model outputs and frequency files,
primarily extracted from `model-files/mfcl/v11/` using
`frqit::cateffpen()`.

## Details

### Preparation workflow

The dataset was constructed through the following steps:

1.  **Time-series mapping**: Created a table mapping calendar
    years/months to model timesteps. MFCL operates on quarterly cycles
    (months 2, 5, 8, 11 representing Q1-Q4), with each quarter assigned
    a sequential timestep (1-268) for RTMB integration.

2.  **Catch extraction**: Extracted catch by fishery/fleet from the MFCL
    frequency file using `frqit::cateffpen()`.

3.  **Unit standardization**: Converted catch in numbers (units==2) to
    thousands of fish to match the RTMB model structure.

4.  **Timestep merging**: Merged extracted catch with the timestep
    mapping table to add calendar-to-timestep correspondence.

5.  **Metadata addition**: Added standardized columns for metric type
    and uncertainty.

6.  **Fleet filtering**: Removed fleet 15 (MFCL index/survey fleet) to
    retain only direct catch observations.

### Data structure

A data frame with 14,868 rows and 8 columns:

- year:

  Calendar year (1952-2018). Integer vector.

- month:

  Fishing quarter month (2, 5, 8, 11 representing Jan-Mar, Apr-Jun,
  Jul-Sep, Oct-Dec quarters approximately). Integer vector.

- ts:

  Model timestep number (1-268 for 67 years × 4 quarters per year).
  Sequential numeric identifier for RTMB dynamics integration. Integer
  vector.

- fishery:

  MFCL fleet index (1-14, excluding 15 which is survey fleet). Integer
  vector.

- metric:

  Data type identifier. Character vector with value "catch".

- units:

  Unit code. Integer: 1 = metric tons, 2 = thousands of fish. Fisheries
  1-7, 15 use units=2; fisheries 8-14 use units=1. Mixed units reflect
  MFCL's integration of catch from different fishery monitoring systems.

- value:

  Catch quantity in specified units (metric tons or thousands of fish).
  Numeric vector.

- se:

  Standard error in log-space, calculated as \\\sqrt{\log(1 + CV^2)}\\,
  which approximates the coefficient of variation (CV) for small SEs.
  Currently fixed at 0.01 for all observations. Numeric vector.

### Catch units by fishery

- **Fisheries 1-7, 15**: Units = 2 (thousands of fish). These fisheries
  report catch in numbers.

- **Fisheries 8-14**: Units = 1 (metric tons). These fisheries report
  catch in weight.

This mixed-unit structure reflects MFCL's original data integration from
multiple fishery monitoring systems with different reporting standards.

## Examples

``` r
if (FALSE) { # \dontrun{
  data(catch_data)
  head(catch_data)
  summary(catch_data)
  # Check available fisheries
  unique(catch_data$fishery)
  # Filter catch for a specific fishery and time period
  catch_data[catch_data$fishery == 1 & catch_data$year == 2000, ]
} # }
```
