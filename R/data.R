#' Catch data prepared for RTMB modeling
#'
#' A dataset containing catch observations from MFCL (Multifan-CL) assessments,
#' standardized and structured for integration with RTMB (R Template Model Builder)
#' stock assessment models. The data represents catch-at-time by fishery/fleet over
#' a 67-year period (1952-2018) at quarterly resolution.
#'
#' @details
#' ## Preparation workflow
#'
#' The dataset was constructed through the following steps:
#'
#' 1. **Time-series mapping**: Created a table mapping calendar years/months to model
#'    timesteps. MFCL operates on quarterly cycles (months 2, 5, 8, 11 representing
#'    Q1-Q4), with each quarter assigned a sequential timestep (1-268) for RTMB integration.
#'
#' 2. **Catch extraction**: Extracted catch by fishery/fleet from the MFCL frequency file
#'    using `frqit::cateffpen()`.
#'
#' 3. **Unit standardization**: Converted catch in numbers (units==2) to thousands of fish
#'    to match the RTMB model structure.
#'
#' 4. **Timestep merging**: Merged extracted catch with the timestep mapping table to
#'    add calendar-to-timestep correspondence.
#'
#' 5. **Metadata addition**: Added standardized columns for metric type and uncertainty.
#'
#' 6. **Fleet filtering**: Removed fleet 15 (MFCL index/survey fleet) to retain only
#'    direct catch observations.
#'
#' ## Data structure
#'
#' A data frame with 14,868 rows and 8 columns:
#'
#' \describe{
#'   \item{year}{Calendar year (1952-2018). Integer vector.}
#'   \item{month}{Fishing quarter month (2, 5, 8, 11 representing Jan-Mar, Apr-Jun,
#'     Jul-Sep, Oct-Dec quarters approximately). Integer vector.}
#'   \item{ts}{Model timestep number (1-268 for 67 years × 4 quarters per year).
#'     Sequential numeric identifier for RTMB dynamics integration. Integer vector.}
#'   \item{fishery}{MFCL fleet index (1-14, excluding 15 which is survey fleet).
#'     Integer vector.}
#'   \item{metric}{Data type identifier. Character vector with value "catch".}
#'   \item{units}{Unit code. Integer: 1 = metric tons, 2 = thousands of fish.
#'     Fisheries 1-7, 15 use units=2; fisheries 8-14 use units=1. Mixed units reflect
#'     MFCL's integration of catch from different fishery monitoring systems.}
#'   \item{value}{Catch quantity in specified units (metric tons or thousands of fish).
#'     Numeric vector.}
#'   \item{se}{Standard error in log-space, calculated as \eqn{\sqrt{\log(1 + CV^2)}},
#'     which approximates the coefficient of variation (CV) for small SEs.
#'     Currently fixed at 0.01 for all observations. Numeric vector.}
#' }
#'
#' ## Catch units by fishery
#'
#' - **Fisheries 1-7, 15**: Units = 2 (thousands of fish). These fisheries report
#'   catch in numbers.
#'
#' - **Fisheries 8-14**: Units = 1 (metric tons). These fisheries report catch
#'   in weight.
#'
#' This mixed-unit structure reflects MFCL's original data integration from multiple
#' fishery monitoring systems with different reporting standards.
#'
#' @format A data frame (implied dimensions: 14,868 × 8)
#'
#' @source Prepared from MFCL (Multifan-CL) model outputs and frequency files, primarily
#'   extracted from `model-files/mfcl/v11/` using `frqit::cateffpen()`.
#'
#' @examples
#' \dontrun{
#'   data(wcpo_bet_catch_data)
#'   head(wcpo_bet_catch_data)
#'   summary(wcpo_bet_catch_data)
#'   # Check available fisheries
#'   unique(wcpo_bet_catch_data$fishery)
#'   # Filter catch for a specific fishery and time period
#'   wcpo_bet_catch_data[wcpo_bet_catch_data$fishery == 1 & wcpo_bet_catch_data$year == 2000, ]
#' }
#'
"wcpo_bet_catch_data"


#' CPUE (Catch-Per-Unit-Effort) index data prepared for RTMB modeling
#'
#' A dataset containing standardized CPUE observations from MFCL (Multifan-CL)
#' index/survey fleet (fishery 15), prepared for integration with RTMB
#' (R Template Model Builder) stock assessment models. CPUE serves as a
#' relative abundance index for model fitting and diagnostics.
#'
#' @details
#' ## Preparation workflow
#'
#' The dataset was constructed from MFCL fishery 15 (survey fleet) through
#' the following steps:
#'
#' 1. **Index extraction**: Extracted fishery 15 (MFCL index/survey fleet) from
#'    the baseline frequency file using `frqit::cateffpen()`.
#'
#' 2. **CPUE calculation**: Calculated CPUE as catch-per-unit-effort
#'    (catch / effort) from the frequency file.
#'
#' 3. **Normalization**: Normalized CPUE to mean across all time periods
#'    (obs = cpue / mean(cpue)) to facilitate model fitting and interpretation.
#'
#' 4. **Uncertainty estimation**: Converted MFCL penalty column to coefficient
#'    of variation (CV). The MFCL penalty column represents the relationship
#'    \eqn{1/(2 \times CV^2)}, which is inverted to recover CV.
#'
#' 5. **Log-space SE**: Calculated log-space standard error as
#'    \eqn{SE_{log} = \sqrt{\log(1 + CV^2)}}, which approximates CV for small SEs.
#'
#' 6. **Timestep merging**: Merged with timestep table to maintain consistent
#'    temporal indexing with catch data.
#'
#' 7. **Standardization**: Standardized to output column structure matching
#'    catch data for seamless integration.
#'
#' ## Data structure
#'
#' A data frame with typically 268 rows (one per quarterly timestep) and 8 columns:
#'
#' \describe{
#'   \item{year}{Calendar year component (1952-2018). Integer vector.}
#'   \item{month}{Fishing quarter month (2, 5, 8, 11 representing Jan-Mar, Apr-Jun,
#'     Jul-Sep, Oct-Dec quarters approximately). Integer vector.}
#'   \item{ts}{Model timestep number (1-268 for 67 years × 4 quarters per year).
#'     Sequential numeric identifier matching catch data. Integer vector.}
#'   \item{fishery}{Fleet identifier. Always 15 for index/survey fleet.
#'     Integer vector.}
#'   \item{metric}{Data type identifier. Character vector with value "cpue".}
#'   \item{units}{Unit code. Integer: 1 = metric tons, 2 = thousands of fish.
#'     Fishery 15 uses units=2 (index in numbers). Integer vector.}
#'   \item{value}{Normalized CPUE observation (cpue / mean(cpue)).
#'     Dimensionless index where geometric mean across timesteps = 1.
#'     Numeric vector.}
#'   \item{se}{Standard error in log-space, calculated as \eqn{\sqrt{\log(1 + CV^2)}}.
#'     Numeric vector capturing observation-specific uncertainty.}
#' }
#'
#' ## Usage in stock assessment
#'
#' The normalized CPUE values serve as relative abundance indices in the
#' likelihood function. Log-space standard errors represent observation
#' uncertainty and are used to weight likelihood contributions. The
#' normalization ensures that the index mean is 1.0, facilitating
#' interpretation of model-estimated abundance trends.
#'
#' @format A data frame (implied dimensions: 268 × 8)
#'
#' @source Prepared from MFCL (Multifan-CL) model outputs, specifically from
#'   fishery 15 (index fleet) in the baseline frequency file,
#'   extracted from `model-files/mfcl/v11/` using `frqit::cateffpen()`.
#'
#' @examples
#' \dontrun{
#'   data(wcpo_bet_cpue_data)
#'   head(wcpo_bet_cpue_data)
#'   summary(wcpo_bet_cpue_data)
#'   # Verify normalized index (geometric mean should be ~1)
#'   exp(mean(log(wcpo_bet_cpue_data$value)))
#'   # Plot CPUE over time
#'   plot(wcpo_bet_cpue_data$ts, wcpo_bet_cpue_data$value, type = "l")
#' }
#'
"wcpo_bet_cpue_data"
