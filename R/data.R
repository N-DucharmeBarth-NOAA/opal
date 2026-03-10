#' West Central Pacific Ocean Bigeye Tuna (WCPO BET) Assessment Data
#'
#' A comprehensive list containing biological parameters, catch and CPUE data,
#' and Bayesian prior specifications prepared for RTMB (R Template Model Builder)
#' stock assessment modeling of western central Pacific bigeye tuna.
#'
#' @details
#' ## Data structure (26 list elements)
#'
#' **Age and time dimensions:**
#' \describe{
#'   \item{age_a}{Vector of integer ages (1:40). Length = 40.}
#'   \item{n_age}{Number of age classes (40).}
#'   \item{n_season}{Number of seasons per year (1).}
#'   \item{n_fishery}{Number of fisheries/fleets (15).}
#'   \item{n_year}{Number of years in model (268, equivalent to 67 years × 4 quarters).}
#'   \item{first_yr}{First year designation in model (1, representing 1952).}
#'   \item{last_yr}{Last year designation in model (268, representing 2018).}
#'   \item{years}{Vector of model years (1:268).}
#' }
#'
#' **Length structure:**
#' \describe{
#'   \item{len_bin_start}{First length bin lower edge (10 cm).}
#'   \item{len_bin_width}{Width of each length bin (2 cm).}
#'   \item{n_len}{Total number of length bins (95).}
#' }
#'
#' **Catch data:**
#' \describe{
#'   \item{first_yr_catch}{First year with catch data (1).}
#'   \item{catch_units_f}{Numeric vector of unit codes by fishery (15 elements).
#'     1 = metric tons (fisheries 8-14), 2 = thousands of fish (fisheries 1-7, 15).}
#'   \item{catch_obs_ysf}{Three-dimensional array of observed catch (268 years × 1 season × 15 fisheries).
#'     Dimensions follow (year, season, fishery). Units follow catch_units_f.}
#' }
#'
#' **CPUE data:**
#' \describe{
#'   \item{cpue_switch}{Indicator for CPUE inclusion (1 = included, typically for fishery 15).}
#'   \item{n_index}{Number of distinct CPUE/survey indices (integer).}
#'   \item{cpue_data}{Tibble with 268 rows and 9 columns containing standardized CPUE observations:
#'     \itemize{
#'       \item{year}: Calendar year (1952-2018)
#'       \item{month}: Quarter month (2, 5, 8, 11 for Q1-Q4)
#'       \item{ts}: Sequential timestep (1-268)
#'       \item{fishery}: Fleet identifier (usually 15 for survey fleet)
#'       \item{index}: Integer index identifier (1, 2, ..., n_index)
#'       \item{metric}: Data type ("cpue")
#'       \item{units}: Unit code (typically 2 for fishery 15)
#'       \item{value}: Normalized CPUE observations
#'       \item{se}: Log-space standard error of CPUE observations
#'     }
#'   }
#' }
#'
#' **Biological parameters:**
#' \describe{
#'   \item{lw_a}{Length-weight allometric coefficient (6.48e-05).}
#'   \item{lw_b}{Length-weight allometric exponent (2.78).}
#'   \item{maturity}{Maturity at length vector (95 elements). Ogive-based maturity
#'     probabilities by length bin.}
#'   \item{fecundity}{Fecundity at length vector (95 elements). Relative fecundity by length bin.}
#'   \item{M}{Natural mortality at age vector (40 elements). Age-specific instantaneous
#'     mortality rates.}
#' }
#'
#' **Age bounds:**
#' \describe{
#'   \item{A1}{Minimum age for credibility interval (1).}
#'   \item{A2}{Maximum age for credibility interval (40).}
#' }
#'
#' **Selectivity:**
#' \describe{
#'   \item{sel_type_f}{Selectivity function type by fishery (15 elements). Integer codes
#'     representing selectivity curve shapes (e.g., 2 = double-normal).}
#' }
#'
#' **Priors:**
#' \describe{
#'   \item{priors}{List of 8 Bayesian prior specifications, each containing:
#'     \itemize{
#'       \item{type}: Prior distribution type ("normal", etc.)
#'       \item{par1}: First parameter (mean or location)
#'       \item{par2}: Second parameter (std dev or scale)
#'       \item{index}: Parameter index in the model
#'     }
#'     Priors cover: log_B0, log_cpue_q, par_sel, log_L1, log_L2, log_k, log_CV1, log_CV2.
#'   }
#' }
#'
#' **Length composition data:**
#' \describe{
#'   \item{lf_minbin}{Integer vector (n_fishery). Per-fishery minimum bin index
#'     for lower-tail aggregation. Bins from 1 to `lf_minbin[f]` are summed
#'     into bin `lf_minbin[f]`. Default is 1 (no aggregation).}
#'   \item{lf_maxbin}{Integer vector (n_fishery). Per-fishery maximum bin index
#'     for upper-tail aggregation. Bins from `lf_maxbin[f]` to n_len are summed
#'     into bin `lf_maxbin[f]`. Default is n_len (no aggregation).}
#' }
#'
#' @format A list with 25 named elements combining dimensions, biological parameters,
#' observed data (catch and CPUE), and Bayesian priors for RTMB model initialization.
#'
#' @source Prepared from MFCL (Multifan-CL) model outputs and frequency files,
#'   with biological parameters and priors derived from assessment modeling.
#'
#' @examples
#' \dontrun{
#'   data(wcpo_bet_data)
#'   str(wcpo_bet_data)
#'   # Access catch data
#'   head(wcpo_bet_data$cpue_data)
#'   # Get biological parameters
#'   wcpo_bet_data$M
#' }
#'
"wcpo_bet_data"

#' West Central Pacific Ocean Bigeye Tuna (WCPO BET) Assessment Parameters
#'
#' A list of estimated parameter values for the WCPO bigeye tuna stock assessment
#' model. Contains growth parameters, recruitment deviations, selectivity curves,
#' and other model parameters fitted using RTMB.
#'
#' @details
#' ## Parameter structure (14 list elements)
#'
#' **Stock parameters:**
#' \describe{
#'   \item{log_B0}{Initial spawning biomass (log-scale): 12}
#'   \item{log_h}{Steepness parameter (log-scale): -0.0513}
#'   \item{log_sigma_r}{Recruitment standard deviation (log-scale): -0.511}
#' }
#'
#' **Observation model parameters:**
#' \describe{
#'   \item{log_cpue_q}{CPUE catchability coefficient (log-scale): 0}
#'   \item{cpue_creep}{CPUE creep adjustment: 0}
#'   \item{log_cpue_tau}{CPUE log-standard deviation (log-scale): -Inf}
#'   \item{log_cpue_omega}{CPUE process error (log-scale): 0}
#' }
#'
#' **Recruitment deviations:**
#' \describe{
#'   \item{rdev_y}{Recruitment deviations by year (268 elements, quarterly timesteps)}
#' }
#'
#' **Selectivity parameters:**
#' \describe{
#'   \item{par_sel}{Selectivity parameters matrix (15 fisheries × 6 columns).
#'     Rows represent different fisheries/fleets, columns represent selectivity
#'     curve parameters (e.g., inflection point, slope for logistic curves or
#'     parameters for double-normal curves).}
#' }
#'
#' **Growth parameters:**
#' \describe{
#'   \item{log_L1}{Log of length at age 1 (von Bertalanffy): 3.431379}
#'   \item{log_L2}{Log of asymptotic length (von Bertalanffy): 5.03333}
#'   \item{log_k}{Log of growth rate parameter (von Bertalanffy): -2.32024}
#'   \item{log_CV1}{Log of CV at young ages: -1.826289}
#'   \item{log_CV2}{Log of CV at old ages: -2.230151}
#' }
#'
#' @format A list with 14 named elements containing estimated parameters from
#' model optimization.
#'
#' @seealso [wcpo_bet_data] for the corresponding data object used in model fitting.
#'
#' @examples
#' \dontrun{
#'   data(wcpo_bet_parameters)
#'   str(wcpo_bet_parameters)
#'   # Access growth parameters
#'   exp(wcpo_bet_parameters$log_L1)  # L1 on real scale
#'   exp(wcpo_bet_parameters$log_L2)  # Linf on real scale
#'   exp(wcpo_bet_parameters$log_k)   # k on real scale
#' }
#'
"wcpo_bet_parameters"

#' West Central Pacific Ocean Bigeye Tuna (WCPO BET) Length Frequency Data
#'
#' Length frequency composition data extracted from MFCL (Multifan-CL) frequency files
#' and prepared for RTMB (R Template Model Builder) stock assessment modeling.
#' Observations are organized by fishery, year, month, and 2 cm length bins.
#'
#' @details
#' ## Data workflow
#'
#' **Preparation steps:**
#' \enumerate{
#'   \item Extract length composition data from MFCL frequency file using `lnfrq()`
#'     - Returns wide-format data with length bins as columns (10, 12, 14, ... 198 cm)
#'     - Observations grouped by fishery-year-month
#'   \item Reshape to long format (one row per fishery-year-month-length bin combination)
#'   \item Merge with timestep table to add model timestep indexing consistent with catch/CPUE
#'   \item Save as long-format data for RTMB model input
#' }
#'
#' ## Output columns
#'
#' \describe{
#'   \item{year}{Calendar year (1952-2018). Integer.}
#'   \item{month}{Fishing quarter month (2, 5, 8, 11 representing Q1-Q4 approx.). Numeric.}
#'   \item{ts}{Model timestep number (1-268 for 67 years × 4 quarters). Integer.}
#'   \item{fishery}{MFCL fleet index (1-15). Numeric.}
#'   \item{bin}{Length bin lower edge (cm). 2 cm bins: 10, 12, 14, ... 198. Numeric.}
#'   \item{value}{Observed count of fish in this length bin. Numeric.}
#'   \item{week}{Temporal indicator within quarter (typically 1 for aggregated data). Numeric.}
#' }
#'
#' ## Notes
#'
#' - MFCL data uses 2 cm length bins (10, 12, 14, ... 198 cm)
#' - Data from `lnfrq()` represents observed length frequency
#' - Sum of values across bins per fishery-year-month instance gives total input sample size
#' - Bin values are left as counts for RTMB to handle (divide by sample size if proportions needed)
#'
#' @format A data.table with 7 columns:
#'   year (integer), month (numeric), ts (integer), fishery (numeric),
#'   bin (numeric), value (numeric), week (numeric).
#'   Multiple rows per fishery-year-month (one per 2 cm length bin).
#'
#' @source Extracted from MFCL (Multifan-CL) model frequency file outputs
#'   and prepared for RTMB model input.
#'
#' @examples
#' \dontrun{
#'   data(wcpo_bet_lf)
#'   str(wcpo_bet_lf)
#'   # View sample data
#'   head(wcpo_bet_lf)
#'   # Get total sample size per fishery-year-month
#'   wcpo_bet_lf[, .(total_n = sum(value)), by = .(year, month, fishery)]
#' }
#'
"wcpo_bet_lf"

#' West Central Pacific Ocean Bigeye Tuna (WCPO BET) Weight Frequency Data
#'
#' Weight frequency composition data extracted from MFCL (Multifan-CL) frequency files
#' and prepared for RTMB (R Template Model Builder) stock assessment modeling.
#' Observations are organized by fishery, year, month, and 1 kg weight bins.
#'
#' @details
#' ## Data workflow
#'
#' **Preparation steps:**
#' \enumerate{
#'   \item Extract weight composition data from MFCL frequency file using `wtfrq()`
#'     - Returns wide-format data with weight bins as columns (1, 2, 3, ... 200 kg)
#'     - Observations grouped by fishery-year-month
#'   \item Reshape to long format (one row per fishery-year-month-weight bin combination)
#'   \item Merge with timestep table to add model timestep indexing consistent with catch/CPUE
#'   \item Save as long-format data for RTMB model input
#' }
#'
#' ## Output columns
#'
#' \describe{
#'   \item{year}{Calendar year (1952-2018). Integer.}
#'   \item{month}{Fishing quarter month (2, 5, 8, 11 representing Q1-Q4 approx.). Numeric.}
#'   \item{ts}{Model timestep number (1-268 for 67 years × 4 quarters). Integer.}
#'   \item{fishery}{MFCL fleet index (1-15). Numeric.}
#'   \item{bin}{Weight bin (kg). 1 kg bins: 1, 2, 3, ... 200. Integer.}
#'   \item{value}{Observations in this weight bin. Numeric.}
#'   \item{week}{Temporal indicator within quarter (typically 1 for aggregated data). Numeric.}
#' }
#'
#' ## Notes
#'
#' - MFCL data uses 1 kg weight bins (1, 2, 3, ... 200 kg)
#' - Data from `wtfrq()` represents observed weight frequency
#' - Sum of values across bins per fishery-year-month instance gives total input sample size
#' - Bin values are left as counts for RTMB to handle (divide by sample size if proportions needed)
#'
#' @format A data.table with 7 columns:
#'   year (integer), month (numeric), ts (integer), fishery (numeric),
#'   bin (integer), value (numeric), week (numeric).
#'   Multiple rows per fishery-year-month (one per 1 kg weight bin).
#'
#' @source Extracted from MFCL (Multifan-CL) model frequency file outputs
#'   and prepared for RTMB model input.
#'
#' @examples
#' \dontrun{
#'   data(wcpo_bet_wf)
#'   str(wcpo_bet_wf)
#'   # View sample data
#'   head(wcpo_bet_wf)
#'   # Get total sample size per fishery-year-month
#'   wcpo_bet_wf[, .(total_n = sum(value)), by = .(year, month, fishery)]
#' }
#'
"wcpo_bet_wf"

#' `opal` model regression baseline
#'
#' A saved reference run for the WCPO bigeye tuna assessment that documents the
#' objective, gradient, and RTMB report output used by the
#' `vignettes/baseline.Rmd` workflow to detect numerical regressions and timing
#' changes before and after refactors.
#'
#' @details
#' ## Baseline contents
#' \describe{
#'   \item{timestamp}{Time when the baseline was created.}
#'   \item{description}{Free-form label ("Pre-refactoring baseline" by default).}
#'   \item{opal_version}{Version string of the installed `opal` package.}
#'   \item{dimensions}{Named list with `n_year`, `n_age`, `n_fishery`, `n_season`,
#'     `n_len`, `n_lf`, `n_wt`, and `n_wf` from the model inputs.}
#'   \item{nll}{Negative log-likelihood at the reference parameter set.}
#'   \item{gradient}{Gradient vector evaluated at `obj$par` for the baseline.}
#'   \item{max_gr}{Maximum absolute gradient entry.}
#'   \item{timing}{List with `fn_elapsed`, `gr_elapsed`, `fn_median`, and
#'     `gr_median` values that summarize runtime measurements used in the
#'     comparison vignette.}
#'   \item{report}{Sublist containing the baseline's `number_ysa`,
#'     `spawning_biomass_y`, catch and harvest rate predictions (`catch_pred_fya`,
#'     `catch_pred_ysf`, `hrate_ysa`, `hrate_ysfa`), selectivity curves
#'     (`sel_fya`), and log-likelihood components (`lp_prior`, `lp_penalty`,
#'     `lp_rec`, `lp_cpue`, `lp_lf`, `lp_wf`).}
#'   \item{par_values}{Vector of `obj$par` values that produced the baseline.}
#' }
#'
#' @format A named list accessible via `data(opal_baseline)` and
#'   packaged as `data/opal_baseline.rda`.
#'
#' @source `vignettes/baseline.Rmd` (the "Save or compare" chunk) which
#'   reproduces the reference run saved to `data/opal_baseline.rda`.
#'
#' @examples
#' \dontrun{
#'   data(opal_baseline)
#'   str(opal_baseline$dimensions)
#'   opal_baseline$timing$fn_median
#' }
#'
"opal_baseline"

#' Opal Baseline Data
#'
#' The full data list used for testing in the opal baseline model.
#'
#' @source \code{vignettes/baseline.Rmd}
"opal_baseline_data"

#' Opal Baseline Parameters
#'
#' The full parameter list used for testing in the opal baseline model.
#'
#' @source \code{vignettes/baseline.Rmd}
"opal_baseline_parameters"

#' Opal Baseline Map
#'
#' The parameter map used for testing in the opal baseline model.
#'
#' @source \code{vignettes/baseline.Rmd}
"opal_baseline_map"

#' 'Opakapaka Stock Assessment Data
#'
#' A comprehensive list containing biological parameters, catch and CPUE data,
#' and Bayesian prior specifications prepared for RTMB (R Template Model Builder)
#' stock assessment modeling. Extracted from a Stock Synthesis 3 (SS3) model
#' and converted to opal data conventions.
#'
#' @details
#' ## Data structure (22 list elements)
#'
#' **Age and time dimensions:**
#' \describe{
#'   \item{age_a}{Vector of integer ages (1:44). Length = 44. Internal opal ages
#'     mapping to SS3 ages 0-43 (opal age \code{i} = SS3 age \code{i - 1}).}
#'   \item{n_age}{Number of age classes (44).}
#'   \item{n_season}{Number of seasons per year (1).}
#'   \item{n_fishery}{Number of fisheries/fleets (3).}
#'   \item{n_year}{Number of years in model (75, representing 1949-2023).}
#'   \item{first_yr}{First year designation in model (1, representing 1949).}
#'   \item{last_yr}{Last year designation in model (75, representing 2023).}
#'   \item{years}{Vector of model years (1:75).}
#' }
#'
#' **Length structure:**
#' \describe{
#'   \item{len_bin_start}{First length bin lower edge (5 cm).}
#'   \item{len_bin_width}{Width of each length bin (5 cm).}
#'   \item{n_len}{Total number of length bins (17).}
#' }
#'
#' **Catch data:**
#' \describe{
#'   \item{first_yr_catch}{First year with catch data (1).}
#'   \item{catch_units_f}{Numeric vector of unit codes by fishery (3 elements).
#'     All fisheries use units = 1 (biomass).}
#'   \item{catch_obs_ysf}{Three-dimensional array of observed catch
#'     (75 years \times 1 season \times 3 fisheries). Dimensions follow
#'     (year, season, fishery). Fleet 3 (ResFish survey) has zero catch.}
#' }
#'
#' **CPUE data:**
#' \describe{
#'   \item{cpue_switch}{Indicator for CPUE inclusion (1 = included).}
#'   \item{cpue_data}{Data.table with 82 rows and 8 columns containing CPUE
#'     observations from two fleets:
#'     \itemize{
#'       \item{year}: Calendar year (1949-2023)
#'       \item{month}: Observation month (7 for Comm fleet, 1 for ResFish)
#'       \item{ts}: Sequential timestep (1-75)
#'       \item{fishery}: Fleet identifier (1 = Comm, 3 = ResFish)
#'       \item{index}: Integer index identifier (1 for Comm, 2 for ResFish)
#'       \item{metric}: Data type ("cpue")
#'       \item{units}: Unit code (1 for biomass)
#'       \item{value}: CPUE observations
#'       \item{se}: Log-space standard error (0.2 for Comm, 0.15 for ResFish)
#'     }
#'     Comm fleet (fleet 1) covers 1949-2023 (75 obs). ResFish survey (fleet 3)
#'     covers 2017-2023 (7 obs).
#'   }
#' }
#'
#' **Biological parameters:**
#' \describe{
#'   \item{lw_a}{Length-weight allometric coefficient (1.75e-05).}
#'   \item{lw_b}{Length-weight allometric exponent (2.99).}
#'   \item{maturity}{Maturity at length vector (17 elements). Logistic maturity
#'     probabilities at each data length bin midpoint (Mat50 = 40.7 cm,
#'     slope = -2.26).}
#'   \item{fecundity}{Fecundity at length vector (17 elements). Fecundity equals
#'     weight at each length bin midpoint (fecundity option: eggs = Wt).}
#'   \item{M}{Natural mortality at age vector (44 elements). Constant
#'     instantaneous rate (0.135) across all ages.}
#' }
#'
#' **Growth reference ages (Schnute parameterization):**
#' \describe{
#'   \item{A1}{Reference age for L1 (1). Corresponds to SS3 age 0.}
#'   \item{A2}{Reference age for L2 (44 = n_age). Corresponds to SS3 age 43.}
#' }
#'
#' **Selectivity:**
#' \describe{
#'   \item{sel_type_f}{Selectivity function type by fishery (3 elements). Integer
#'     codes: 1 = logistic (fleets 1-2), 24 = double-normal (fleet 3).}
#' }
#'
#' **Priors:**
#' \describe{
#'   \item{priors}{List of 8 Bayesian prior specifications, each containing:
#'     \itemize{
#'       \item{type}: Prior distribution type ("none" = uninformative)
#'       \item{par1}: First parameter (NA when type = "none")
#'       \item{par2}: Second parameter (NA when type = "none")
#'       \item{index}: Parameter index in the model
#'     }
#'     Priors cover: log_B0, log_cpue_q, par_sel, log_L1, log_L2, log_k,
#'     log_CV1, log_CV2.
#'   }
#' }
#'
#' ## Fleet definitions
#'
#' \describe{
#'   \item{Fleet 1 (Comm)}{Commercial fishery. Catch fleet with logistic
#'     selectivity. CPUE and length composition data available 1949-2023.}
#'   \item{Fleet 2 (Non_comm)}{Non-commercial fishery. Catch fleet with logistic
#'     selectivity (fixed parameters). Catch data only, no CPUE or compositions.}
#'   \item{Fleet 3 (ResFish)}{Research fishery survey. Survey fleet with
#'     double-normal selectivity. CPUE and length composition data available
#'     2017-2023.}
#' }
#'
#' ## Growth conversion
#'
#' Growth parameters were converted from the SS3 standard von Bertalanffy
#' parameterization (L_at_Amin = 6.0 at age 0, Linf = 67.5 at age 999,
#' k = 0.242) to the opal Schnute parameterization with A1 = 1 and A2 = n_age.
#' The conversion is exact to machine precision (~1e-14). The VB growth
#' coefficient k is unchanged between parameterizations.
#'
#' @format A list with 22 named elements combining dimensions, biological
#' parameters, observed data (catch and CPUE), and Bayesian priors for RTMB
#' model initialization.
#'
#' @source Extracted from Stock Synthesis 3 (SS3 v3.30.19) model input files
#'   (\code{data_echo.ss_new}, \code{control.ss_new}). Historical data only
#'   (1949-2023), forecast period excluded.
#'
#' @seealso [opaka_parameters] for starting parameter values,
#'   [opaka_lf] for length frequency composition data.
#'
#' @examples
#' \dontrun{
#'   data(opaka_data)
#'   str(opaka_data)
#'   # Access catch data
#'   head(opaka_data$cpue_data)
#'   # Get biological parameters
#'   opaka_data$M
#'   # Check dimensions
#'   cat("Ages:", opaka_data$n_age, "\n")
#'   cat("Years:", opaka_data$n_year, "\n")
#'   cat("Fleets:", opaka_data$n_fishery, "\n")
#'   cat("Length bins:", opaka_data$n_len, "\n")
#' }
#'
"opaka_data"

#' 'Opakapaka Stock Assessment Parameters
#'
#' A list of parameter values for the opakapaka stock assessment model,
#' extracted from a fitted Stock Synthesis 3 (SS3) model and converted to
#' opal conventions. Contains growth parameters, recruitment deviations,
#' selectivity curves, and other model parameters suitable as starting
#' values for RTMB optimization.
#'
#' @details
#' ## Parameter structure (14 list elements)
#'
#' **Stock parameters:**
#' \describe{
#'   \item{log_B0}{Log initial recruitment/biomass (SS3's LN(R0)): 5.656}
#'   \item{log_h}{Steepness parameter (log-scale): log(0.76) = -0.274}
#'   \item{log_sigma_r}{Recruitment standard deviation (log-scale):
#'     log(0.52) = -0.654}
#' }
#'
#' **Observation model parameters:**
#' \describe{
#'   \item{log_cpue_q}{CPUE catchability coefficients (log-scale). Vector of
#'     length 2: -3.773 (Comm fleet), -6.253 (ResFish survey).}
#'   \item{cpue_creep}{CPUE creep adjustment: 0}
#'   \item{log_cpue_tau}{CPUE observation error (log-scale): -Inf (disabled)}
#'   \item{log_cpue_omega}{CPUE process error (log-scale): 0}
#' }
#'
#' **Recruitment deviations:**
#' \describe{
#'   \item{rdev_y}{Recruitment deviations by year (75 elements, annual timesteps
#'     1949-2023). Extracted from SS3 estimated main recruitment deviations.}
#' }
#'
#' **Selectivity parameters:**
#' \describe{
#'   \item{par_sel}{Selectivity parameters matrix (3 fisheries \times 6 columns).
#'     \itemize{
#'       \item{Fleet 1 (Comm)}: Logistic. Cols 1-2 = inflection (36.1 cm),
#'         95\% width (4.08 cm). Cols 3-6 unused.
#'       \item{Fleet 2 (Non_comm)}: Logistic (fixed). Cols 1-2 = inflection
#'         (40 cm), 95\% width (11 cm). Cols 3-6 unused.
#'       \item{Fleet 3 (ResFish)}: Double-normal (SS3 pattern 24). Cols 1-6 =
#'         peak (20.5), top_logit (-1.06), ascend_se (-0.34), descend_se (4.11),
#'         start_logit (-999), end_logit (-1.40).
#'     }
#'   }
#' }
#'
#' **Growth parameters (Schnute parameterization, log-scale):**
#' \describe{
#'   \item{log_L1}{Log of length at reference age A1 = 1 (opal internal),
#'     corresponding to SS3 age 0: log(6.0) = 1.792}
#'   \item{log_L2}{Log of length at reference age A2 = 44 (opal internal),
#'     corresponding to SS3 age 43: log(67.498) = 4.212}
#'   \item{log_k}{Log of von Bertalanffy growth coefficient:
#'     log(0.242) = -1.419}
#'   \item{log_CV1}{Log of CV at young ages: log(0.085) = -2.465}
#'   \item{log_CV2}{Log of CV at old ages: log(0.085) = -2.465. Equal to
#'     log_CV1, giving constant CV across all lengths.}
#' }
#'
#' @format A list with 14 named elements containing parameters extracted from
#' SS3 model optimization and converted to opal conventions.
#'
#' @source Extracted from Stock Synthesis 3 parameter file (\code{ss.par})
#'   and control file (\code{control.ss_new}). Growth parameters converted
#'   from standard VB to Schnute parameterization with A1 = 1, A2 = n_age.
#'
#' @seealso [opaka_data] for the corresponding data object,
#'   [opaka_lf] for length frequency composition data.
#'
#' @examples
#' \dontrun{
#'   data(opaka_parameters)
#'   str(opaka_parameters)
#'   # Access growth parameters on real scale
#'   exp(opaka_parameters$log_L1)   # L1 = 6.0 cm
#'   exp(opaka_parameters$log_L2)   # L2 = 67.5 cm
#'   exp(opaka_parameters$log_k)    # k = 0.242
#'   # Steepness
#'   exp(opaka_parameters$log_h)    # h = 0.76
#' }
#'
"opaka_parameters"

#' 'Opakapaka Length Frequency Data
#'
#' Length frequency composition data extracted from a Stock Synthesis 3 (SS3)
#' model and prepared for RTMB (R Template Model Builder) stock assessment
#' modeling. Observations are organized by fishery, year, month, and 5 cm
#' length bins.
#'
#' @details
#' ## Data workflow
#'
#' **Preparation steps:**
#' \enumerate{
#'   \item Parse length composition observations from SS3 data file
#'     (\code{data_echo.ss_new})
#'   \item Extract proportions and sample sizes by fleet-year
#'   \item Convert proportions to pseudo-counts (proportion \times Nsamp)
#'   \item Reshape to long format (one row per fishery-year-length bin)
#'   \item Truncate to historical period (1949-2023)
#' }
#'
#' ## Output columns
#'
#' \describe{
#'   \item{year}{Calendar year (1949-2023). Integer.}
#'   \item{month}{Observation month (1 for all observations). Integer.}
#'   \item{ts}{Model timestep number (1-75). Integer.}
#'   \item{fishery}{Fleet index: 1 = Comm (75 obs), 3 = ResFish (7 obs).
#'     Numeric.}
#'   \item{bin}{Length bin lower edge (cm). 5 cm bins: 5, 10, 15, ..., 85.
#'     Numeric.}
#'   \item{value}{Pseudo-count (proportion \times Nsamp) in this length bin.
#'     Numeric.}
#'   \item{week}{Temporal indicator (always 1). Numeric.}
#' }
#'
#' ## Fleet coverage
#'
#' \describe{
#'   \item{Fleet 1 (Comm)}{75 observation-years (1949-2023). Sample sizes
#'     range from ~23 to ~755.}
#'   \item{Fleet 2 (Non_comm)}{No length composition data available.}
#'   \item{Fleet 3 (ResFish)}{7 observation-years (2017-2023). Fixed sample
#'     size of 60 per year.}
#' }
#'
#' ## Notes
#'
#' \itemize{
#'   \item Data length bins are 5 cm wide (5, 10, 15, ..., 85 cm), giving 17 bins
#'   \item Values are pseudo-counts derived from SS3 proportions \times Nsamp
#'   \item Sum of values across bins per fishery-year gives the effective sample size
#'   \item Very small pseudo-counts (< 1e-10) represent SS3's numerical zeros
#'     and should be treated as absent observations
#' }
#'
#' @format A data.table with 7 columns:
#'   year (integer), month (integer), ts (integer), fishery (numeric),
#'   bin (numeric), value (numeric), week (numeric).
#'   Multiple rows per fishery-year (one per 5 cm length bin). Total rows = 1394
#'   (82 observation-years \times 17 bins).
#'
#' @source Extracted from Stock Synthesis 3 data file (\code{data_echo.ss_new}).
#'   Historical period only (1949-2023).
#'
#' @seealso [opaka_data] for the main data object,
#'   [opaka_parameters] for starting parameter values.
#'
#' @examples
#' \dontrun{
#'   data(opaka_lf)
#'   str(opaka_lf)
#'   # View sample data
#'   head(opaka_lf)
#'   # Get effective sample size per fishery-year
#'   opaka_lf[, .(total_n = sum(value)), by = .(year, fishery)]
#'   # Count observations by fleet
#'   opaka_lf[, uniqueN(paste(year, month)), by = fishery]
#' }
#'
"opaka_lf"

#' Opakapaka SS3 OM/EM truth and EM output (extracted)
#'
#' A small list containing Stock Synthesis (SS3) output extracted from an
#' operating model (OM) and an estimation model (EM) pair using
#' `r4ss::SS_output()`. Intended as a truth / comparison object for
#' simulation-experiment workflows and diagnostics.
#'
#' @details
#' The object is a list with two elements: `om` and `em`. Each element is a
#' sub-list that contains `sprseries` (spawning biomass time series) and `cpue`
#' (standardized index output). The contained data frames have the following
#' structures:
#'
#' \describe{
#'   \item{om$sprseries}{A data.frame with 100 rows and 2 columns:
#'     \describe{
#'       \item{Yr}{Integer year (e.g., 1949...)}
#'       \item{SSB}{Numeric spawning stock biomass}
#'     }
#'   }
#'   \item{om$cpue}{A data.frame with 107 rows and 7 columns:
#'     \describe{
#'       \item{Yr}{Integer year}
#'       \item{Fleet}{Integer fleet id}
#'       \item{Fleet_name}{Character fleet label}
#'       \item{Vuln_bio}{Numeric vulnerable biomass predicted by SS3}
#'       \item{Obs}{Observed CPUE (integer in OM)}
#'       \item{Exp}{Expected CPUE (numeric)}
#'       \item{SE}{Standard error (numeric)}
#'     }
#'   }
#'   \item{em$sprseries}{A data.frame with 100 rows and 2 columns: same format as `om$sprseries`.}
#'   \item{em$cpue}{A data.frame with 107 rows and 7 columns: same column names as `om$cpue`,
#'     with `Obs` shown as numeric in the EM output shown here.}
#' }
#'
#' @format A list with two named elements, `om` and `em`, each containing the
#'   SS3 `sprseries` and `cpue` outputs as data.frames suitable for plotting
#'   and comparison of true vs estimated indices and spawning biomass trajectories.
#'
#' @source Produced by calling `r4ss::SS_output()` on OM and EM model directories
#'   and collecting the `sprseries` and `cpue` components.
#'
#' @examples
#' \dontrun{
#'   data(opaka_truth)
#'   str(opaka_truth)
#'   plot(opaka_truth$om$sprseries$Yr, opaka_truth$om$sprseries$SSB, type = "l")
#'   lines(opaka_truth$em$sprseries$Yr, opaka_truth$em$sprseries$SSB, col = "red")
#'   # Compare CPUE
#'   head(opaka_truth$om$cpue)
#'   head(opaka_truth$em$cpue)
#' }

"opaka_truth"