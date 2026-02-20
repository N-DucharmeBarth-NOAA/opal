#' West Central Pacific Ocean Bigeye Tuna (WCPO BET) Assessment Data
#'
#' A comprehensive list containing biological parameters, catch and CPUE data,
#' and Bayesian prior specifications prepared for RTMB (R Template Model Builder)
#' stock assessment modeling of western central Pacific bigeye tuna.
#'
#' @details
#' ## Data structure (25 list elements)
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
#'   \item{cpue_data}{Tibble with 268 rows and 8 columns containing standardized CPUE observations:
#'     \itemize{
#'       \item{year}: Calendar year (1952-2018)
#'       \item{month}: Quarter month (2, 5, 8, 11 for Q1-Q4)
#'       \item{ts}: Sequential timestep (1-268)
#'       \item{fishery}: Fleet identifier (usually 15 for survey fleet)
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
