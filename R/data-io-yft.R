#' Indian Ocean yellowfin tuna (simulated) assessment data
#'
#' Model inputs for a quarterly, single-region assessment of simulated Indian
#' Ocean yellowfin tuna, converted from the IATTC analyst team's Stock
#' Synthesis 3 (SS3) model for the Spatial Stock Assessment Simulation
#' Experiment (replicate 4). Used by the "Comparing opal with Stock Synthesis"
#' article.
#'
#' @details
#' The model runs 256 quarterly time steps (`n_year = 256`, `n_season = 1`)
#' and SS3 ages 1--28, so the opal age index equals the SS3 age. There are 20
#' fisheries: 16 extraction fleets with catch in numbers, and the longline
#' CPUE index split into four quarterly series (fisheries 17--20), each with
#' its own catchability. Length compositions (5 cm bins) are multinomial for
#' fisheries 1--14 and 16, with SS3's `mult_by_lencomp_N` carried as
#' `lf_var_adjust`; [prep_lf_data()] has already been applied.
#'
#' Biology is SS3's age-based end-of-run biology: M, weight, maturity and
#' spawning potential at age, with M scaled to the quarter. Growth uses SS3's
#' `GrowthModel 3` (age-specific K): `k_mult_ages` gives the ages with K
#' multipliers, which are in `io_yft_parameters$log_k_mult`, and
#' `growth_plus_group_adjust` applies SS3's plus-group length approximation.
#' Selectivity is length-based (`sel_type_f`: 2 = double normal for fisheries
#' 1--16, 1 = logistic for the index series). `sel_basis_f` is not set, so
#' every curve is evaluated on length; adding `sel_basis_f = rep(2L, 20)`
#' evaluates them on age instead, as SS3 does (see the "Length- or age-based
#' selectivity" article). `priors` holds the default
#' priors from [get_priors()].
#'
#' @format A named list of 62 elements in the form [opal_obj()] expects:
#'   dimensions, catch (`catch_obs_ysf`), the index (`cpue_data`), biology,
#'   growth settings, length compositions (`lf_*`), the recruitment bias
#'   adjustment (`bias_adj_y`) and priors.
#' @source Converted from SS3 v3.24 files with the SS3-to-opal converter in
#'   <https://github.com/N-DucharmeBarth/wcpfc-bet-single-region>
#'   (`code/rtmb/scripts/04-export-io-yft-for-opal.r`). SS3 model by M.N.
#'   Maunder, H. Xu and C. Minte-Vera (IATTC); data simulated from the
#'   operating model described in Goethel et al. (2024)
#'   \doi{10.1111/faf.12819}.
#' @seealso [io_yft_parameters], [io_yft_map], [io_yft_reference]
#' @examples
#' data(io_yft_data)
#' io_yft_data$n_year
#' head(io_yft_data$cpue_data)
"io_yft_data"

#' Indian Ocean yellowfin tuna (simulated) starting parameters
#'
#' Starting values for [io_yft_data], at the SS3 solution: SS3's unfished
#' spawning biomass (`log_B0`), recruitment deviations (`rdev_y`), steepness,
#' sigma_R and growth (including the age-specific K multipliers `log_k_mult`).
#' Selectivity (`par_sel`) starts from opal curves fitted by least squares to
#' SS3's selectivity-at-age projected onto length.
#'
#' @format A named list of 18 parameter vectors and matrices, including
#'   `log_B0`, `log_cpue_q` (4), `rdev_y` (256), `par_sel` (20 x 6) and
#'   `log_k_mult` (9).
#' @source As for [io_yft_data].
#' @seealso [io_yft_map]
"io_yft_parameters"

#' Indian Ocean yellowfin tuna (simulated) parameter map
#'
#' The map for fit "C2" of [io_yft_data]: `log_B0`, the 256 recruitment
#' deviations, the four quarterly catchabilities and 49 selectivity
#' parameters are estimated (310 in total). Growth, natural mortality,
#' steepness and sigma_R are fixed at SS3's values.
#'
#' Selectivity estimates the opal equivalents of the parameters SS3
#' estimated: columns 1--4 of the double normal for fisheries 1--14 and both
#' logistic parameters for the index; fisheries 15--16 mirror 14, and the
#' index series 18--20 mirror 17. Nine parameters the length data do not
#' inform (the plateau and descending width of fisheries 4, 5, 7 and 13, and
#' the descending width of fishery 9) are fixed at their starting values.
#'
#' @format A named list of factors, as for the `map` argument of [opal_obj()].
#' @source As for [io_yft_data].
"io_yft_map"

#' Indian Ocean yellowfin tuna (simulated) SS3 and truth reference values
#'
#' Quantities from the SS3 model and from the simulation operating model,
#' for comparison with opal fits to [io_yft_data].
#'
#' @format A named list:
#' \describe{
#'   \item{fleet_names}{Names of the 20 opal fisheries.}
#'   \item{ss3}{SS3 results: `version`; unfished spawning biomass `ssb0` (t);
#'     `trajectory` (`ts`, `ssb`, `static_dep`); `recdevs` (`ts`, `recdev`);
#'     `cpue` (`ts`, `fishery`, `obs`, `expected`, `se`), each quarterly
#'     series scaled to mean one; `lf` (`fishery`, `bin`, `expected`),
#'     expected length compositions aggregated per fishery weighted by sample
#'     size; `sel_age` (`fishery`, `age`, `sel`); `sel_pattern` (`fleet`,
#'     `name`, `pattern`, `special`), SS3's `age_selex` pattern per fleet
#'     (20 = double normal, 12 = logistic, 15 = mirror of fleet `special`);
#'     `sel_par` (`fleet`, `pattern`, `par`, `value`, `phase`, `lo`, `hi`,
#'     `label`), SS3's estimated age selectivity parameters P1--P6, with
#'     `phase > 0` for the ones SS3 estimated (see
#'     [convert_ss3_selex_to_rtmb()] with `basis = "age"`); and `length_at_age`
#'     (`age`, `length`, start of the time step).}
#'   \item{truth}{The operating model (SPM, 159 cells on a 5 degree grid),
#'     summed over the whole domain: `b0`, unfished equilibrium spawning
#'     biomass (t), and `trajectory` (`ts`, `ssb`, `static_dep`).}
#'   \item{schnute_length_at_age}{Length-at-age (`age`, `length`) of the
#'     single-k Schnute curve used before opal supported age-specific K.}
#'   \item{history}{B0 and terminal depletion of three opal fits: Schnute
#'     growth with SS3's selectivity, SS3 growth with SS3's selectivity, and
#'     SS3 growth with opal's estimated selectivity (C2).}
#'   \item{source}{Provenance of the export: repository, script, commit,
#'     date and package versions.}
#' }
#' @source SS3 Report.sso and data.ss_new expected values; operating model
#'   outputs from
#'   <https://github.com/aaronmberger-nwfsc/Spatial-Assessment-Modeling-Workshop>.
#'   Exported by `code/rtmb/scripts/04-export-io-yft-for-opal.r` in
#'   <https://github.com/N-DucharmeBarth/wcpfc-bet-single-region>.
"io_yft_reference"
