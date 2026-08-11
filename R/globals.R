#' @importFrom RTMB MakeADFun
#' @importFrom utils data head
NULL

parse_number <- function(x) {
  suppressWarnings(as.numeric(gsub("[^0-9.+-]+", "", as.character(x))))
}

get_sliced_afs <- function(...) {
  stop("get_sliced_afs() is not currently implemented in opal.", call. = FALSE)
}

get_selectivity_v1 <- function(...) {
  stop("get_selectivity_v1() is not currently implemented in opal.", call. = FALSE)
}

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    ".",
    "0", "1", "2", "3", "4", "5", "6", "30",
    "aerial_cov", "aerial_survey", "age", "Age", "age_freq",
    "Australian", "Bin", "CaptureCov", "CaptureSwitch", "CaptureYear",
    "catch", "catch_units_f", "catch_UA", "cdiff", "chain", "change",
    "cmax", "cmin", "Cohort", "cohort1", "cohort2", "Comps",
    "cpue", "data_csv1", "data_par1", "dBin", "fishery", "Fishery",
    "Group", "GTs", "HSPs", "id", "ifishery", "index", "iter",
    "iteration", "L1", "L2", "len_lower", "len_mid", "len_upper",
    "Length", "length_freq", "length_mean", "length_sd", "lf_fishery_f",
    "lf_n_f", "lf_obs_flat", "lf_obs_ints", "lf_obs_prop", "LL1",
    "LL1_case", "log_lf_tau", "lp__", "MaxAge", "MinAge", "Model",
    "N", "n_len", "name", "nC", "nK", "NPOPS", "obj", "output",
    "paly", "P_t_20", "par", "parameter", "pars", "POPs", "RecAge",
    "RecYear", "RelAge", "RelYear", "removal", "SD", "season",
    "Season", "sim", "spawning_potential", "Surf", "surf_case",
    "tag_recaptures", "tag_releases", "tag_reporting", "total",
    "troll", "UA", "V1", "V2", "value", "Var1", "Var2", "Var3",
    "weight", "year", "Year"
  ))
}
