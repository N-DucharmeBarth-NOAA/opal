#' Mean length at age
#' 
#' Mean length (cm) at age by year and season.
#'
#' @format a \code{tibble} containing 33 fields including:
#' \describe{
#'   \item{Year}{The year.}
#'   \item{Season}{The season.}
#'   \item{0-30}{The mean length for each age for each year and season.}
#'  }
#' 
"length_mean"

#' Standard deviation of mean length at age
#' 
#' Standard deviation (SD) of mean length (cm) at age.
#'
#' @format a \code{tibble} containing 2 fields including:
#' \describe{
#'   \item{Age}{Age in years.}
#'   \item{SD}{Standard deviation of length at age.}
#'  }
#' 
"length_sd"

#' Catch
#' 
#' Standard deviation of mean length (cm) by age.
#'
#' @format a \code{tibble} containing 2 fields including:
#' \describe{
#'   \item{Age}{The year.}
#'   \item{SD}{The unscaled index of abundance.}
#'  }
#' 
"catch"

#' Catch UA
#' 
#' Catch UA.
#'
#' @format a \code{tibble} containing 2 fields including:
#' \describe{
#'   \item{Age}{The year.}
#'   \item{SD}{The unscaled index of abundance.}
#'  }
#' 
"catch_UA"

#' Parent offspring pairs
#' 
#' Parent offspring pairs (POPs).
#'
#' @format a \code{tibble} containing 2 fields including:
#' \describe{
#'   \item{Age}{The year.}
#'   \item{SD}{The unscaled index of abundance.}
#'  }
#' 
"POPs"

#' Parent offspring pairs
#' 
#' Parent offspring pairs (POPs).
#'
#' @format a \code{tibble} containing 2 fields including:
#' \describe{
#'   \item{Age}{The year.}
#'   \item{SD}{The unscaled index of abundance.}
#'  }
#' 
"POPs_v1"

#' Distribution of age given length in genotyped adults
#' 
#' paly(a,l,y) distribution of age given length in genotyped adults if not aged.
#' 
#' @format a \code{tibble} containing 33 fields including:
#' 
"paly"

#' Half sibling pairs
#' 
#' Half sibling pairs (HSPs).
#'
#' @format a \code{tibble} containing 2 fields including:
#' \describe{
#'   \item{Age}{The year.}
#'   \item{SD}{The unscaled index of abundance.}
#'  }
#' 
"HSPs"

#' Gene tagging
#' 
#' Gene tagging (GT) observations.
#'
#' @format a \code{tibble} containing 2 fields including:
#' \describe{
#'   \item{Age}{The year.}
#'   \item{SD}{The unscaled index of abundance.}
#'  }
#' 
"GTs"

#' Aerial survey
#'
#' @format a \code{tibble} containing 3 fields including:
#' \describe{
#'   \item{Year}{The year.}
#'   \item{Unscaled_Index}{The unscaled index of abundance.}
#'   \item{CV}{The coefficient of variation.}
#'  }
#' 
"aerial_survey"

#' Aerial survey covariance
#'
#' @format a 20 x 20 \code{matrix}.
#' 
"aerial_cov"

#' CPUE
#'
#' @format a 20 x 20 \code{matrix}.
#' 
"cpue"

#' Troll
#'
#' @format a 20 x 20 \code{matrix}.
#' 
"troll"

#' AFs
#'
#' @format a 20 x 20 \code{matrix}.
#' 
"age_freq"

#' Proportion at length in the catch
#' 
#' Proportion at length in the catch observations for the LL1, LL2, LL3, and LL4 
#' fisheries
#'
#' @format a \code{tibble} containing 113 columns with length bins from 32 cm 
#' to 250 cm in 2 cm bin widths:
#' \describe{
#'   \item{Fishery}{The fishery.}
#'   \item{Year}{The year.}
#'   \item{N}{The effective sample size.}
#'   \item{32}{The proportion in the 32 cm bin.}
#'   \item{250}{The proportion in the 250 cm bin.}
#'  }
#'
"length_freq"

#' Tag release observations
#'
#' @format a \code{tibble} containing 5 fields including:
#' \describe{
#'   \item{Cohort}{Cohort.}
#'   \item{Group}{Tagger group.}
#'   \item{1}{The number of releases of age 1 fish.}
#'   \item{2}{The number of releases of age 2 fish.}
#'   \item{3}{The number of releases of age 3 fish.}
#'  }
#' 
"tag_releases"

#' Tag recapture observations
#'
#' @format a \code{tibble} containing 10 fields including:
#' \describe{
#'   \item{Cohort}{Cohort.}
#'   \item{Group}{Tagger group.}
#'   \item{RelAge}{The release age of the fish.}
#'   \item{1}{Recapture age 1.}
#'   \item{2}{Recapture age 2.}
#'   \item{3}{Recapture age 3.}
#'   \item{4}{Recapture age 4.}
#'   \item{5}{Recapture age 5.}
#'   \item{6}{Recapture age 6.}
#'   \item{7}{Recapture age 7.}
#'  }
#' 
"tag_recaptures"

#' Tag reporting rates
#'
#' @format a \code{tibble} containing 10 fields including:
#' \describe{
#'   \item{Year}{Cohort.}
#'   \item{Surf}{Tagger group.}
#'   \item{LL1}{The release age of the fish.}
#'   \item{1}{Recapture age 1.}
#'   \item{2}{Recapture age 2.}
#'   \item{3}{Recapture age 3.}
#'   \item{4}{Recapture age 4.}
#'   \item{5}{Recapture age 5.}
#'   \item{6}{Recapture age 6.}
#'   \item{7}{Recapture age 7.}
#'   \item{8}{Recapture age 8.}
#'  }
#' 
"tag_reporting"

#' Outputs from a single model run
#'
#' The \code{labrep1} data set is a \code{list} containing the \code{lab.rep} 
#' output file from a \code{base22sqrt} model run.
#' 
#' @format a list containing 100 fields including:
#' \describe{
#'   \item{scenario_number}{The scenario number for this labrep file.}
#'   \item{n}{The.}
#'   \item{years}{The first and last model years (i.e., 1931 and 2022).}
#'   \item{ObjF}{The objective function value.}
#'   \item{lnlike}{The log likelihood values.}
#'   \item{penal}{The prior and penalty values.}
#'   \item{Hhigh}{The.}
#'   \item{H2003}{The.}
#'   \item{H2004}{The.}
#'   \item{sigma.cpue}{The.}
#'   \item{tag.var.factor}{The.}
#'   \item{gtOD.factor}{The.}
#'   \item{tau.aerial}{The.}
#'   \item{tau.troll}{The.}
#'   \item{sigma.r}{Recruitment standard deviation.}
#'   \item{rec_AC}{The.}
#'   \item{AC_penalty}{The.}
#'   \item{q_AC}{The.}
#'   \item{res.stats}{The.}
#'   \item{B0}{The.}
#'   \item{R0}{The.}
#'   \item{alpha}{The.}
#'   \item{beta}{The.}
#'   \item{rho}{The.}
#'   \item{steep}{Steepness of the stock recruit relationship.}
#'   \item{depletion}{The.}
#'   \item{psi}{The.}
#'   \item{lnq}{The.}
#'   \item{qhsp}{The.}
#'   \item{qgt}{The.}
#'   \item{M}{The natural mortality at age.}
#'   \item{omega}{The.}
#'   \item{TOTbio}{The.}
#'   \item{B10_plus}{The.}
#'   \item{Sbio}{The spawning biomass.}
#'   \item{phi}{The.}
#'   \item{phsp}{The.}
#'   \item{Recruitment}{The.}
#'   \item{Rdev}{The.}
#'   \item{Reps}{The.}
#'   \item{ages}{The first and last model ages (i.e., 0 and 30).}
#'  }
#' 
"data_labrep1"

#' Parameter values from a single model run
#'
#' The \code{par1} data set is a \code{list} containing the \code{.par} 
#' output file from a \code{base22sqrt} model run.
#' 
#' @format a list containing 100 fields including:
#' \describe{
#'   \item{ln_B0}{The}
#'   \item{lnq}{the}
#'   \item{lnqhsp}{the}
#'   \item{deltalnq08}{the}
#'   \item{steep}{Steepness of the stock recruit relationship.}
#'   \item{sigma_r}{Recruitment standard deviation.}
#'   \item{sigma_cpue}{the}
#'   \item{m0}{Natural mortality at age 0.}
#'   \item{m4}{Natural mortality at age 4.}
#'   \item{m10}{Natural mortality at age 10.}
#'   \item{m30}{Natural mortality at age 30.}
#'   \item{cpue_omega}{the}
#'   \item{psi}{the}
#'   \item{Reps}{the}
#'   \item{ln_bdiff}{the}
#'   \item{ln_sel_aerial}{the}
#'   \item{tau_aerial}{the}
#'   \item{tau_troll}{the}
#'   \item{tag_H_factor}{the}
#'   \item{par_sels_init_i}{the initial selectivity by age for each fishery.}
#'   \item{par_sels_change_i}{the}
#'   \item{par_log_hstar_i}{the}
#'  }
#' 
"data_par1"

#' Tables read in from csv files
#'
#' The \code{csv1} data set is a \code{list} containing inputs read in from csv
#' files.
#' 
#' @format a list containing 100 fields including:
#' \describe{
#'   \item{length_mean}{The mean length at age by year and season.}
#'   \item{length_sd}{The standard devation (SD) of the mean length by age.}
#'   \item{catch}{The catch in tonnes by year for each fishery.}
#'   \item{catch_UA}{The unaccounted catch in tonnes by year for each fishery.}
#'   \item{scenarios_surface}{the}
#'   \item{scenarios_LL1}{the}
#'   \item{sel_change_sd}{the}
#'   \item{POPs}{the parent offspring pairs (POP) data.}
#'   \item{HSPs}{the half sibling paris (HSP) data.}
#'   \item{GTs}{the gene tagging (GT) data.}
#'   \item{aerial_surv}{the aerial survey data by year.}
#'   \item{aerial_cov}{the covariance matrix for the aerial survey data.}
#'   \item{cpue}{the catch per unit effort (CPUE) data by year.}
#'   \item{af}{the age frequency data for each fishery by year and age.}
#'   \item{lf}{the length frequency (LF) data for each fishery by year and 2 cm length bin.}
#'  }
#' 
"data_csv1"
