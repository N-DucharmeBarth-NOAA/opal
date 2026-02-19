#' Set up the data input file
#' 
#' Set up the data input file to be passed to \code{MakeADFun}. This function 
#' runs data cross validation tests and appends several inputs to the data list 
#' including model dimensions and processed inputs:
#' 
#' * \code{n_year}: dervied from \code{first_yr} and \code{last_yr}
#' * \code{n_season}: set to 2
#' * \code{n_length}: not in use
#' * \code{n_age}: derived from \code{min_age} and \code{max_age}
#' * \code{n_fishery}: set to 6
#' * \code{age_a}: sequence of modeled ages derived from \code{min_age} and \code{max_age}
#' * \code{length_mu_ysa}: derived from the \code{length_mean} input
#' * \code{length_sd_a}: derived from the \code{length_sd} input
#' * \code{dl_yal}: derived from \code{length_mu_ysa} and \code{length_sd_a}
#' * \code{weight_fya}: derived from \code{length_mu_ysa} and \code{length_sd_a}
#' * \code{catch_obs_ysf}: derived from \code{catch}, \code{catch_UA}, \code{scenarios_LL1}, and \code{scenarios_surf}
#' * \code{sel_change_year_fy}: derived from \code{sel_change_sd_fy}
#' 
#' This function produces the data input file to be passed to \code{MakeADFun}.
#' 
#' @param data_in a \code{list} containing the data inputs.
#' @return a \code{list} ready to be passed to \code{MakeADFun}.
#' @importFrom testthat expect_identical expect_equal
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom RTMB pnorm
#' @export
#' 
get_data <- function(data_in) {

  # Dimensions ----
  
  data_in$first_yr <- 1931
  data_in$n_year <- length(data_in$first_yr:data_in$last_yr)
  data_in$n_season <- 2
  data_in$min_age <- 0
  data_in$max_age <- 30
  data_in$n_age <- length(data_in$min_age:data_in$max_age)
  data_in$n_length <- 1
  data_in$n_fishery <- 6
  data_in$age_a <- data_in$min_age:data_in$max_age
  
  fsh <- data.frame(ifishery = 1:6, 
                    fishery = c("LL1", "LL2", "LL3", "LL4", "Indonesia", "Australia"),
                    season = c(2, 2, 1, 1, 1, 1))
  
  # Length ----
  
  data_in$length_mu_ysa <- get_length_at_age(length_mean = sbt::length_mean)
  
  data_in$length_sd_a <- sbt::length_sd$SD
  names(data_in$length_sd_a) <- sbt::length_sd$Age
  
  expect_identical(dim(data_in$length_mu_ysa), 
                   as.integer(c(data_in$n_year, data_in$n_season, data_in$n_age)), 
                   info = "Dimension error in length_sd_a")
  expect_identical(length(data_in$length_sd_a), 
                   as.integer(data_in$n_age), 
                   info = "Dimension error in length_sd_a")
  
  # Weight ----
  
  data_in$weight_fya <- get_weight_at_age(length_mu_ysa = data_in$length_mu_ysa, length_sd_a = data_in$length_sd_a)
  
  # Catch ----
  
  data_in$first_yr_catch <- min(sbt::catch$Year)
  data_in$first_yr_catch_f <- c(1952, 1969, 1954, 1953, 1976, 1952)
  data_in$n_catch <- nrow(sbt::catch)
  data_in$catch_year <- sbt::catch$Year
  
  scenarios_LL1 <- data_in$scenarios_LL1 %>%
    select(Year, LL1_case = data_in$catch_LL1_case + 2) %>%
    mutate(fishery = "LL1")
  
  scenarios_surf <- data_in$scenarios_surf %>%
    select(Year, surf_case = data_in$catch_surf_case + 2) %>%
    mutate(fishery = "Australia")
  
  catch_UA <- sbt::catch_UA %>%
    pivot_longer(cols = -Year, names_to = "fishery", values_to = "UA")
  
  catch <- sbt::catch %>%
    pivot_longer(cols = -Year, names_to = "fishery") %>%
    full_join(scenarios_LL1, by = join_by("Year", "fishery")) %>%
    full_join(scenarios_surf, by = join_by("Year", "fishery")) %>%
    replace(is.na(.), 1) %>%
    full_join(catch_UA, by = join_by("Year", "fishery")) %>%
    replace(is.na(.), 0) %>%
    mutate(value = value * LL1_case, value = value * surf_case) %>%
    mutate(value = value + UA) %>%
    left_join(fsh, by = join_by("fishery")) %>%
    select(Year, season, ifishery, value) %>%
    pivot_wider(names_from = ifishery, values_from = value, values_fill = 0)

  data_in$catch_obs_ysf <- array(data = 0, 
                                 dim = c(data_in$n_catch, data_in$n_season, data_in$n_fishery),
                                 dimnames = list(Year = data_in$catch_year, Season = 1:2, Fishery = fsh$fishery))
  data_in$catch_obs_ysf[,1,3:6] <- as.matrix(catch %>% filter(season == 1) %>% select(`3`, `4`, `5`, `6`))
  data_in$catch_obs_ysf[,2,1:2] <- as.matrix(catch %>% filter(season == 2) %>% select(`1`, `2`))

  expect_identical(dim(data_in$catch_obs_ysf), 
                   as.integer(c(data_in$n_catch, data_in$n_season, data_in$n_fishery)), 
                   info = "Dimension error in catch_obs_ysf")
  
  # Selectivity ----
  
  expect_equal(length(data_in$sel_min_age_f), 7, info = "Dimension error in sel_min_age_f.")
  expect_equal(length(data_in$sel_max_age_f), 7, info = "Dimension error in sel_max_age_f.")
  expect_equal(length(data_in$sel_end_f), 7, info = "Dimension error in sel_end_f.")
  
  sel_change_year_fy <- array(0, dim = c(7, data_in$n_year))
  dimnames(sel_change_year_fy) <- list(fishery = c("LL1", "LL2", "LL3", "LL4", "Indonesia", "Australia", "CPUE"), 
                                       year = data_in$first_yr:data_in$last_yr)
  sel_change_year_fy[1,] <- ifelse(colnames(sel_change_year_fy) %in% data_in$sel_LL1_yrs, 1, 0)
  sel_change_year_fy[2,] <- ifelse(colnames(sel_change_year_fy) %in% data_in$sel_LL2_yrs, 1, 0)
  sel_change_year_fy[3,] <- ifelse(colnames(sel_change_year_fy) %in% data_in$sel_LL3_yrs, 1, 0)
  sel_change_year_fy[4,] <- ifelse(colnames(sel_change_year_fy) %in% data_in$sel_LL4_yrs, 1, 0)
  sel_change_year_fy[5,] <- ifelse(colnames(sel_change_year_fy) %in% data_in$sel_Ind_yrs, 1, 0)
  sel_change_year_fy[6,] <- ifelse(colnames(sel_change_year_fy) %in% data_in$sel_Aus_yrs, 1, 0)
  sel_change_year_fy[7,] <- ifelse(colnames(sel_change_year_fy) %in% data_in$sel_CPUE_yrs, 1, 0)
  data_in$sel_change_year_fy <- sel_change_year_fy
  
  # POPs ----

  data_in$pop_obs <- sbt::POPs %>%
    filter(Comps > 0) %>%
    mutate(Cohort = Cohort - data_in$first_yr + 1) %>%
    mutate(CaptureYear = CaptureYear - data_in$first_yr + 1) %>%
    mutate(CaptureCov = ifelse(CaptureSwitch == 0, CaptureCov - data_in$min_age, CaptureCov)) %>%
    # mutate(CaptureCov = ifelse(CaptureSwitch == 1, CaptureCov, CaptureCov)) %>%
    select(Cohort, CaptureYear, CaptureCov, CaptureSwitch, NPOPS, Comps) %>%
    as.matrix()
  
  # paly ----

  paly <- sbt::paly
  xbins <- dim(paly)[1]
  xages <- as.character(dimnames(paly)[[2]])
  xyrs <- as.character(dimnames(paly)[[3]])
  xpaly <- array(0, dim = c(xbins, data_in$n_age, data_in$n_year))
  dimnames(xpaly)[[2]] <- as.character(data_in$age_a)
  dimnames(xpaly)[[3]] <- as.character(data_in$first_yr:data_in$last_yr) 
  xpaly[, 1:min(as.numeric(xages)),] <- 0 
  xpaly[, xages, xyrs] <- paly
  data_in$paly <- xpaly

  # HSPs ----
  
  data_in$hsp_obs <- sbt::HSPs %>%
    mutate(cohort1 = cohort1 - data_in$first_yr + 1) %>% 
    mutate(cohort2 = cohort2 - data_in$first_yr + 1) %>% 
    rowwise() %>%
    mutate(cmin = min(cohort1, cohort2), cmax = max(cohort1, cohort2)) %>%
    mutate(cdiff = cmax - cmin) %>%
    select(cmin, cmax, cdiff, nC, nK) %>%
    as.matrix()
  
  # Gene tagging (GT) ----
  
  data_in$gt_obs <- sbt::GTs %>%
    mutate(RelYear = RelYear - data_in$first_yr + 1) %>% 
    mutate(RecYear = RecYear - data_in$first_yr + 1) %>% 
    mutate(RelAge = RelAge + 1) %>% # change to index
    as.matrix()
  data_in$gt_nscan <- data_in$gt_obs[, 5]
  data_in$gt_nrec <- data_in$gt_obs[, 6]
  
  # Aerial surveys ----
  
  data_in$aerial_years <- sbt::aerial_survey$Year - data_in$first_yr + 1
  data_in$aerial_obs <- sbt::aerial_survey$Unscaled_Index
  data_in$aerial_cv <- sbt::aerial_survey$CV
  data_in$aerial_cov <- sbt::aerial_cov
  
  expect_equal(dim(data_in$aerial_cov), rep(length(data_in$aerial_years), 2), info = "Dimension error in aerial_cov.")
  
  # Troll surveys ----
  
  data_in$troll_years <- sbt::troll$Year - data_in$first_yr + 1
  data_in$troll_obs <- sbt::troll$Median
  data_in$troll_sd <- sbt::troll$SD
  
  # CPUE ----
  
  data_in$cpue_years <- sbt::cpue$Year - data_in$first_yr + 1
  data_in$cpue_obs <- sbt::cpue$CPUE / mean(sbt::cpue$CPUE)
  data_in$cpue_sd <- numeric(length(sbt::cpue$CPUE))

  # Age-frequency ----
  
  if (is.null(data_in$af_data)) {
    af_data <- sbt::age_freq
  } else {
    af_data <- data_in$af_data
  }
  
  data_in$af_min_age <- af_data$MinAge
  data_in$af_max_age <- af_data$MaxAge
  
  Pt <- data_in$scenarios_surf %>%
    select(Year, Pt = P_t_20) %>%
    mutate(Fishery = 6)
  
  af1 <- af_data %>% 
    pivot_longer(cols = !c(1:5), names_to = "Age") %>%
    mutate(Age = as.numeric(Age)) %>%
    mutate(Age = ifelse(Age > data_in$max_age, data_in$max_age, Age)) %>% # Ages go up to 40 so need to aggregate into max_age of 30
    mutate(Age = ifelse(Age < MinAge, MinAge, Age)) %>%
    mutate(Age = ifelse(Age > MaxAge, MaxAge, Age)) %>%
    mutate(Age = factor(Age, levels = 0:30)) %>%
    group_by(Fishery, Year, N, Age) %>%
    summarise(value = sum(value, na.rm = TRUE), .groups = "drop") %>%
    ungroup() %>%
    # complete(Age, nesting(year, month, area), fill = list(count = 0)) %>%
    pivot_wider(names_from = Age, names_expand = TRUE, values_fill = 0)# %>% 
    # left_join(Pt, by = join_by("Fishery", "Year"))
  
  for (i in 1:nrow(af1)) {
    if (af1$Fishery[i] == 6 & af1$Year[i] >= 1992 & data_in$catch_surf_case > 0) {
      obs2 <- as.numeric(af1[i, "2"])
      obs3 <- as.numeric(af1[i, "3"])
      pp <- Pt$Pt[Pt$Year == af1$Year[i]]
      af1[i, "2"] <- af1[i, "2"] * (1 - pp) # obs_age_freq_ija(in,irec,2)*= (1.-Pt(iy));
      af1[i, "3"] <- (1 - pp) * (obs3 + pp * obs2) # obs_age_freq_ija(in,irec,3) = (1.-Pt(iy))*(obs3 + Pt(iy)* obs2);
      af1[i, "4"] <- af1[i, "4"] + pp * (obs3 + pp * obs2) # obs_age_freq_ija(in,irec,4)+= Pt(iy)*(obs3 + Pt(iy)* obs2);
    }
  }
  
  # vector Pt(1992,last_yr); // to adjust age comp of surface fishery
  # if(iff==6 && surf_case > 0){
  #   // adjusts catch in weight
  #   catch_fy(iff)(1992,last_yr) = elem_prod(catch_fy(iff)(1992,last_yr), column(surf_scen,2*surf_case-1+2));
  #   // to adjust age composition of surface fishery
  #   Pt = column(surf_scen,2*surf_case+2);
  # }
  # if (iff == 6 && iy >= 1992 && surf_case > 0) {
  #   obs2 = obs_age_freq_ija(in,irec,2);
  #   obs3 = obs_age_freq_ija(in,irec,3);
  #   obs_age_freq_ija(in,irec,2)*= (1.-Pt(iy));
  #   obs_age_freq_ija(in,irec,3) = (1.-Pt(iy))*(obs3 + Pt(iy)* obs2);  
  #   obs_age_freq_ija(in,irec,4)+= Pt(iy)*(obs3 + Pt(iy)* obs2);  
  # }
  
  data_in$n_af <- nrow(af1)
  data_in$af_year <- af1$Year - data_in$first_yr + 1
  data_in$af_fishery <- af1$Fishery
  data_in$af_obs <- af1 %>% select(-Fishery, -Year, -N) %>% as.matrix()
  data_in$af_obs[is.na(data_in$af_obs)] <- 0
  data_in$af_n <- af1$N
  
  # Age-length key ----
  
  data_in$dl_yal <- get_dl(length_mu_ysa = data_in$length_mu_ysa, length_sd_a = data_in$length_sd_a)
  
  # MOVE TO ITS OWN FUNCTION LIKE get_dl - COMPARE WITH output from get_dl
  min_len <- 86
  bin_width <- 4
  nbins <- 25
  
  alk_ysal <- array(NA, dim = c(data_in$n_year, 2, data_in$n_age, nbins))
  
  for (y in 1:data_in$n_year) {
    for (s in 1:2) {
      for (a in 1:data_in$n_age) {
        mu_len <- data_in$length_mu_ysa[y, s, a]
        cumhld <- 0.0
        for (l in 1:(nbins - 1)) {
          bin_max_len <- min_len + bin_width * (l - 1)
          cum <- pnorm((bin_max_len - mu_len) / data_in$length_sd_a[a])
          alk_ysal[y, s, a, l] <- cum - cumhld
          cumhld <- cum
        }
        alk_ysal[y, s, a, nbins] <- 1.0 - cumhld
      }
    }
  }
  
  data_in$alk_ysal <- alk_ysal
  
  # y <- 40; s <- 1; a <- 6
  # ll <- seq(from = min_len, by = bin_width, length.out = nbins)
  # plot(ll, alk_ysal[y, s, a, ])
  # lines(data_in$dl_yal[y, a, ])
  # // FUNCTION get_age_length_key
  # // int is, iy, ia, ib;
  # // double cum,cumhld;
  # // double mu_len,bin_max_len;
  # // for (is=1; is<=2; is++) {
  #   //   for (iy=first_yr; iy<=last_yr; iy++) {
  #     //     mean_len_age(is,iy)=input_len_age(is,iy)(0,last_age);
  #     //     for (ia=0; ia<=last_age; ia++) {
  #       //       mu_len=mean_len_age(is,iy,ia);
  #       //       cumhld=0.;
  #       //       for (ib=1; ib<nbins; ib++) {
  #         //         bin_max_len                   = min_len+bin_width*(ib-1);
  #         //         cum                           = cumd_norm( (bin_max_len-mu_len)/std_len(ia));
  #         //         lenage_dist_syal(is,iy,ia,ib) = cum-cumhld;
  #         //         cumhld                        = cum;
  #         //       }
  #       //       lenage_dist_syal(is,iy,ia,nbins)= 1.-cumhld;
  #       //     }
  #     //   }
  #   // }
  
  # Length-frequency (LF) ----

  get_lf_obs <- function(lf_data, nbins = 25, min_len = 86, bin_width = 4) {
    obs_len_freq_il <- matrix(0, nrow = nrow(lf_data), ncol = nbins)
    ncol_lf <- ncol(lf_data) - 3
    for (irec in 1:nrow(lf_data)) {
      kbin <- 1
      mod_bin_wid <- min_len + bin_width * (kbin - 1)
      for (i in seq_len(ncol_lf)) {
        obs_bin_wid <- 32 + 2 * (i - 1);
        if (obs_bin_wid > mod_bin_wid && kbin < nbins) {
          kbin <- kbin + 1
          mod_bin_wid <- min_len + bin_width * (kbin - 1)
        }
        obs_len_freq_il[irec, kbin] = obs_len_freq_il[irec, kbin] + as.numeric(lf_data[irec, -c(1:3)][i])
      }
      obs_len_freq_il[irec,] = obs_len_freq_il[irec,] / sum(obs_len_freq_il[irec,])
    }
    return(obs_len_freq_il)
  }
  
  # if (is.null(data_in$lf_data)) {
    lf_data <- sbt::length_freq %>% filter(Fishery != 7)
  # } else {
    # lf_data <- data_in$lf_data
  # }
  obs_len_freq_il <- get_lf_obs(lf_data, nbins = 25, min_len, bin_width)
  ll <- seq(from = min_len, by = bin_width, length.out = nbins + 1) - 1
  
  lf1 <- lf_data %>%
    pivot_longer(cols = !c(1:3), names_to = "Bin") %>%
    mutate(Bin = as.numeric(Bin)) %>%
    mutate(dBin = cut(Bin, breaks = ll, include.lowest = TRUE, right = FALSE)) %>%
    filter(!is.na(dBin)) %>%
    group_by(Fishery, Year, N, dBin) %>%
    summarise(value = sum(value, na.rm = TRUE)) %>%
    ungroup() %>%
    pivot_wider(names_from = dBin) %>%
    left_join(fsh, by = join_by(Fishery == ifishery)) %>%
    relocate(Fishery, Year, season, N) %>%
    select(-fishery)

  data_in$n_lf <- nrow(lf_data)
  data_in$lf_year <- lf_data$Year - data_in$first_yr + 1
  data_in$lf_fishery <- lf_data$Fishery
  data_in$lf_season <- lf1$season
  data_in$lf_obs <- obs_len_freq_il
  data_in$lf_n <- lf_data$N
  
  # Cohort slice the LFs ----
  
  afs1 <- get_sliced_afs(data = data_in, lf_data = lf_data)
  data_in$lf_slices <- afs1$lf_slices
  data_in$af_sliced <- afs1$af_sliced
  
  sliced_ysfa <- array(data = 0, 
                       dim = c(data_in$n_year, data_in$n_season, data_in$n_fishery, data_in$n_age),
                       dimnames = list(Year = data_in$first_yr:data_in$last_yr, Season = 1:2, Fishery = fsh$fishery, Age = data_in$age_a))
  
  sliced_ysfa[data_in$lf_year[data_in$lf_fishery == 1], 2, 1,] <- data_in$af_sliced[data_in$lf_fishery == 1,]
  sliced_ysfa[data_in$lf_year[data_in$lf_fishery == 2], 2, 2,] <- data_in$af_sliced[data_in$lf_fishery == 2,]
  sliced_ysfa[data_in$lf_year[data_in$lf_fishery == 3], 1, 3,] <- data_in$af_sliced[data_in$lf_fishery == 3,]
  sliced_ysfa[data_in$lf_year[data_in$lf_fishery == 4], 1, 4,] <- data_in$af_sliced[data_in$lf_fishery == 4,]
  sliced_ysfa[data_in$af_year[data_in$af_fishery == 5], 1, 5,] <- data_in$af_obs[data_in$af_fishery == 5,]
  sliced_ysfa[data_in$af_year[data_in$af_fishery == 6], 1, 6,] <- data_in$af_obs[data_in$af_fishery == 6,]
  data_in$af_sliced_ysfa <- sliced_ysfa
  
  # CPUE LFs ----
  
  cpue_lfs <- sbt::length_freq %>% filter(Fishery == 7)
  obs_len_freq_il <- get_lf_obs(cpue_lfs, nbins = 25) 
  ll <- seq(from = min_len, by = bin_width, length.out = nbins + 1) - 1
  
  lf1 <- cpue_lfs %>%
    pivot_longer(cols = !c(1:3), names_to = "Bin") %>%
    mutate(Bin = as.numeric(Bin)) %>%
    mutate(dBin = cut(Bin, breaks = ll, include.lowest = TRUE, right = FALSE)) %>%
    filter(!is.na(dBin)) %>%
    group_by(Fishery, Year, N, dBin) %>%
    summarise(value = sum(value, na.rm = TRUE)) %>%
    ungroup() %>%
    pivot_wider(names_from = dBin) %>%
    # left_join(fsh2, by = join_by(Fishery == ifishery)) %>%
    mutate(ifishery = 7, fishery = "CPUE", season = 2) %>%
    relocate(Fishery, Year, season, N) %>%
    select(-fishery)
  
  data_in$cpue_lfs <- obs_len_freq_il
  data_in$cpue_n <- cpue_lfs$N
  
  # Tagging ----
  
  data_in$tag_shed_immediate <- c(0.9737, 0.9608, 1, 1, 0.9342, 0.9666)
  data_in$tag_shed_continuous <- c(0.0391, 0.0492, 0.0672, 0.0925, 0.0885, 0.1601)
  
  data_in$tag_rep_rates_ya <- sbt::tag_reporting %>%
    filter(LL1 == data_in$catch_LL1_case, Surf == data_in$catch_surf_case) %>%
    select(-c(1:3)) %>%
    as.matrix()
  
  data_in$scenarios_LL1 %>%
    select(Year, LL1_case = data_in$catch_LL1_case + 2) %>%
    mutate(fishery = "LL1")
  
  scenarios_surf <- data_in$scenarios_surf
  
  df <- sbt::tag_releases %>% 
    pivot_longer(3:5, names_to = "Age") %>%
    arrange(Cohort, Group, Age)
  
  tag_rel_age <- df %>% 
    filter(value > 0) %>% 
    group_by(Cohort) %>% 
    summarise(min_age = min(Age), max_age = max(Age))
  data_in$tag_rel_min_age <- as.numeric(tag_rel_age$min_age)
  data_in$tag_rel_max_age <- as.numeric(tag_rel_age$max_age)

  a1 <- length(unique(df$Cohort))
  a2 <- length(unique(df$Group))
  a3 <- length(unique(df$Age))
  data_in$tag_release_cta <- array(NA, dim = c(a1, a2, a3))
  
  for (i1 in 1:a1) {
    for (i2 in 1:a2) {
      for (i3 in 1:a3) {
        data_in$tag_release_cta[i1, i2, i3] <- df %>% 
          filter(Cohort == unique(df$Cohort)[i1], Group == unique(df$Group)[i2], Age == unique(df$Age)[i3]) %>% 
          select(value) %>% 
          as.numeric()
      }
    }
  }
  
  df <- sbt::tag_recaptures %>% 
    pivot_longer(4:10, names_to = "RecAge") %>%
    arrange(Cohort, Group, RelAge, RecAge)
  
  a1 <- length(unique(df$Cohort))
  a2 <- length(unique(df$Group))
  a3 <- length(unique(df$RelAge))
  a4 <- length(unique(df$RecAge))
  data_in$tag_recap_ctaa <- array(NA, dim = c(a1, a2, a3, a4))
  
  for (i1 in 1:a1) {
    for (i2 in 1:a2) {
      for (i3 in 1:a3) {
        for (i4 in 1:a4) {
          data_in$tag_recap_ctaa[i1, i2, i3, i4] <- df %>% 
            filter(Cohort == unique(df$Cohort)[i1], Group == unique(df$Group)[i2],
                   RelAge == unique(df$RelAge)[i3], RecAge == unique(df$RecAge)[i4]) %>% 
            select(value) %>% 
            as.numeric()
        }
      }
    }
  }
  
  tag_rec_age <- df %>% 
    filter(value > 0) %>% 
    group_by(Cohort) %>% 
    summarise(max_age = max(RecAge))
  data_in$tag_recap_max_age <- as.numeric(tag_rec_age$max_age)
  
  data_in$min_K <- 1989 - data_in$first_yr # first tagged cohorts (1989)
  max_K <- 1994 - data_in$first_yr # last tagged cohorts (1994)
  data_in$n_K <- max_K - data_in$min_K + 1 # number of cohorts (6)
  data_in$n_T <- 6 # number of tagger groups (6)
  data_in$n_I <- 3 # number of rel ages being included (3)
  data_in$n_J <- 7 # number of recap ages being included  (7)

  # Reducing output ----
  
  data_in$catch_UR_on <- NULL
  data_in$scenarios_LL1 <- data_in$scenarios_surf <- NULL
  data_in$catch_LL1_case <- data_in$catch_surf_case <- NULL
  data_in$nf <- NULL
  data_in$sel_LL1_yrs <- data_in$sel_LL2_yrs <- data_in$sel_LL3_yrs <- data_in$sel_LL4_yrs <- NULL
  data_in$sel_Ind_yrs <- data_in$sel_Aus_yrs <- data_in$sel_CPUE_yrs <- NULL
  data_in$n_lf <- NULL
  
  return(data_in)
}
