#' Simulated inputs for the assessment-tools tutorial
#'
#' A small two-fishery, two-season assessment with CPUE, length compositions,
#' and weight compositions. Observations are simulated from the selected
#' likelihood with known parameters; these are not data for a real stock.
#' @param family Composition likelihood for both streams.
#' @param seed Simulation seed, restored after generating the observations.
#' @return A list with `data`, `parameters`, `map`, and `truth`. Three fixed
#'   effects are estimated: log spawning output and two log catchabilities.
#' @family assessment workflow
#' @export
opal_example_inputs <- function(family = c("multinomial", "Dirichlet", "Dirichlet-multinomial"), seed = 712L) {
  family <- match.arg(family)
  switch <- match(family, c("multinomial", "Dirichlet", "Dirichlet-multinomial"))
  .opal_with_seed(seed, {
    years <- 2001:2012
    d <- list(n_year = 12L, n_season = 2L, n_age = 5L, n_fishery = 2L,
      n_len = 3L, n_wt = 3L, n_index = 2L, min_age = 1L, max_age = 5L,
      first_yr = 2001L, first_yr_catch = 2001L, last_yr = 2012L, years = years,
      A1 = 1, A2 = 5, len_lower = c(20,40,60), len_upper = c(40,60,80),
      len_mid = c(30,50,70), len_bin_start = 20, len_bin_width = 20,
      M = rep(.3,5), maturity = c(0,.2,.5,.8,1), fecundity = c(0,2,5,10,15),
      lw_a = 1e-5, lw_b = 3, catch_weight_units = "model weight units",
      catch_obs_ysf = array(rep(seq(30,70,length.out=12),4),c(12,2,2)),
      catch_units_f = c(1L,1L), removal_switch_f = c(0L,0L), sel_type_f = c(1L,1L),
      sel_fa_external = rbind(c(.2,.5,.8,1,1),c(.1,.4,.7,1,1)),
      cpue_switch = 1L, cpue_data = data.frame(ts=rep(1:24,2), fishery=rep(1:2,each=24),
        index=rep(1:2,each=24), units=1L, value=1, se=.1),
      wf_rebin_matrix = diag(3), bias_adj_y = rep(0,12))
    p <- list(log_B0=log(1e5),log_h=log(.75),log_sigma_r=log(.4),
      log_cpue_q=log(c(.8,1.2)),cpue_creep=c(0,0),log_cpue_tau=rep(log(.08),2),
      log_cpue_omega=c(0,0),log_lf_tau=rep(log(1),2),log_wf_tau=rep(log(1),2),
      log_L1=log(30),log_L2=log(65),log_k=log(.2),log_CV1=log(.15),log_CV2=log(.1),
      par_sel=matrix(0,2,6),rdev_y=rep(0,12))
    if (switch == 3L) p$log_lf_tau <- p$log_wf_tau <- rep(log(30),2)
    for (type in c("lf","wf")) {
      fields <- list(switch=switch,n_f=c(12L,12L),fishery_f=1:2,year_fi=list(1:12,1:12),
        n_fi=list(rep(120,12),rep(120,12)),n_int_fi=list(rep(120L,12),rep(120L,12)),
        minbin=c(1L,1L),maxbin=c(3L,3L),obs_flat=rep(40,72),obs_ints=rep(40L,72),obs_prop=rep(1/3,72))
      names(fields) <- paste0(type,"_",names(fields))
      d[names(fields)] <- fields
      d[[paste0("n_",type)]] <- 24L
    }
    map <- lapply(p,function(z) factor(rep(NA,length(z))))
    map$log_B0 <- map$log_cpue_q <- NULL
    d$priors <- list(log_B0=list(type="normal",par1=log(1e5),par2=.3,index=1L))
    object <- RTMB::MakeADFun(cmb(opal_model,d),p,map=map,silent=TRUE)
    report <- object$report(object$par)
    d$cpue_data$value <- exp(stats::rnorm(48,log(report$cpue_pred),report$cpue_sigma))
    for (type in c("lf","wf")) {
      predictions <- report[[paste0(type,"_pred")]]
      observation <- unlist(lapply(predictions,function(z) unlist(lapply(seq_len(nrow(z)),function(i) {
        prob <- z[i,]
        if (switch > 1L) {
          alpha <- prob * if (switch == 2L) 120 else 30
          prob <- stats::rgamma(3,shape=alpha); prob <- prob/sum(prob)
        }
        if (switch == 2L) prob else as.numeric(stats::rmultinom(1,120,prob))
      }))),use.names=FALSE)
      if (switch == 2L) {
        d[[paste0(type,"_obs_prop")]] <- observation
      } else {
        d[[paste0(type,"_obs_flat")]] <- observation
        d[[paste0(type,"_obs_ints")]] <- as.integer(observation)
        d[[paste0(type,"_obs_prop")]] <- observation/120
      }
    }
    list(data=d,parameters=p,map=map,truth=p)
  })
}
