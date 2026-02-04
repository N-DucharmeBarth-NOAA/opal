

# Nicholas Ducharme-Barth
# 2026/01/06
# R code to make a baseline stock synthesis model version of the mfcl/v11 model
# Make sure rates and time-steps are adjusted for a seasonal model

# Copyright (c) 2026 Nicholas Ducharme-Barth
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

#_____________________________________________________________________________________________________________________________
# load packages
	library(data.table)
	library(magrittr)
	library(FLR4MFCL)
    library(frqit)
    library(r4ss)

#_____________________________________________________________________________________________________________________________
# define paths
	proj_dir = this.path::this.proj()
	dir_model = file.path(proj_dir,"model-files")
    dir_base_mfcl = file.path(dir_model,"mfcl","v11")
    dir_base_stock_synthesis = file.path(dir_model,"ss3","00-swpo-mls-base-file")
    dir_helper_fns_ss3 = file.path(proj_dir,"code","ss3","helper-fns")
    dir_helper_fns_mfcl = file.path(proj_dir,"code","mfcl","helper-fns")

#________________________________________________________________________________________________________________________________________________________________________________________________________
# source helper functions
    sapply(file.path(dir_helper_fns_ss3,(list.files(dir_helper_fns_ss3))),source)
    sapply(file.path(dir_helper_fns_mfcl,(list.files(dir_helper_fns_mfcl))),source)

#_____________________________________________________________________________________________________________________________
# read in baseline mfcl files
    base_frq = parse_frq(file.path(dir_base_mfcl,"bet.frq"))
    base_ini = read.MFCLIni(file.path(dir_base_mfcl,"bet.ini"), nseasons=4)
    base_par = read.MFCLPar(file.path(dir_base_mfcl,"10.par"), first.yr=1952)
    base_rep = read.MFCLRep(file.path(dir_base_mfcl,"plot-10.par.rep"))

#_____________________________________________________________________________________________________________________________
# run new version of stock synthesis
    # file.copy(from=paste0(proj_dir,"/executables/stock-synthesis/3.30.24.1/ss3_win.exe"),to=dir_base_stock_synthesis)
    # run(dir=dir_base_stock_synthesis,exe="ss3_win.exe")

#_____________________________________________________________________________________________________________________________
# read in baseline stock synthesis files
    tmp_starter = SS_readstarter(file=file.path(dir_base_stock_synthesis,"starter.ss_new"),verbose=FALSE)
    tmp_ctl = SS_readctl(file=file.path(dir_base_stock_synthesis,"control.ss_new"),datlist = file.path(dir_base_stock_synthesis,"data_echo.ss_new"))
    tmp_data = SS_readdat(file=file.path(dir_base_stock_synthesis,"data_echo.ss_new"))
    tmp_forecast = SS_readforecast(file=file.path(dir_base_stock_synthesis,"forecast.ss_new"),verbose=FALSE)

#_____________________________________________________________________________________________________________________________
# create new directory for stock synthesis bet files
    dir_bet_stock_synthesis_base = file.path(dir_model,"ss3","xx-test-matrix-error")
    dir.create(dir_bet_stock_synthesis_base,recursive=TRUE)


#_____________________________________________________________________________________________________________________________
# update starter
    tmp_starter$init_values_src = 0
    tmp_starter$last_estimation_phase = 10

    tmp_starter$min_age_summary_bio = 12
    tmp_starter$depl_basis = 5
    tmp_starter$SPR_basis = 4
    tmp_starter$F_std_units = 4
    tmp_starter$F_age_range = c(12,40)
    SS_writestarter(tmp_starter,dir=dir_bet_stock_synthesis_base,overwrite=TRUE)

#_____________________________________________________________________________________________________________________________
# update forecast
    tmp_forecast$SPRtarget = 0.3
    tmp_forecast$Btarget = 0.3
    SS_writeforecast(tmp_forecast,dir=dir_bet_stock_synthesis_base,overwrite=TRUE)

#_____________________________________________________________________________________________________________________________
# update data
    tmp_data$Comments = "#C No comments"
    tmp_data$styr = 1
    tmp_data$endyr = (diff(range(cateffpen(base_frq)$year,na.rm=TRUE))+1)*4
    tmp_data$nseas = 1
    tmp_data$months_per_seas = 12
    tmp_data$spawn_month = 1
    tmp_data$spawn_seas = 1 # deprecated
    tmp_data$Nsexes = 1
    tmp_data$Nages = unname(dimensions(base_ini)["agecls"])
    tmp_data$Nareas = 1
    tmp_data$Nfleets = n_fisheries(base_frq)
    tmp_data$Nfleet = tmp_data$Nfleets - 1
    tmp_data$Nsurveys = 1

# define fleet info
    tmp_fleetname = c("F01_LL.NORTH",      # Northern longline (Fisheries 1 & 2)
                        "F02_LL.US",         # US longline (Fishery 3)
                        "F03_LL.OFFSHORE",   # Offshore longline (Fisheries 5 & 6)
                        "F04_LL.EQUAT",      # Equatorial longlines (Fisheries 4, 8 & 9)
                        "F05_LL.WEST",       # Western longline (Fishery 7)
                        "F06_LL.SOUTH",      # Southern longline (Fisheries 11, 12 & 29)
                        "F07_LL.AU",         # Australian longline (Fisheries 10 & 27)
                        "F08_PS.ASS",        # Associated purse seines (Fisheries 13, 15, 25, 30, & 24)
                        "F09_PS.UNA",        # Unassociated purse seines (14, 16, 26, & 31)
                        "F10_DOM.MISC",      # Domestic miscellaneous fisheries (Fisheries 17, 23, & 22)
                        "F11_DOM.HL",        # Domestic handline (Fishery 18)
                        "F12_PS.JP.NORTH",   # Northern Japanese purse-seine (Fishery 19)
                        "F13_PL.JP.NORTH",   # Northern Japanese pole-and-line (Fishery 20)
                        "F14_PL.EQUAT",      # Equatorial pole-and-line (Fisheries 21, 22, & 28)
                        "S01_INDEX")         # Index fishery (Fisheries 33-41)
    tmp_fleetinfo = data.frame(type=c(rep(1,tmp_data$Nfleets-1),3),
                               surveytiming=c(rep(-1,tmp_data$Nfleets-1),1),
                               area=rep(1,tmp_data$Nfleets),
                               units=c(2,2,2,2,2,2,2,1,1,1,1,1,1,1,2),
                               need_catch_mult=rep(0,tmp_data$Nfleets),
                               fleetname=tmp_fleetname)
    tmp_data$fleetinfo = tmp_fleetinfo
    tmp_data$fleetnames = tmp_fleetinfo$fleetname
    tmp_data$surveytiming = tmp_fleetinfo$surveytiming
    tmp_data$units_of_catch = tmp_fleetinfo$units
    tmp_data$areas = tmp_fleetinfo$area

    tmp_data$fleetinfo1 = rbind(tmp_data$surveytiming,tmp_data$areas,c(rep(1,tmp_data$Nfleets-1),3))
    rownames(tmp_data$fleetinfo1) = c("surveytiming","areas","type")
    colnames(tmp_data$fleetinfo1) = tmp_data$fleetnames

    tmp_data$fleetinfo2 = rbind(tmp_data$units_of_catch,rep(0,tmp_data$Nfleets))
    rownames(tmp_data$fleetinfo2) = c("units","need_catch_mult")
    colnames(tmp_data$fleetinfo2) = tmp_data$fleetnames

    # define catch
    # no equilibrium catch for any fishery
    # remember catch in number needs to be divided by 1000 (stock synthesis records catch in 1000s of fish)
    ts_dt = data.table(expand.grid(year=min(range(cateffpen(base_frq)$year,na.rm=TRUE)):max(range(cateffpen(base_frq)$year,na.rm=TRUE)),month=c(2,5,8,11))) %>%
            .[order(year,month)] %>%
            .[,ts:=1:.N]
    
    
    cateffpen_dt = as.data.table(cateffpen(base_frq)) %>%
                   .[,units:=tmp_fleetinfo$units[fishery]] %>%
                   .[units==2,catch:=catch/1000] %>%
                   merge(.,ts_dt,by=c("year","month"))

    tmp_catch.list = as.list(rep(NA,tmp_data$Nfleets-1))
    for(i in seq_along(tmp_catch.list)){
        tmp_catch.list[[i]] = cateffpen_dt[fishery==i,.(ts,fishery,catch)] %>%
                         setnames(.,c("ts","fishery"),c("year","fleet")) %>%
                         .[,seas:=1] %>%
                         .[,catch_se:=0.01] %>%
                         .[,.(year,seas,fleet,catch,catch_se)] %>%
                         .[catch>0]
        # add equilibrium catch (0)
            tmp_catch.list[[i]] = rbind(data.table(year=-999,seas=1,fleet=i,catch=0,catch_se=0.01),tmp_catch.list[[i]])
    }
    tmp_catch = rbindlist(tmp_catch.list) %>%
                as.data.frame(.)
    tmp_data$catch = tmp_catch

    # CPUEinfo
    tmp_cpueinfo = data.frame(fleet=1:tmp_data$Nfleets,
                              units= -tmp_fleetinfo$units+2,
                              errtype=rep(0,tmp_data$Nfleets),
                              SD_report=rep(0,tmp_data$Nfleets))
    rownames(tmp_cpueinfo) = tmp_fleetname
    tmp_data$CPUEinfo = tmp_cpueinfo

    # define CPUE
    # penalty column appears to be cv
    # use fishery 15
    tmp_cpue = cateffpen_dt[fishery==15,.(cpue=catch/effort,penalty=penalty),by=.(fishery,year,month)] %>%
               .[,obs:=cpue/mean(cpue)] %>%
               .[,cv:=1/sqrt(2*penalty)] %>%
               .[,se_log:=sqrt(log(1+cv^2))] %>%
                merge(.,ts_dt,by=c("year","month")) %>%
                .[,year:=NULL] %>%
                .[,month:=1] %>%
                setnames(.,c("ts","fishery"),c("year","index")) %>%
               .[,.(year,month,index,obs,se_log)] %>%
               as.data.frame(.)
    tmp_data$CPUE = tmp_cpue

    # define population length structure
    mfcl_bin_lower = seq(from=lf_range(base_frq)[3],by=lf_range(base_frq)[4],length.out=lf_range(base_frq)[2])
    mfcl_bin_mid = seq(from=lf_range(base_frq)[3],by=lf_range(base_frq)[4],length.out=lf_range(base_frq)[2]) + 0.5*lf_range(base_frq)[4]
    mfcl_bin_upper = seq(from=lf_range(base_frq)[3],by=lf_range(base_frq)[4],length.out=lf_range(base_frq)[2]) + lf_range(base_frq)[4]

    tmp_data$minimum_size = 10
    tmp_data$maximum_size = 200
    tmp_data$binwidth = 2
    tmp_data$lbin_vector_pop = seq(from=tmp_data$minimum_size,to=tmp_data$maximum_size,by=tmp_data$binwidth)
    tmp_data$N_lbinspop = length(tmp_data$lbin_vector_pop)

    # define length info
    tmp_data$use_lencomp = 1
    tmp_leninfo = data.table(mintailcomp=rep(-1e-04,nrow(tmp_data$fleetinfo))) %>%
                  .[,addtocomp:=1e-04] %>%
                  .[,combine_M_F:=0] %>%
                  .[,CompressBins:=0] %>%
                  .[,CompError:=0] %>%
                  .[,ParmSelect:=0] %>%
                  .[,minsamplesize:=0.001] %>%
                  as.data.frame(.)
    rownames(tmp_leninfo) = tmp_fleetname
    tmp_data$len_info = tmp_leninfo
    tmp_data$N_lbins = length(mfcl_bin_lower)
    tmp_data$lbin_vector = mfcl_bin_lower

    # define length composition data
    tmp_lencomp = as.data.table(lnfrq(base_frq)) %>%
                  na.omit(.) %>%
                  as.data.frame(.)
    tmp_lencomp_a = tmp_lencomp[,c("year","month","fishery")]
    colnames(tmp_lencomp_a) = c("year","month","fleet")

    tmp_lencomp_a = as.data.table(tmp_lencomp_a) %>%
                    merge(.,ts_dt,by=c("year","month")) %>%
                    .[,month:=1] %>%
                    .[,.(ts,month,fleet)] %>%
                    setnames(.,c("ts"),c("year")) %>%
                    as.data.frame(.)

    tmp_lencomp_b = tmp_lencomp[,-c(1:4)]
    tmp_lencomp_c = tmp_lencomp_b
    tmp_lencomp_c[tmp_lencomp_c>0] = 0
    colnames(tmp_lencomp_c) = paste0("m",colnames(tmp_lencomp_b))
    colnames(tmp_lencomp_b) = paste0("f",colnames(tmp_lencomp_b))
    
    tmp_lencomp_a$sex = 0
    tmp_lencomp_a$part = 0
    tmp_lencomp_a$Nsamp = rowSums(tmp_lencomp_b)
    # tmp_lencomp_a$Nsamp = sapply(rowSums(tmp_lencomp_b),function(x)min(c(x,1000))) # apply MFCL like filter to maximum sample size

    # exclude males

    tmp_lencomp = as.data.table(cbind(tmp_lencomp_a,tmp_lencomp_b))

    tmp_data$lencomp = as.data.frame(tmp_lencomp)

    # update generalized size comp
    # tmp_data$N_sizefreq_methods_rd = 0
    # tmp_data$N_sizefreq_methods = 0
    tmp_data$N_sizefreq_methods_rd = 1
    tmp_data$N_sizefreq_methods = 1
    
    tmp_data$units_per_method = 2
    tmp_data$scale_per_method = 1
    tmp_data$mincomp_per_method = 0.001

    
    tmp_wtcomp_init = as.data.table(wtfrq(base_frq)) %>%
                      melt(.,id.vars=c("year","month","week","fishery")) %>%
                      setnames(.,"variable","bin") %>%
                      .[,bin:=as.numeric(as.character(bin))]
    
    # need to re-aggregate comp data beginning at old bin 44kg to make 2 kg
    # need to re-aggregate comp data beginning at old bin 180 to make 4 kg
    old_wt_bins = seq(from=44,to=uniqueN(tmp_wtcomp_init$bin),by=1)
    new_wt_bins = rep(NA,length(old_wt_bins))
    for(i in seq_along(old_wt_bins)){
        # if(i %% 2 != 0){
        #     new_wt_bins[i] = old_wt_bins[i]
        # } else {
        #     new_wt_bins[i] = old_wt_bins[i-1]
        # }
        if(old_wt_bins[i]>=44){
            new_wt_bins[i] = floor(old_wt_bins[i]/2)*2
        }
    }

    wt_bins_dt = data.table(bin=old_wt_bins,new_bin=new_wt_bins)
    # make everything bigger than 138 the 138 bin
    wt_bins_dt[bin>=138,new_bin:=138]

    # wt_at_pop_len = lw_params(base_ini)[1]*tmp_data$lbin_vector_pop^lw_params(base_ini)[2]
    # t_dt=as.data.table(cbind(tmp_data$lbin_vector_pop,wt_at_pop_len))
    # colnames(t_dt) = c("len","wt")
    # t_dt$bin = floor(t_dt$wt)
    # t_dt[bin<1,bin:=1]

    tmp_wtcomp_init_a = tmp_wtcomp_init[bin<44]
    tmp_wtcomp_init_b = tmp_wtcomp_init[bin>=44] %>%
                        merge(.,wt_bins_dt,by="bin") %>%
                        .[,.(value=sum(value)),by=.(year,month,week,fishery,new_bin)] %>%
                        setnames(.,"new_bin","bin")
    tmp_wtcomp = rbind(tmp_wtcomp_init_a,tmp_wtcomp_init_b) %>%
                        dcast(.,year+month+week+fishery~bin) %>%
                  na.omit(.) %>%
                  as.data.frame(.)
    # tmp_wtcomp = tmp_wtcomp_init_a %>%
    #                      dcast(.,year+month+week+fishery~bin) %>%
    #                na.omit(.) %>%
    #                as.data.frame(.)
    tmp_wtcomp_a = tmp_wtcomp[,c("year","month","fishery")]
    tmp_wtcomp_a = as.data.table(tmp_wtcomp_a) %>%
                    merge(.,ts_dt,by=c("year","month")) %>%
                    .[,month:=1] %>%
                    .[,.(ts,month,fishery)] %>%
                    setnames(.,c("ts","fishery"),c("year","fleet")) %>%
                    as.data.frame(.)
    
    tmp_wtcomp_a = cbind(rep(1,nrow(tmp_wtcomp_a)),tmp_wtcomp_a)
    colnames(tmp_wtcomp_a) = c("method","year","month","fleet")
    
    tmp_wtcomp_b = tmp_wtcomp[,-c(1:4)]
    all_new_wt_bins = as.numeric(colnames(tmp_wtcomp_b))
    all_new_wt_bins[1] = -all_new_wt_bins[1] # set first bin to be inclusive of all smaller weights
    tmp_wtcomp_c = tmp_wtcomp_b
    tmp_wtcomp_c[tmp_wtcomp_c>0] = 0
    colnames(tmp_wtcomp_c) = paste0("m",colnames(tmp_wtcomp_b))
    colnames(tmp_wtcomp_b) = paste0("f",colnames(tmp_wtcomp_b))
    
    tmp_wtcomp_a$sex = 0
    tmp_wtcomp_a$part = 0
    tmp_wtcomp_a$Nsamp = rowSums(tmp_wtcomp_b)
    # tmp_wtcomp_a$Nsamp = sapply(rowSums(tmp_wtcomp_b),function(x)min(c(x,1000))) # apply MFCL like filter to maximum sample size


    # exclude males
    tmp_wtcomp = as.data.table(cbind(tmp_wtcomp_a,tmp_wtcomp_b))
    
    tmp_data$Nobs_per_method = nrow(tmp_wtcomp)
    tmp_data$sizefreq_bins_list = list(all_new_wt_bins)
    tmp_data$nbins_per_method = length(all_new_wt_bins)
    tmp_data$sizefreq_data_list = list(as.data.frame(tmp_wtcomp))

    # misc
    tmp_data$comp_tail_compression = rep(-0.0001,tmp_data$Nfleets)
    tmp_data$add_to_comp = rep(0.0001,tmp_data$Nfleets)
    tmp_data$max_combined_lbin = rep(0,tmp_data$Nfleets)    

    SS_writedat(tmp_data, outfile=file.path(dir_bet_stock_synthesis_base,"data.ss"), overwrite=TRUE)

#_____________________________________________________________________________________________________________________________
# update control
    tmp_ctl$Comments = "#C No comments"
    tmp_ctl$nseas = tmp_data$nseas
    tmp_ctl$N_areas = tmp_data$N_areas
    tmp_ctl$Nages = tmp_data$Nages
    tmp_ctl$Nsexes = tmp_data$Nsexes
    tmp_ctl$Npopbins = tmp_data$N_lbinspop
    tmp_ctl$Nfleets = tmp_data$Nfleets
    tmp_ctl$fleetnames = tmp_data$fleetnames

    # recruitment timing and distribution
    tmp_ctl$recr_dist_method = 4
    tmp_ctl$Block_Design[[1]] = rep(tmp_data$styr-1,2)

    # don't use steepness because starting at unfished
    tmp_ctl$Use_steep_init_equi = 0

    # natural mortality
    tmp_M_dt = data.table(M = m_at_age(base_rep)) %>%
            .[,age:=seq(from=1,by=1,length.out=length(m_at_age(base_rep)))]
    tmp_M = matrix(tmp_M_dt$M,nrow=1)
    # add age 0 M
    tmp_M = cbind(matrix(rep(tmp_M_dt$M[1],1),nrow=1),tmp_M)
    rownames(tmp_M) = "natM1"
    colnames(tmp_M) = paste0("Age_",c(0,tmp_M_dt$age))
    tmp_ctl$natM = as.data.frame(tmp_M)

    # growth
    tmp_laa_dt = as.data.table(mean_laa(base_rep)) %>%
                 .[,.(age,season,value)] %>%
                 .[,dec_age:=as.numeric(age)+(as.numeric(season))/4] %>%
                 .[,model_age:=dec_age*4]
    
    tmp_sdlaa_dt = as.data.table(sd_laa(base_rep)) %>%
                 .[,.(age,season,value)] %>%
                 .[,dec_age:=as.numeric(age)+(as.numeric(season))/4] %>%
                 .[,model_age:=dec_age*4]
    
    tmp_ctl$Growth_Age_for_L1 = 1
    tmp_ctl$Growth_Age_for_L2 = 40
    tmp_ctl$CV_Growth_Pattern = 2 # make CV a function of length at age; like MFCL

    # fecundity
    tmp_ctl$fecundity_option = 3 # makes fecundity a function of weight
    tmp_ctl$maturity_option = 6 # reads in mat at length
    # Interpolate maturity to population length bins
        mat_interp = approx(
            x = mfcl_bin_mid,
            y = as.vector(mat_at_length(base_par)),
            xout = tmp_data$lbin_vector_pop,
            method = "linear",
            rule = 2
        )$y
    tmp_ctl$Length_Maturity = as.data.frame(matrix(mat_interp,nrow=1,ncol=length(tmp_data$lbin_vector_pop)))
    colnames(tmp_ctl$Length_Maturity) = as.character(tmp_data$lbin_vector_pop)  

    # add MG params
    # only grab female parameters + last 2 rows
    nrow_mg_params = nrow(tmp_ctl$MG_parms)
    tmp_MG_parms = tmp_ctl$MG_parms[c(grep("_Fem_",rownames(tmp_ctl$MG_parms),fixed=TRUE),nrow_mg_params+(-1:0)),]

    # change growth L1, L2, k, sd1, sd2
    tmp_MG_parms$INIT[grep("L_at_Amin_",rownames(tmp_MG_parms),fixed=TRUE)] = tmp_laa_dt[dec_age==0.25]$value
    tmp_MG_parms$PRIOR[grep("L_at_Amin_",rownames(tmp_MG_parms),fixed=TRUE)] = tmp_laa_dt[dec_age==0.25]$value
    
    tmp_MG_parms$INIT[grep("L_at_Amax_",rownames(tmp_MG_parms),fixed=TRUE)] = tmp_laa_dt[dec_age==10]$value
    tmp_MG_parms$PRIOR[grep("L_at_Amax_",rownames(tmp_MG_parms),fixed=TRUE)] = tmp_laa_dt[dec_age==10]$value
    
    tmp_MG_parms$INIT[grep("VonBert_K_",rownames(tmp_MG_parms),fixed=TRUE)] = growth(base_par)[3,1]
    tmp_MG_parms$PRIOR[grep("VonBert_K_",rownames(tmp_MG_parms),fixed=TRUE)] = growth(base_par)[3,1]
    tmp_MG_parms$HI[grep("VonBert_K_",rownames(tmp_MG_parms),fixed=TRUE)] = 0.99
    
    tmp_MG_parms$INIT[grep("CV_young_",rownames(tmp_MG_parms),fixed=TRUE)] = tmp_sdlaa_dt[dec_age==1]$value
    tmp_MG_parms$PRIOR[grep("CV_young_",rownames(tmp_MG_parms),fixed=TRUE)] = tmp_sdlaa_dt[dec_age==1]$value
    tmp_MG_parms$HI[grep("CV_young_",rownames(tmp_MG_parms),fixed=TRUE)] = 20
    
    tmp_MG_parms$INIT[grep("CV_old_",rownames(tmp_MG_parms),fixed=TRUE)] = tmp_sdlaa_dt[dec_age==10]$value
    tmp_MG_parms$PRIOR[grep("CV_old_",rownames(tmp_MG_parms),fixed=TRUE)] = tmp_sdlaa_dt[dec_age==10]$value
    tmp_MG_parms$HI[grep("CV_old_",rownames(tmp_MG_parms),fixed=TRUE)] = 20

    # change wtlen1 & wtlen2
    tmp_MG_parms$INIT[grep("Wtlen_1_",rownames(tmp_MG_parms),fixed=TRUE)] = lw_params(base_ini)[1]
    tmp_MG_parms$PRIOR[grep("Wtlen_1_",rownames(tmp_MG_parms),fixed=TRUE)] = lw_params(base_ini)[1]
    
    tmp_MG_parms$INIT[grep("Wtlen_2_",rownames(tmp_MG_parms),fixed=TRUE)] = lw_params(base_ini)[2]
    tmp_MG_parms$PRIOR[grep("Wtlen_2_",rownames(tmp_MG_parms),fixed=TRUE)] = lw_params(base_ini)[2]

    # mat 50% & mat slope
    tmp_mat_at_length = as.data.frame(cbind(len = mfcl_bin_mid, mat = mat_at_length(base_par)))

    fit_logistic = nls(mat ~ 1/(1 + exp(mat_slope*(len-mat50))),
                 data = tmp_mat_at_length,
                 start = list(mat_slope = -1, mat50 = 180))
    

    tmp_MG_parms$INIT[grep("Mat50%_Fem_",rownames(tmp_MG_parms),fixed=TRUE)] = summary(fit_logistic)$parameters["mat50","Estimate"]
    tmp_MG_parms$PRIOR[grep("Mat50%_Fem_",rownames(tmp_MG_parms),fixed=TRUE)] = summary(fit_logistic)$parameters["mat50","Estimate"]
    
    tmp_MG_parms$INIT[grep("Mat_slope_",rownames(tmp_MG_parms),fixed=TRUE)] = summary(fit_logistic)$parameters["mat_slope","Estimate"]
    tmp_MG_parms$PRIOR[grep("Mat_slope_",rownames(tmp_MG_parms),fixed=TRUE)] = summary(fit_logistic)$parameters["mat_slope","Estimate"]
    
    # fec a & b; set both to 1
    tmp_MG_parms$INIT[grep("Eggs_alpha_",rownames(tmp_MG_parms),fixed=TRUE)] = 1
    tmp_MG_parms$PRIOR[grep("Eggs_alpha_",rownames(tmp_MG_parms),fixed=TRUE)] = 1
    
    tmp_MG_parms$INIT[grep("Eggs_beta_",rownames(tmp_MG_parms),fixed=TRUE)] = 1
    tmp_MG_parms$PRIOR[grep("Eggs_beta_",rownames(tmp_MG_parms),fixed=TRUE)] = 1
    
    tmp_ctl$MG_parms = tmp_MG_parms

    # update virgin recruitment
    tmp_ctl$SR_parms["SR_LN(R0)","INIT"] = 15
    tmp_ctl$SR_parms["SR_LN(R0)","PRIOR"] = 15
    tmp_ctl$SR_parms["SR_LN(R0)","HI"] = 20

    # update steepness
    tmp_ctl$SR_parms["SR_BH_steep","INIT"] = 0.8
    tmp_ctl$SR_parms["SR_BH_steep","PRIOR"] = 0.8

    # update sigmaR
    tmp_ctl$SR_parms["SR_sigmaR","INIT"] = 0.6
    tmp_ctl$SR_parms["SR_sigmaR","PRIOR"] = 0.6

    # recruitment deviation setup
    tmp_ctl$MainRdevYrFirst = tmp_data$styr + 80 # start main rec devs in 1972 but could be later (1992)
    tmp_ctl$MainRdevYrLast = tmp_data$endyr - 1 # late dev for terminal year
    tmp_ctl$recdev_early_start = -40
    tmp_ctl$last_early_yr_nobias_adj = tmp_data$styr + 40 # initial values
    tmp_ctl$first_yr_fullbias_adj = tmp_data$endyr - 16 # initial values
    tmp_ctl$last_yr_fullbias_adj = tmp_data$endyr - 1 # initial values
    tmp_ctl$first_recent_yr_nobias_adj = tmp_data$endyr - 1
    tmp_ctl$max_bias_adj = 0.6

    # remove init F
    tmp_ctl$init_F = NULL
    tmp_ctl$F_Method = 3
    tmp_ctl$F_iter = 3

    # catchability options; define for each survey
    tmp_ctl$Q_options = data.frame(fleet=15,link=1,link_info=0,extra_se=0,biasadj=0,float=1)
    rownames(tmp_ctl$Q_options) = "S01_INDEX"

    tmp_ctl$Q_parms = tmp_ctl$Q_parms[1,]
    rownames(tmp_ctl$Q_parms) = "LnQ_base_S01_INDEX"

    # size based selectivity
    tmp_size_selex_types = as.data.frame(matrix(0,nrow=length(tmp_ctl$fleetnames),ncol=4))
    rownames(tmp_size_selex_types) = tmp_ctl$fleetnames
    colnames(tmp_size_selex_types) = colnames(tmp_ctl$size_selex_types)
    tmp_size_selex_types$Pattern = 24
    tmp_size_selex_types$Pattern[c(11,15)] = 1
    tmp_ctl$size_selex_types = tmp_size_selex_types

    # age based selectivity
    tmp_age_selex_types = as.data.frame(matrix(0,nrow=length(tmp_ctl$fleetnames),ncol=4))
    rownames(tmp_age_selex_types) = tmp_ctl$fleetnames
    colnames(tmp_age_selex_types) = colnames(tmp_ctl$age_selex_types)
    tmp_age_selex_types$Pattern = 11
    tmp_ctl$age_selex_types = tmp_age_selex_types

    # define selex shapes
    # LO, HI, INIT, PHASE 
    tmp_size_selex_parms = rbind(
    # Fleet 1-7: longline fleets
    make_size_selex_par_24(peak=c(10.1,200,100,2),
                          top_logit=c(-7,7,-5,-3),
                          ascend_se=c(-7,7,0,4),
                          descend_se=c(-7,7,0,5),
                          start_logit=c( -999.0,9.0,-9,-4),
                          end_logit=c( -999.0,9.0,9,6)),
    make_size_selex_par_24(peak=c(10.1,200,100,2),
                          top_logit=c(-7,7,-5,-3),
                          ascend_se=c(-7,7,0,4),
                          descend_se=c(-7,7,0,5),
                          start_logit=c( -999.0,9.0,-9,-4),
                          end_logit=c( -999.0,9.0,9,6)),
    make_size_selex_par_24(peak=c(10.1,200,100,2),
                          top_logit=c(-7,7,-5,-3),
                          ascend_se=c(-7,7,0,4),
                          descend_se=c(-7,7,0,5),
                          start_logit=c( -999.0,9.0,-9,-4),
                          end_logit=c( -999.0,9.0,9,6)),
    make_size_selex_par_24(peak=c(10.1,200,100,2),
                          top_logit=c(-7,7,-5,-3),
                          ascend_se=c(-7,7,0,4),
                          descend_se=c(-7,7,0,5),
                          start_logit=c( -999.0,9.0,-9,-4),
                          end_logit=c( -999.0,9.0,9,6)),
    make_size_selex_par_24(peak=c(10.1,200,100,2),
                          top_logit=c(-7,7,-5,-3),
                          ascend_se=c(-7,7,0,4),
                          descend_se=c(-7,7,0,5),
                          start_logit=c( -999.0,9.0,-9,-4),
                          end_logit=c( -999.0,9.0,9,6)),
    make_size_selex_par_24(peak=c(10.1,200,100,2),
                          top_logit=c(-7,7,-5,-3),
                          ascend_se=c(-7,7,0,4),
                          descend_se=c(-7,7,0,5),
                          start_logit=c( -999.0,9.0,-9,-4),
                          end_logit=c( -999.0,9.0,9,6)),
    make_size_selex_par_24(peak=c(10.1,200,100,2),
                          top_logit=c(-7,7,-5,-3),
                          ascend_se=c(-7,7,0,4),
                          descend_se=c(-7,7,0,5),
                          start_logit=c( -999.0,9.0,-9,-4),
                          end_logit=c( -999.0,9.0,9,6)),
    # Fleet 8-9: purse seine fleets
    make_size_selex_par_24(peak=c(10.1,200,40,2),
                          top_logit=c(-7,7,-5,-3),
                          ascend_se=c(-7,7,0,4),
                          descend_se=c(-7,7,0,5),
                          start_logit=c( -999.0,9.0,-9,-4),
                          end_logit=c( -999.0,9.0,-9,-6)),
    make_size_selex_par_24(peak=c(10.1,200,40,2),
                          top_logit=c(-7,7,-5,-3),
                          ascend_se=c(-7,7,0,4),
                          descend_se=c(-7,7,0,5),
                          start_logit=c( -999.0,9.0,-9,-4),
                          end_logit=c( -999.0,9.0,-9,-6)),
    # Fleet 10: domestic misc fisheries
    make_size_selex_par_24(peak=c(10.1,200,40,2),
                          top_logit=c(-7,7,-5,-3),
                          ascend_se=c(-7,7,0,4),
                          descend_se=c(-7,7,0,5),
                          start_logit=c( -999.0,9.0,0,6),
                          end_logit=c( -999.0,9.0,-9,-6)),
    # Fleet 11: domestic handline
    make_size_selex_par_1(inflection=c(5,200,100,2),
                          width=c(0,500,10,3)),
    # Fleet 12-14: PL and JP PS
    make_size_selex_par_24(peak=c(10.1,200,40,2),
                          top_logit=c(-7,7,-5,-3),
                          ascend_se=c(-7,7,0,4),
                          descend_se=c(-7,7,0,5),
                          start_logit=c( -999.0,9.0,-9,-4),
                          end_logit=c( -999.0,9.0,-9,-6)),
    make_size_selex_par_24(peak=c(10.1,200,40,2),
                          top_logit=c(-7,7,-5,-3),
                          ascend_se=c(-7,7,0,4),
                          descend_se=c(-7,7,0,5),
                          start_logit=c( -999.0,9.0,-9,-4),
                          end_logit=c( -999.0,9.0,-9,-6)),
    make_size_selex_par_24(peak=c(10.1,200,40,2),
                          top_logit=c(-7,7,-5,-3),
                          ascend_se=c(-7,7,0,4),
                          descend_se=c(-7,7,0,5),
                          start_logit=c( -999.0,9.0,-9,-4),
                          end_logit=c( -999.0,9.0,-9,-6)),
    # Fleet 15: index fishery
    make_size_selex_par_1(inflection=c(5,200,100,2),
                          width=c(0,500,10,3))
   )
    tmp_ctl$size_selex_parms = as.data.frame(tmp_size_selex_parms)

    tmp_ctl$age_selex_parms = as.data.frame(rbindlist(lapply(rep(tmp_ctl$Nages,length(tmp_ctl$fleetnames)),make_age_selex_par_11)))

    # variance adjustment list
    tmp_ctl$Variance_adjustment_list = expand.grid(Factor=1:7,Fleet=1:length(tmp_ctl$fleetnames),Value=0)
    tmp_ctl$Variance_adjustment_list$Value[which(tmp_ctl$Variance_adjustment_list$Factor%in%c(4,7))] = 1

    # adjust ess for comps
    length_var_adj = 1/rep(20000,15)
    weight_var_adj = 1/rep(20000,15)
    length_var_adj[c(1,4,5,6,15)] = 1/40000
    weight_var_adj[c(1,4,5,6,15)] = 1/40000

    # length_var_adj = 1/rep(20,15)
    # weight_var_adj = 1/rep(20,15)
    # length_var_adj[c(1,4,5,6,15)] = 1/40
    # weight_var_adj[c(1,4,5,6,15)] = 1/40

    tmp_ctl$Variance_adjustment_list$Value[which(tmp_ctl$Variance_adjustment_list$Factor%in%c(4))] = length_var_adj
    tmp_ctl$Variance_adjustment_list$Value[which(tmp_ctl$Variance_adjustment_list$Factor%in%c(7))] = weight_var_adj

    # lambdas
    tmp_lambdas_surv= data.table(like_comp=rep(1,length(tmp_ctl$fleetnames)),fleet=1:length(tmp_ctl$fleetnames)) %>%
                      .[,phase:=1] %>%
                      .[,value:=0] %>%
                      .[,sizefreq_method:=0] %>%
                      .[fleet %in% c(15),value:=1]
    tmp_lambdas_lf = data.table(like_comp=rep(4,length(tmp_ctl$fleetnames)),fleet=1:length(tmp_ctl$fleetnames)) %>%
                      .[,phase:=1] %>%
                      .[,value:=0] %>%
                      .[,sizefreq_method:=0] %>%
                      .[fleet %in% unique(tmp_lencomp$fleet),value:=1]
    tmp_lambdas_gs = data.table(like_comp=rep(6,length(unique(tmp_wtcomp$fleet))),fleet=unique(tmp_wtcomp$fleet)) %>%
                      .[,phase:=1] %>%
                      .[,value:=1] %>%
                      .[,sizefreq_method:=1]
    # tmp_ctl$lambdas = as.data.frame(tmp_lambdas_surv)
    # tmp_ctl$lambdas = as.data.frame(rbind(tmp_lambdas_surv,tmp_lambdas_lf))
    tmp_ctl$lambdas = as.data.frame(rbind(tmp_lambdas_surv,tmp_lambdas_lf,tmp_lambdas_gs))
    tmp_ctl$N_lambdas = nrow(tmp_ctl$lambdas)

    # stddev reporting
    tmp_ctl$more_stddev_reporting = 2
    tmp_ctl$stddev_reporting_specs = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0)

    SS_writectl(tmp_ctl,file.path(dir_bet_stock_synthesis_base,"control.ss"),overwrite = TRUE)

#_____________________________________________________________________________________________________________________________
# run new version of stock synthesis
    file.copy(from=paste0(proj_dir,"/executables/stock-synthesis/3.30.24.1/ss3_win.exe"),to=dir_bet_stock_synthesis_base)
    run(dir=dir_bet_stock_synthesis_base,exe="ss3_win.exe",show_in_console = TRUE,skipfinished=FALSE)
