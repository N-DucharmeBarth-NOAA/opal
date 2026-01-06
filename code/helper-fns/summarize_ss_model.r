

# Nicholas Ducharme-Barth
# 2024/01/31
# Summarize model & save outputs to dir

# Copyright (c) 2024 Nicholas Ducharme-Barth
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

# dir = to_dir

weighted_mean <- function(x, w) {
  # Input validation
  if(length(x) != length(w)) {
    stop("Lengths of x and w must be equal")
  }
  if(any(w < 0)) {
    stop("Weights must be non-negative")
  }
  if(sum(w) == 0) {
    stop("Sum of weights must be positive")
  }
  
  # Calculate weighted mean: sum(x*w)/sum(w)
  total_w <- sum(w)
  weighted_sum <- sum(x * w)
  
  # Return the weighted mean
  return(weighted_sum / total_w)
}

summarize_ss_model = function(dir,verbose=TRUE,overwrite_html=FALSE)
{
    A = proc.time()
    # make summary file
    tmp_dt = data.table(id=tail(strsplit(dir,"/")[[1]],n=1),
                        number=strsplit(tail(strsplit(dir,"/")[[1]],n=1),"-")[[1]][1],
                        name=paste0(strsplit(tail(strsplit(dir,"/")[[1]],n=1),"-")[[1]][-1],collapse="-"))
    if(length(strsplit(tail(strsplit(dir,"/")[[1]],n=1),"-")[[1]])==3){
        tmp_dt$tag = strsplit(tail(strsplit(dir,"/")[[1]],n=1),"-")[[1]][3]
    } else {
        tmp_dt$tag = as.character(NA)
    }

    # check if Report.sso exists
    if(!file.exists(file.path(dir,"Report.sso"))|!file.exists(file.path(dir,"CompReport.sso"))){
        if(file.exists(file.path(dir,"Report.sso"))){
            file.remove(file.path(dir,"Report.sso"))
        }
        if(file.exists(file.path(dir,"CompReport.sso"))){
            file.remove(file.path(dir,"CompReport.sso"))
        }
        make_ss_output(dir)
    }

    exe_path = readLines(file.path(dir,"executable.txt"))
    exec_version = strsplit(exe_path,"/")[[1]][1]

    # read report file
    tmp_report = SS_output(dir,verbose = FALSE,printstats = FALSE)
    tmp_flt = tmp_report$definitions

    # grab likelihood table
    if(exec_version=="3.24U"){
        tmp_lik = as.data.table(tmp_report$likelihoods_used) %>%
              .[,id:=rep(tmp_dt$id,.N)] %>%
              .[,Label:=rownames(tmp_report$likelihoods_used)] %>%
              .[,.(id,Label,values,lambdas)]
        tmp_lik = rbind(tmp_lik[1:7],data.table(id=tmp_dt$id,Label="InitEQ_Regime",values=0,lambdas=1),tmp_lik[8:12]) 
    } else {
        tmp_lik = as.data.table(tmp_report$likelihoods_used) %>%
              .[,id:=rep(tmp_dt$id,.N)] %>%
              .[,Label:=rownames(tmp_report$likelihoods_used)] %>%
              .[,.(id,Label,values,lambdas)]
    }

    tmp_lik_flt = as.data.table(tmp_report$likelihoods_by_fleet) %>%
              .[,id:=rep(tmp_dt$id,.N)]

    # grab parameter table
    tmp_est_par = as.data.table(tmp_report$estimated_non_dev_parameters) %>%
              .[,Label:=rownames(tmp_report$estimated_non_dev_parameters)] %>%
              .[,id:=rep(tmp_dt$id,.N)] %>%
              .[,.(id,Label,Value,Parm_StDev,Phase,Min,Max,Init,Status,Afterbound)]
    tmp_rec_dev = as.data.table(tmp_report$recruitpars) %>%
              .[,Label:=rownames(tmp_report$recruitpars)] %>%
              .[,id:=rep(tmp_dt$id,.N)] %>%
              .[,.(id,Label,Value,Parm_StDev,type,Yr)]
    tmp_all_par = as.data.table(tmp_report$parameters) %>%
              .[,id:=rep(tmp_dt$id,.N)] %>%
              .[,.(id,Label,Value,Parm_StDev,Active_Cnt,Phase,Min,Max,Init,Status)]
                
    # grab derived quantities
    tmp_quants = as.data.table(tmp_report$derived_quants) %>%
                 .[,split_1:=sapply(Label,function(x)strsplit(x,"_")[[1]][1])] %>%
                 .[,split_2:=sapply(Label,function(x)strsplit(x,"_")[[1]][2])]
        # ssb
        tmp_ssb = tmp_quants[split_1=="SSB"&split_2%in%(tmp_report$startyr:tmp_report$endyr)] %>%
                  .[,id:=rep(tmp_dt$id,.N)] %>%
                  .[,yr:=as.numeric(split_2)] %>%
                  .[,type:="ssb"] %>%
                  .[,.(id,type,yr,Value,StdDev)]
        # f
        tmp_f = tmp_quants[split_1=="F"&split_2%in%(tmp_report$startyr:tmp_report$endyr)] %>%
                  .[,id:=rep(tmp_dt$id,.N)] %>%
                  .[,yr:=as.numeric(split_2)] %>%
                  .[,type:="f"] %>%
                  .[,.(id,type,yr,Value,StdDev)]
        # rec
        tmp_rec = tmp_quants[split_1=="Recr"&split_2%in%(tmp_report$startyr:tmp_report$endyr)] %>%
                  .[,id:=rep(tmp_dt$id,.N)] %>%
                  .[,yr:=as.numeric(split_2)] %>%
                  .[,type:="rec"] %>%
                  .[,.(id,type,yr,Value,StdDev)]
    
    if(exec_version=="3.24U"){
        tmp_kobe = as.data.table(tmp_report$Kobe) %>%
               .[,id:=rep(tmp_dt$id,.N)] %>%
               .[,.(id,Year,B.Bmsy,F.Fmsy)] 
    } else {
        tmp_kobe = as.data.table(tmp_report$Kobe) %>%
               .[,id:=rep(tmp_dt$id,.N)] %>%
               setnames(.,"Yr","Year") %>%
               .[,.(id,Year,B.Bmsy,F.Fmsy)]
    }
    
    tmp_dep = as.data.table(tmp_report$Dynamic_Bzero)%>%
               .[,id:=rep(tmp_dt$id,.N)] %>%
               .[,.(id,Yr,Era,SSB,SSB_nofishing)]
    tmp_srr = as.data.table(tmp_report$recruit) %>%
              .[,id:=rep(tmp_dt$id,.N)] %>%
              .[,.(id,Yr,SpawnBio,bias_adjusted,pred_recr)]

    # grab bio: M & growth
    if(exec_version=="3.24U"){
        tmp_bio = as.data.table(tmp_report$Growth_Parameters) %>%
              .[,id:=rep(tmp_dt$id,.N)] %>%
              setnames(.,"Gender","Sex")
              .[,.(id,Count,Yr,Sex,Platoon,A1,A2,L_a_A1,L_a_A2,K,A_a_L0,Linf,CVmin,CVmax,natM_amin,natM_max,M_age0,M_nages,WtLen1,WtLen2,Mat1,Mat2,Fec1,Fec2)]
    } else {
        tmp_bio = as.data.table(tmp_report$Growth_Parameters) %>%
              .[,id:=rep(tmp_dt$id,.N)] %>%
              .[,.(id,Count,Yr,Sex,Platoon,A1,A2,L_a_A1,L_a_A2,K,A_a_L0,Linf,CVmin,CVmax,natM_amin,natM_max,M_age0,M_nages,WtLen1,WtLen2,Mat1,Mat2,Fec1,Fec2)]
    
    }

    # grab selex
    tmp_len_selex = as.data.table(tmp_report$sizeselex) %>%
                    .[Factor=="Lsel"] %>%
                    .[,Factor:=NULL] %>%
                    .[,Label:=NULL] %>%
                    melt(.,id.vars=c("Fleet","Yr","Sex")) %>%
                    .[Yr==tmp_report$endyr] %>%
                    .[,id:=rep(tmp_dt$id,.N)] %>%
                    .[,Fleet_name:=tmp_flt$Fleet_name[Fleet]] %>%
                    .[,.(id,Fleet,Fleet_name,Yr,Sex,variable,value)]

    tmp_age_selex = as.data.table(tmp_report$ageselex) %>%
                    .[Factor=="Asel"] %>%
                    .[,Factor:=NULL] %>%
                    .[,Label:=NULL] %>%
                    melt(.,id.vars=c("Fleet","Yr","Seas","Sex","Morph")) %>%
                    .[Yr==tmp_report$endyr] %>%
                    .[,id:=rep(tmp_dt$id,.N)] %>%
                    .[,Fleet_name:=tmp_flt$Fleet_name[Fleet]] %>%
                    .[,.(id,Fleet,Fleet_name,Yr,Sex,variable,value)]

    # grab cpue fit
    tmp_cpue = as.data.table(tmp_report$cpue) %>% .[,.(Fleet,Fleet_name,Time,Obs,Exp,SE,Dev,Use)] %>%
               .[,id:=rep(tmp_dt$id,.N)] %>%
               .[,.(id,Fleet,Fleet_name,Time,Obs,Exp,SE,Dev,Use)]

    # grab lf fit
    tmp_len = as.data.table(tmp_report$lendbase) %>% .[,.(Yr.S,Fleet,Used,Kind,Sex,Bin,Obs,Exp,effN,Nsamp_in,Nsamp_adj)] %>%
              .[,.(Obs=sum(Obs),Exp=sum(Exp),effN=sum(effN),Nsamp_in=sum(Nsamp_in),Nsamp_adj=sum(Nsamp_adj)),by=.(Fleet,Used,Kind,Sex,Bin)] %>%
              .[,Dev:=Obs-Exp] %>%
              .[,id:=rep(tmp_dt$id,.N)] %>%
              .[,Fleet_name:=tmp_flt$Fleet_name[Fleet]] %>%
              .[,.(id,Fleet,Fleet_name,Used,Kind,Sex,Bin,Obs,Exp,Dev,effN,Nsamp_in,Nsamp_adj)]

    tmp_len_obs_time = as.data.table(tmp_report$lendbase) %>% .[,.(Yr.S,Fleet,Used,Kind,Sex,Bin,Obs)] %>%
              .[,id:=rep(tmp_dt$id,.N)] %>%
              .[,Fleet_name:=tmp_flt$Fleet_name[Fleet]] %>%
              .[,.(id,Fleet,Fleet_name,Yr.S,Used,Kind,Sex,Bin,Obs)]
    
    tmp_len_exp_time = as.data.table(tmp_report$lendbase) %>% .[,.(Yr.S,Fleet,Used,Kind,Sex,Bin,Exp)] %>%
              .[,id:=rep(tmp_dt$id,.N)] %>%
              .[,Fleet_name:=tmp_flt$Fleet_name[Fleet]] %>%
              .[,.(id,Fleet,Fleet_name,Yr.S,Used,Kind,Sex,Bin,Exp)] %>%
              .[,.(Exp_Mean=weighted_mean(x=Bin,w=Exp)),by=.(id,Fleet,Fleet_name,Yr.S,Used,Kind,Sex)]
    
    tmp_size = as.data.table(tmp_report$sizedbase) %>% .[,.(Yr.S,Fleet,Used,Kind,Sex,Bin,Obs,Exp,effN,Nsamp_in,Nsamp_adj)] %>%
              .[,.(Obs=sum(Obs),Exp=sum(Exp),effN=sum(effN),Nsamp_in=sum(Nsamp_in),Nsamp_adj=sum(Nsamp_adj)),by=.(Fleet,Used,Kind,Sex,Bin)] %>%
              .[,Dev:=Obs-Exp] %>%
              .[,id:=rep(tmp_dt$id,.N)] %>%
              .[,Fleet_name:=tmp_flt$Fleet_name[Fleet]] %>%
              .[,.(id,Fleet,Fleet_name,Used,Kind,Sex,Bin,Obs,Exp,Dev,effN,Nsamp_in,Nsamp_adj)]

    tmp_size_obs_time = as.data.table(tmp_report$sizedbase) %>% .[,.(Yr.S,Fleet,Used,Kind,Sex,Bin,Obs)] %>%
              .[,id:=rep(tmp_dt$id,.N)] %>%
              .[,Fleet_name:=tmp_flt$Fleet_name[Fleet]] %>%
              .[,.(id,Fleet,Fleet_name,Yr.S,Used,Kind,Sex,Bin,Obs)]
    
    tmp_size_exp_time = as.data.table(tmp_report$sizedbase) %>% .[,.(Yr.S,Fleet,Used,Kind,Sex,Bin,Exp)] %>%
              .[,id:=rep(tmp_dt$id,.N)] %>%
              .[,Fleet_name:=tmp_flt$Fleet_name[Fleet]] %>%
              .[,.(id,Fleet,Fleet_name,Yr.S,Used,Kind,Sex,Bin,Exp)] %>%
              .[,.(Exp_Mean=weighted_mean(x=Bin,w=Exp)),by=.(id,Fleet,Fleet_name,Yr.S,Used,Kind,Sex)]

    
    # tmp_age = as.data.table(tmp_report$agedbase) %>% .[,.(Yr.S,Fleet,Used,Kind,Sex,Bin,Obs,Exp,effN)] %>%
    #           .[,.(Obs=sum(Obs),Exp=sum(Exp),effN=sum(effN)),by=.(Fleet,Used,Kind,Sex,Bin)] %>%
    #           .[,Dev:=Obs-Exp] %>%
    #           .[,id:=rep(tmp_dt$id,.N)] %>%
    #           .[,Fleet_name:=as.vector(as.matrix(tmp_flt[2,]))[-1][Fleet]] %>%
    #           .[,.(id,Fleet,Fleet_name,Used,Kind,Sex,Bin,Obs,Exp,Dev,effN)]

    # tmp_con = as.data.table(tmp_report$condbase) %>% .[,.(Yr.S,Fleet,Used,Kind,Sex,Bin,Obs,Exp,effN)] %>%
    #           .[,.(Obs=sum(Obs),Exp=sum(Exp),effN=sum(effN)),by=.(Fleet,Used,Kind,Sex,Bin)] %>%
    #           .[,Dev:=Obs-Exp] %>%
    #           .[,id:=rep(tmp_dt$id,.N)] %>%
    #           .[,Fleet_name:=as.vector(as.matrix(tmp_flt[2,]))[-1][Fleet]] %>%
    #           .[,.(id,Fleet,Fleet_name,Used,Kind,Sex,Bin,Obs,Exp,Dev,effN)]

    # populate summary file
    tmp_dt$n_warn = tmp_report$Nwarnings
    tmp_dt$par_status = mean(tmp_est_par$Status=="OK")
    tmp_dt$par_afterbound = mean(tmp_est_par$Afterbound=="OK")
    tmp_dt$nfleets = tmp_report$nfleets
    tmp_dt$startyr = tmp_report$startyr
    tmp_dt$endyr = tmp_report$endyr
    tmp_dt$n_seasons = tmp_report$nseasons
    tmp_dt$ss_version = tmp_report$SS_versionNumeric
    tmp_dt$n_par_active = tmp_report$N_estimated_parameters
    tmp_dt$n_par_active_nd = nrow(tmp_report$estimated_non_dev_parameters)
    tmp_dt$n_par_total = nrow(tmp_report$parameters)
    tmp_dt$log_det_hessian = tmp_report$log_det_hessian
    tmp_dt$mgc = tmp_report$maximum_gradient_component
    tmp_dt$nll = tmp_lik[Label=="TOTAL"]$values
    
    tmp_dt = cbind(tmp_dt,matrix(data=tmp_lik[Label!="TOTAL",.(values,Label)]$values,nrow=1,dimnames=list(1,tmp_lik[Label!="TOTAL",.(values,Label)]$Label)))

    # trend in recruitment
    set.seed(123)
    if(length(tmp_rec_dev[type=="Main_RecrDev"]$Value)>0){
        tmp_ltrend = try(funtimes::notrend_test(tmp_rec_dev[type=="Main_RecrDev"]$Value)$p.value,silent=TRUE)
        if(class(tmp_ltrend)=="try-error"){
            tmp_dt$ltrend = NA
        } else {
            tmp_dt$ltrend = tmp_ltrend
        }
        tmp_dt$runs = randtests::runs.test(tmp_rec_dev[type=="Main_RecrDev"]$Value, threshold = 0, alternative = "two.sided")$p.value
        tmp_dt$ar = lmtest::dwtest(tmp_rec_dev[type=="Main_RecrDev"]$Value ~ 1, order.by = 1:length(tmp_rec_dev[type=="Main_RecrDev"]$Value), alternative = "two.sided")$p.value
        tmp_dt$sdr = sd(tmp_rec_dev[type=="Main_RecrDev"]$Value)
    } else {
        tmp_dt$ltrend = NA
        tmp_dt$runs = NA
        tmp_dt$ar = NA
        tmp_dt$sdr = NA
    }
    
    tmp_dt$sigmaR = tmp_report$sigma_R_in
    tmp_dt$steepness = tmp_all_par[Label=="SR_BH_steep"]$Value
    if(exec_version=="3.24U"){
        tmp_dt$M_m = mean(as.data.table(tmp_report$Growth_Parameters)[Gender==2]$M_nages)
        tmp_dt$L0_m = mean(as.data.table(tmp_report$Growth_Parameters)[Gender==2]$L_a_A1)
        tmp_dt$Linf_m = mean(as.data.table(tmp_report$Growth_Parameters)[Gender==2]$Linf)
        tmp_dt$k_m = mean(as.data.table(tmp_report$Growth_Parameters)[Gender==1]$K)
        tmp_dt$M_f = mean(as.data.table(tmp_report$Growth_Parameters)[Gender==1]$M_nages)
        tmp_dt$L0_f = mean(as.data.table(tmp_report$Growth_Parameters)[Gender==1]$L_a_A1)
        tmp_dt$Linf_f = mean(as.data.table(tmp_report$Growth_Parameters)[Gender==1]$Linf)
        tmp_dt$k_f = mean(as.data.table(tmp_report$Growth_Parameters)[Gender==1]$K)
    } else{
        tmp_dt$M_m = mean(as.data.table(tmp_report$Growth_Parameters)[Sex==2]$M_nages)
        tmp_dt$L0_m = mean(as.data.table(tmp_report$Growth_Parameters)[Sex==2]$L_a_A1)
        tmp_dt$Linf_m = mean(as.data.table(tmp_report$Growth_Parameters)[Sex==2]$Linf)
        tmp_dt$k_m = mean(as.data.table(tmp_report$Growth_Parameters)[Sex==1]$K)
        tmp_dt$M_f = mean(as.data.table(tmp_report$Growth_Parameters)[Sex==1]$M_nages)
        tmp_dt$L0_f = mean(as.data.table(tmp_report$Growth_Parameters)[Sex==1]$L_a_A1)
        tmp_dt$Linf_f = mean(as.data.table(tmp_report$Growth_Parameters)[Sex==1]$Linf)
        tmp_dt$k_f = mean(as.data.table(tmp_report$Growth_Parameters)[Sex==1]$K)
    }
    
    tmp_dt$F_2016 = mean(tmp_quants[Label=="F_2016"]$Value)
    tmp_dt$F_2022 = mean(tmp_quants[Label=="F_2022"]$Value)
    tmp_dt$SSB_virgin = mean(tmp_quants[Label=="SSB_Virgin"]$Value)
    tmp_dt$SSB_Initial = mean(tmp_quants[Label=="SSB_Initial"]$Value)
    tmp_dt$SSB_2016 = mean(tmp_quants[Label=="SSB_2016"]$Value)
    tmp_dt$SSB_2022 = mean(tmp_quants[Label=="SSB_2022"]$Value)
    if(exec_version=="3.24U"){
        tmp_dt$MSY = mean(tmp_quants[Label=="TotYield_MSY"]$Value)
    } else {
        tmp_dt$MSY = mean(tmp_quants[Label=="Dead_Catch_MSY"]$Value)
    }
    tmp_dt$SSB_MSY = mean(tmp_quants[Label=="SSB_MSY"]$Value)
    if(exec_version=="3.24U"){
        tmp_dt$F_MSY = mean(tmp_quants[Label=="Fstd_MSY"]$Value)
    } else {
        tmp_dt$F_MSY = mean(tmp_quants[Label=="annF_MSY"]$Value)
    }
    tmp_dt$B_Bmsy_2016 = mean(tmp_kobe[Year==2016]$B.Bmsy)
    tmp_dt$B_Bmsy_2022 = mean(tmp_kobe[Year==2022]$B.Bmsy)
    tmp_dt$F_Fmsy_2016 = mean(tmp_kobe[Year==2016]$F.Fmsy)
    tmp_dt$F_Fmsy_2022 = mean(tmp_kobe[Year==2022]$F.Fmsy)
    tmp_dt$R_virgin = mean(tmp_quants[Label=="Recr_Virgin"]$Value)
    tmp_dt$R_initial = mean(tmp_quants[Label=="Recr_Initial"]$Value)
    tmp_dt$dep_st = tmp_report$current_depletion
    tmp_dt$dep_dyn_2016 = mean(tmp_dep[Yr==2016,.(SSB/SSB_nofishing),by=Yr]$V1)
    tmp_dt$dep_dyn_2022 = mean(tmp_dep[Yr==2022,.(SSB/SSB_nofishing),by=Yr]$V1)

    if(file.exists(file.path(dir,"retro.csv"))){
        stop("Functionality for summarizing retrospectives not yet built.")
    } else {
        tmp_dt$mohns_rho_ssb = NA
        tmp_dt$mohns_rho_f = NA
    }
    if(file.exists(file.path(dir,"hindcast.csv"))){
        stop("Functionality for summarizing hindcasting not yet built.")
    } else {
        tmp_dt$average_xval_cpue_mase = NA
        tmp_dt$average_xval_len_mase = NA
    }

    # extra quantities
    tmp_catch = as.data.table(tmp_report$catch) %>%
                .[Exp>0,.(Fleet,Time,Obs,Exp,sel_bio,sel_num,F,Like,Mult)] %>%
                .[,id:=rep(tmp_dt$id,.N)] %>%
                .[,.(id,Fleet,Time,Obs,Exp,sel_bio,sel_num,F,Like,Mult)]
    
    if("sex_ratio_area:_1" %in% colnames(tmp_report$endgrowth)){
        tmp_bio_age = as.data.table(tmp_report$endgrowth) %>%
                  setnames(.,c("Sex","Real_Age","M","Len_Beg","SD_Beg","Wt_Beg","Len_Mat","Mat*Fecund","sex_ratio_area:_1"),c("sex","age","m","len","len_sd","wt","mat","mat_fec","sr")) %>%
                  .[,id:=rep(tmp_dt$id,.N)] %>%
                  .[,.(id,sex,age,m,len,len_sd,wt,mat,mat_fec,sr)]
    } else {
        tmp_bio_age = as.data.table(tmp_report$endgrowth) %>%
                  setnames(.,c("Sex","Real_Age","M","Len_Beg","SD_Beg","Wt_Beg","Len_Mat","Mat*Fecund"),c("sex","age","m","len","len_sd","wt","mat","mat_fec")) %>%
                  .[,id:=rep(tmp_dt$id,.N)] %>%
                  .[,sr:=NA] %>%
                  .[,.(id,sex,age,m,len,len_sd,wt,mat,mat_fec,sr)]
    }

    tmp_mature_n = as.data.table(tmp_report$natage) %>%
                    setnames(.,"Beg/Mid","B_M") %>%
                    .[,Area:=NULL] %>%
                    .[,Bio_Pattern:=NULL] %>%
                    .[Sex==1] %>%
                    .[,Sex:=NULL] %>%
                    .[,BirthSeas:=NULL] %>%
                    .[,Settlement:=NULL] %>%
                    .[,Platoon:=NULL] %>%
                    .[,Morph:=NULL] %>%
                    .[,Yr:=NULL] %>%
                    .[,Seas:=NULL] %>%
                    .[B_M=="B"] %>%
                    .[,B_M:=NULL] %>%
                    .[,id:=rep(tmp_dt$id,.N)] %>%
                   melt(.,id.vars=c("id","Time","Era")) %>%
                   setnames(.,c("Time","Era","variable","value"),c("time","era","age","n")) %>%
                   .[,age:=as.integer(as.character(age))] %>%
                   merge(.,tmp_bio_age[sex==1,.(age,mat,mat_fec)],by="age") %>%
                   .[,n_mature:=n*mat] %>%
                   .[,n_pup_produced:=n*mat_fec] %>%
                   .[,mat:=NULL] %>%
                   .[,mat_fec:=NULL] %>%
                   .[,.(n=sum(n),n_mature=sum(n_mature),n_pup_produced=sum(n_pup_produced)),by=.(id,time,era)]
    
    if(file.exists(file.path(dir,"control.ss_new")) & file.exists(file.path(dir,"data_echo.ss_new"))){
        tmp_ctl = SS_readctl(file=file.path(dir,"control.ss_new"),datlist = file.path(dir,"data_echo.ss_new"))
    } else {
        tmp_ctl = SS_readctl(file=file.path(dir,"control.ss"),datlist = file.path(dir,"data.ss"))
    }
    

    tmp_varadj = as.data.table(tmp_ctl$Variance_adjustment_list) %>%
                .[,id:=rep(tmp_dt$id,.N)] %>%
                .[(factor==1&value!=0)|(factor%in%c(4,7)&(value!=1))] %>%
                .[,.(id,factor,fleet,value)]
    
    

    if(!(1 %in% unique(tmp_varadj$factor)) & length(grep("Q_extraSD",rownames(tmp_ctl$Q_parms)))>0){
        extra_var_idx = grep("Q_extraSD",rownames(tmp_ctl$Q_parms))
        extra_var_val = extra_var_fleet = rep(NA,length(extra_var_idx))
        for(i in seq_along(extra_var_idx)){
           extra_var_fleet[i] = as.numeric(gsub(")","",strsplit(rownames(tmp_ctl$Q_parms)[extra_var_idx[i]],"(",fixed=TRUE)[[1]][2],fixed=TRUE))
           extra_var_val[i] = tmp_ctl$Q_parms$INIT[extra_var_idx[i]]
        }
        tmp_varadj = rbind(tmp_varadj,data.table(id=rep(unique(tmp_varadj$id),length(extra_var_idx)),factor=rep(1,length(extra_var_idx)),fleet=extra_var_fleet,value=extra_var_val))
        tmp_varadj = tmp_varadj[order(id,fleet,factor)]
    }
    
    tmp_lambdas = as.data.table(tmp_ctl$lambdas) %>%
                .[,id:=rep(tmp_dt$id,.N)] %>%
                .[value>0] %>%
                .[,.(id,like_comp, fleet, phase, value, sizefreq_method)]


    # write outputs
    summarize_fleets(dir,tmp_report)
    fwrite(tmp_lik,file=file.path(dir,"lik_tbl.csv"))
    fwrite(tmp_lik_flt,file=file.path(dir,"lik_flt_tbl.csv"))
    fwrite(tmp_est_par,file=file.path(dir,"est_par.csv"))
    fwrite(tmp_rec_dev,file=file.path(dir,"rec_dev.csv"))
    fwrite(tmp_ssb,file=file.path(dir,"ssb.csv"))
    fwrite(tmp_f,file=file.path(dir,"f.csv"))
    fwrite(tmp_rec,file=file.path(dir,"rec.csv"))
    fwrite(tmp_kobe,file=file.path(dir,"kobe.csv"))
    fwrite(tmp_dep,file=file.path(dir,"dyn_dep.csv"))
    fwrite(tmp_srr,file=file.path(dir,"srr.csv"))
    fwrite(tmp_bio,file=file.path(dir,"bio.csv"))
    fwrite(tmp_len_selex,file=file.path(dir,"selex_l.csv"))
    fwrite(tmp_age_selex,file=file.path(dir,"selex_a.csv"))
    fwrite(tmp_cpue,file=file.path(dir,"cpue.csv"))
    fwrite(tmp_len,file=file.path(dir,"comp_len.csv"))
    fwrite(tmp_size,file=file.path(dir,"comp_size.csv"))
    fwrite(tmp_dt,file=file.path(dir,"summary.csv"))

    fwrite(tmp_catch,file=file.path(dir,"catch.csv"))
    fwrite(tmp_bio_age,file=file.path(dir,"bio_age.csv"))
    fwrite(tmp_mature_n,file=file.path(dir,"mature_n.csv"))

    fwrite(tmp_varadj,file=file.path(dir,"varadj.csv"))
    fwrite(tmp_lambdas,file=file.path(dir,"lambdas.csv"))

    fwrite(tmp_len_obs_time,file=file.path(dir,"comp_len_obs_time.csv"))
    fwrite(tmp_len_exp_time,file=file.path(dir,"comp_len_exp_time.csv"))
    fwrite(tmp_size_obs_time,file=file.path(dir,"comp_size_obs_time.csv"))
    fwrite(tmp_size_exp_time,file=file.path(dir,"comp_size_exp_time.csv"))

    # save output for html table
    if(!file.exists(file.path(dir,"html_parameters.txt"))|overwrite_html){
        write.table(tmp_report$parameters,file=file.path(dir,"html_parameters.txt"))
    }
    if(!file.exists(file.path(dir,"html_parameters_rownames.txt"))|overwrite_html){
        writeLines(row.names(tmp_report$parameters),con=file.path(dir,"html_parameters_rownames.txt"))
    }
    if(!file.exists(file.path(dir,"html_recruitpars.txt"))|overwrite_html){
            write.table(tmp_report$recruitpars,file=file.path(dir,"html_recruitpars.txt"))
    }
    if(!file.exists(file.path(dir,"html_recruitpars_rownames.txt"))|overwrite_html){
        writeLines(row.names(tmp_report$recruitpars),con=file.path(dir,"html_recruitpars_rownames.txt"))
    }
    if(!file.exists(file.path(dir,"html_estimated_non_dev_parameters.txt"))|overwrite_html){
        write.table(tmp_report$estimated_non_dev_parameters,file=file.path(dir,"html_estimated_non_dev_parameters.txt"))
    }
    if(!file.exists(file.path(dir,"html_estimated_non_dev_parameters_rownames.txt"))|overwrite_html){
        writeLines(row.names(tmp_report$estimated_non_dev_parameters),con=file.path(dir,"html_estimated_non_dev_parameters_rownames.txt"))
    }
    if(!file.exists(file.path(dir,"html_derived_quants.txt"))|overwrite_html){
        write.table(tmp_report$derived_quants,file=file.path(dir,"html_derived_quants.txt"))
    }
    if(!file.exists(file.path(dir,"html_derived_quants_rownames.txt"))|overwrite_html){
    writeLines(row.names(tmp_report$derived_quants),con=file.path(dir,"html_derived_quants_rownames.txt"))
    }    
  
    B=proc.time()
    time = round((B-A)[3],digits=2)
    if(verbose){
        print(paste0("Model run summarized in ",time," seconds."))
        output_files = list.files(dir)[grep(".csv",list.files(dir),fixed=TRUE)]
        output_files = c(output_files,list.files(dir)[grep("html_",list.files(dir),fixed=TRUE)])
        print(paste0(length(output_files)," files produced taking up ",round(sum(unname(sapply(file.path(dir,output_files),file.size)))/1024^2,digits=2)," MB."))
    }
}
