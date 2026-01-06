#V3.30.23.1;_safe;_compile_date:_Dec  5 2024;_Stock_Synthesis_by_Richard_Methot_(NOAA)_using_ADMB_13.2
#C file created using an r4ss function
#C file write time: 2025-01-10  23:43:44
# for all year entries except rebuilder; enter either: actual year, -999 for styr, 0 for endyr, neg number for rel. endyr
1 # Benchmarks: 0=skip; 1=calc F_spr,F_btgt,F_msy; 2=calc F_spr,F0.1,F_msy; 3=add F_Blimit; 
2 # Do_MSY: 1= set to F(SPR); 2=calc F(MSY); 3=set to F(Btgt) or F0.1; 4=set to F(endyr); 5=calc F(MEY) with MSY_unit options
# if Do_MSY=5, enter MSY_Units; then list fleet_ID, cost/F, price/mt, include_in_Fmey_scaling; # -fleet_ID to fill; -9999 to terminate
0.3 # SPR target (e.g. 0.40)
0.3 # Biomass target (e.g. 0.40)
#_Bmark_years: beg_bio, end_bio, beg_selex, end_selex, beg_relF, end_relF, beg_recr_dist, end_recr_dist, beg_SRparm, end_SRparm (enter actual year, or values of 0 or -integer to be rel. endyr)
 -10 -1 -10 -1 -10 -1 0 0 0 0
#  2012 2021 2012 2021 2012 2021 2022 2022 2022 2022
# value <0 convert to endyr-value; except -999 converts to start_yr; must be >=start_yr and <=endyr
1 #Bmark_relF_Basis: 1 = use year range; 2 = set relF same as forecast below
#
0 # Forecast: -1=none; 0=simple_1yr; 1=F(SPR); 2=F(MSY) 3=F(Btgt) or F0.1; 4=Ave F (uses first-last relF yrs); 5=input annual F scalar
# where none and simple require no input after this line; simple sets forecast F same as end year F
1 # N forecast years 
1 # Fmult (only used for Do_Forecast==5) such that apical_F(f)=Fmult*relF(f)
#_Fcast_years for averaging:  beg_selex, end_selex, beg_relF, end_relF, beg_mean recruits, end_recruits  (enter actual year, or values of 0 or -integer to be rel. endyr)
 0 0 0 0 0 0
#  0 0 0 0 0 0
0 # Forecast selectivity (0=fcast selex is mean from year range; 1=fcast selectivity from time-vary parms). NOTE: logic reverses in new format
# A revised protocol for the Fcast_yr specification is available and recommended. Template is below.
#
#-12345  # code to invoke new format for expanded fcast year controls
# biology and selectivity vectors are updated annually in the forecast according to timevary parameters, so check end year of blocks and dev vectors
# input in this section directs creation of means over historical years to override any time_vary changes
# Factors implemented so far: 1=M, 4=recr_dist, 5=migration, 10=selectivity, 11=rel_F, 12=recruitment
# rel_F and Recruitment also have additional controls later in forecast.ss
# input as list: Factor, method (0, 1), st_yr, end_yr
# Terminate with -9999 for Factor
# st_yr and end_yr input can be actual year; <=0 sets rel. to timeseries endyr; Except -999 for st_yr sets to first year if time series
# Method = 0 (or omitted) continue using time_vary parms; 1  use mean of derived factor over specified year range
# Factor method st_yr end_yr 
#-9999 0 0 0
#
0 # Control rule method (0: none; 1: ramp does catch=f(SSB), buffer on F; 2: ramp does F=f(SSB), buffer on F; 3: ramp does catch=f(SSB), buffer on catch; 4: ramp does F=f(SSB), buffer on catch) 
# values for top, bottom and buffer exist, but not used when Policy=0
0.001 # Control rule inflection for constant F (as frac of Bzero, e.g. 0.40); must be > control rule cutoff, or set to -1 to use Bmsy/SSB_unf 
0.0001 # Control rule cutoff for no F (as frac of Bzero, e.g. 0.10) 
1 # Buffer:  enter Control rule target as fraction of Flimit (e.g. 0.75), negative value invokes list of [year, scalar] with filling from year to YrMax 
#
2 #_N forecast loops (1=OFL only; 2=ABC; 3=get F from forecast ABC catch with allocations applied)
1 # First forecast loop with stochastic recruitment
0 # Forecast base recruitment:  0= spawn_recr; 1=mult*spawn_recr_fxn; 2=mult*VirginRecr; 3=deprecated; 4=mult*mean_over_yr_range
# for option 4, set phase for fore_recr_devs to -1 in control to get constant mean in MCMC, else devs will be applied
1 # Value multiplier is ignored
0 # not used
#
2023  # FirstYear for caps and allocations (should be after years with fixed inputs) 
0 # stddev of log(realized catch/target catch) in forecast (set value>0.0 to cause active impl_error)
0 # Do West Coast gfish rebuilder output: 0=no; 1=yes 
2022 # Rebuilder:  first year catch could have been set to zero (Ydecl)(-1 to set to 1999)
2022 # Rebuilder:  year for current age structure (Yinit) (-1 to set to endyear+1)
1 # fleet relative F:  1=use mean over year range; 2=read seas, fleet, alloc list below
# Note that fleet allocation values is used directly as F if Do_Forecast=4 
2 # basis for fcast catch tuning and for fcast catch caps and allocation  (2=deadbio; 3=retainbio; 5=deadnum; 6=retainnum); NOTE: same units for all fleets
# Conditional input if relative F choice = 2
# enter list of:  season, fleet, relF; if used, terminate with season=-9999
# 1 1 0.00144531
# 1 2 0.000223395
# 1 3 0.00949959
# 1 5 0.00321863
# 1 6 0.0197125
# 1 7 0.0187449
# 1 8 0.00685789
# 1 9 0.0004785
# 1 10 0.148345
# 1 11 0.0382822
# 1 12 0.007172
# 1 13 0.000103256
# 1 14 0.00934297
# 2 1 0.00316748
# 2 2 0.0115051
# 2 3 0.00935086
# 2 5 0.0155025
# 2 6 0.023386
# 2 7 0.0289318
# 2 8 0.000792651
# 2 11 0.0490599
# 2 12 0.0235391
# 2 13 0.0259645
# 2 14 0.0361299
# 3 1 0.00276413
# 3 2 0.0122637
# 3 3 0.00291678
# 3 4 0.000259814
# 3 5 0.0247256
# 3 6 0.0358132
# 3 7 0.0143015
# 3 8 0.000146682
# 3 11 0.0597925
# 3 12 0.0554372
# 3 13 0.0059169
# 3 14 0.060017
# 4 1 0.00202688
# 4 5 0.0013976
# 4 6 0.0942678
# 4 7 0.010559
# 4 8 0.00262102
# 4 11 0.0463735
# 4 12 0.065616
# 4 14 0.0120278
# -9999 0 0  # terminator for list of relF
# enter list of: fleet number, max annual catch for fleets with a max; terminate with fleet=-9999
1 0
2 0
3 0
4 0
5 0
6 0
7 0
8 0
9 0
10 0
11 0
12 0
13 0
14 0
-9999 -1
# enter list of area ID and max annual catch; terminate with area=-9999
1 0
-9999 -1
# enter list of fleet number and allocation group assignment, if any; terminate with fleet=-9999
-9999 -1
#_if N allocation groups >0, list year, allocation fraction for each group 
# list sequentially because read values fill to end of N forecast
# terminate with -9999 in year field 
# no allocation groups
#
2 # basis for input Fcast catch: -1=read basis with each obs; 2=dead catch; 3=retained catch; 99=input apical_F; NOTE: bio vs num based on fleet's catchunits
#enter list of Fcast catches or Fa; terminate with line having year=-9999
#_year seas fleet catch(or_F)
-9999 1 1 0 
#
999 # verify end of input 
