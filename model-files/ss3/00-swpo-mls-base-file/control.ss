#C No comments
#C file created using an r4ss function
#C file write time: 2025-01-10  23:43:50
#
0 # 0 means do not read wtatage.ss; 1 means read and usewtatage.ss and also read and use growth parameters
1 #_N_Growth_Patterns
1 #_N_platoons_Within_GrowthPattern
4 # recr_dist_method for parameters
1 # not yet implemented; Future usage:Spawner-Recruitment; 1=global; 2=by area
1 # number of recruitment settlement assignments 
0 # unused option
# for each settlement assignment:
#_GPattern	month	area	age
1	1	1	0	#_recr_dist_pattern1
#
#_Cond 0 # N_movement_definitions goes here if N_areas > 1
#_Cond 1.0 # first age that moves (real age at begin of season, not integer) also cond on do_migration>0
#_Cond 1 1 1 2 4 10 # example move definition for seas=1, morph=1, source=1 dest=2, age1=4, age2=10
#
1 #_Nblock_Patterns
1 #_blocks_per_pattern
#_begin and end years of blocks
1951 1951
#
# controls for all timevary parameters 
1 #_env/block/dev_adjust_method for all time-vary parms (1=warn relative to base parm bounds; 3=no bound check)
#
# AUTOGEN
1 1 1 1 1 # autogen: 1st element for biology, 2nd for SR, 3rd for Q, 4th reserved, 5th for selex
# where: 0 = autogen all time-varying parms; 1 = read each time-varying parm line; 2 = read then autogen if parm min==-12345
#
# setup for M, growth, maturity, fecundity, recruitment distibution, movement
#
4 #_natM_type:_0=1Parm; 1=N_breakpoints;_2=Lorenzen;_3=agespecific;_4=agespec_withseasinterpolate;_5=Maunder_M;_6=Age-range_Lorenzen
#_ #_Age_natmort_by sex x growthpattern
#_Age_0	Age_1	Age_2	Age_3	Age_4	Age_5	Age_6	Age_7	Age_8	Age_9	Age_10
0.707867	0.4718	0.4007	0.3766	0.3669	0.3629	0.3612	0.3605	0.3602	0.36	0.36	#_natM1
1 # GrowthModel: 1=vonBert with L1&L2; 2=Richards with L1&L2; 3=age_specific_K_incr; 4=age_specific_K_decr;5=age_specific_K_each; 6=NA; 7=NA; 8=growth cessation
0.25 #_Age(post-settlement)_for_L1;linear growth below this
10 #_Growth_Age_for_L2 (999 to use as Linf)
-999 #_exponential decay for growth above maxage (value should approx initial Z; -999 replicates 3.24; -998 to not allow growth above maxage)
0 #_placeholder for future growth feature
#
0 #_SD_add_to_LAA (set to 0.1 for SS2 V1.x compatibility)
2 #_CV_Growth_Pattern:  0 CV=f(LAA); 1 CV=F(A); 2 SD=F(LAA); 3 SD=F(A); 4 logSD=F(A)
1 #_maturity_option:  1=length logistic; 2=age logistic; 3=read age-maturity matrix by growth_pattern; 4=read age-fecundity; 5=disabled; 6=read length-maturity
1 #_First_Mature_Age
3 #_fecundity option:(1)eggs=Wt*(a+b*Wt);(2)eggs=a*L^b;(3)eggs=a*Wt^b; (4)eggs=a+b*L; (5)eggs=a+b*W
0 #_hermaphroditism option:  0=none; 1=female-to-male age-specific fxn; -1=male-to-female age-specific fxn
1 #_parameter_offset_approach (1=none, 2= M, G, CV_G as offset from female-GP1, 3=like SS2 V1.x)
#
#_growth_parms
#_LO	HI	INIT	PRIOR	PR_SD	PR_type	PHASE	env_var&link	dev_link	dev_minyr	dev_maxyr	dev_PH	Block	Block_Fxn
    5	     100	    88.3822	    88.3822	 99	0	 -3	0	0	0	0	  0	0	0	#_L_at_Amin_Fem_GP_1 
   50	     600	    210.956	    210.956	 99	0	 -4	0	0	0	0	  0	0	0	#_L_at_Amax_Fem_GP_1 
 0.01	    0.99	     0.8351	     0.8351	 99	0	 -5	0	0	0	0	  0	0	0	#_VonBert_K_Fem_GP_1 
 0.01	      20	     9.9314	     9.9314	 99	0	 -2	0	0	0	0	  0	0	0	#_CV_young_Fem_GP_1  
 0.01	      20	    11.0517	    11.0517	 99	0	 -3	0	0	0	0	  0	0	0	#_CV_old_Fem_GP_1    
   -3	       3	5.39942e-07	5.39942e-07	 99	0	 -3	0	0	0	0	  0	0	0	#_Wtlen_1_Fem_GP_1   
   -3	       5	    3.58378	    3.58378	 99	0	 -3	0	0	0	0	  0	0	0	#_Wtlen_2_Fem_GP_1   
    1	     300	    181.205	    181.205	 99	0	 -3	0	0	0	0	  0	0	0	#_Mat50%_Fem_GP_1    
 -200	       3	  -0.197458	  -0.197458	 99	0	 -3	0	0	0	0	  0	0	0	#_Mat_slope_Fem_GP_1 
   -3	      20	          1	          1	0.5	6	 -3	0	0	0	0	0.5	0	0	#_Eggs_alpha_Fem_GP_1
   -3	       3	          1	          1	0.5	6	 -3	0	0	0	0	0.5	0	0	#_Eggs_beta_Fem_GP_1 
  0.1	      10	          1	          1	  1	6	 -1	0	0	0	0	  0	0	0	#_CohortGrowDev      
1e-06	0.999999	        0.5	        0.5	0.5	0	-99	0	0	0	0	  0	0	0	#_FracFemale_GP_1    
#_no timevary MG parameters
#
#_seasonal_effects_on_biology_parms
0 0 0 0 0 0 0 0 0 0 #_femwtlen1,femwtlen2,mat1,mat2,fec1,fec2,Malewtlen1,malewtlen2,L1,K
#_ LO HI INIT PRIOR PR_SD PR_type PHASE
#_Cond -2 2 0 0 -1 99 -2 #_placeholder when no seasonal MG parameters
#
3 #_Spawner-Recruitment; 2=Ricker; 3=std_B-H; 4=SCAA;5=Hockey; 6=B-H_flattop; 7=survival_3Parm;8=Shepard_3Parm
0 # 0/1 to use steepness in initial equ recruitment calculation
0 # future feature: 0/1 to make realized sigmaR a function of SR curvature
#_LO	HI	INIT	PRIOR	PR_SD	PR_type	PHASE	env-var	use_dev	dev_mnyr	dev_mxyr	dev_PH	Block	Blk_Fxn # parm_name
   2	  15	5.49073	  6	  99	0	 1	0	0	0	0	0	0	0	#_SR_LN(R0)  
 0.2	0.99	    0.8	0.8	1000	6	-2	0	0	0	0	0	0	0	#_SR_BH_steep
0.05	 1.9	    0.2	0.2	1000	6	-4	0	0	0	0	0	0	0	#_SR_sigmaR  
  -4	   4	      0	  0	  99	0	-1	0	0	0	0	0	0	0	#_SR_regime  
   0	   0	      0	  0	  99	0	-1	0	0	0	0	0	0	0	#_SR_autocorr
#_no timevary SR parameters
1 #do_recdev:  0=none; 1=devvector (R=F(SSB)+dev); 2=deviations (R=F(SSB)+dev); 3=deviations (R=R0*dev; dev2=R-f(SSB)); 4=like 3 with sum(dev2) adding penalty
1952 # first year of main recr_devs; early devs can preceed this era
2021 # last year of main recr_devs; forecast devs start in following year
1 #_recdev phase
1 # (0/1) to read 13 advanced options
1942 #_recdev_early_start (0=none; neg value makes relative to recdev_start)
2 #_recdev_early_phase
0 #_forecast_recruitment phase (incl. late recr) (0 value resets to maxphase+1)
1 #_lambda for Fcast_recr_like occurring before endyr+1
1962 #_last_yr_nobias_adj_in_MPD; begin of ramp
2006 #_first_yr_fullbias_adj_in_MPD; begin of plateau
2021 #_last_yr_fullbias_adj_in_MPD
2021 #_end_yr_for_ramp_in_MPD (can be in forecast to shape ramp, but SS sets bias_adj to 0.0 for fcast yrs)
0.5 #_max_bias_adj_in_MPD (-1 to override ramp and set biasadj=1.0 for all estimated recdevs)
0 #_period of cycles in recruitment (N parms read below)
-5 #min rec_dev
5 #max rec_dev
0 #_read_recdevs
#_end of advanced SR options
#
#_placeholder for full parameter lines for recruitment cycles
# read specified recr devs
#_Yr Input_value
#
#Fishing Mortality info
0.2 # F ballpark
-2010 # F ballpark year (neg value to disable)
3 # F_Method:  1=Pope; 2=instan. F; 3=hybrid (hybrid is recommended)
4 # max F or harvest rate, depends on F_Method
5 # N iterations for tuning F in hybrid method (recommend 3 to 7)
#
#_initial_F_parms; count = 0
#
#_Q_setup for fleets with cpue or survey data
#_fleet	link	link_info	extra_se	biasadj	float  #  fleetname
   15	1	0	1	0	1	#_S01_INDEX.1-4
-9999	0	0	0	0	0	#_terminator   
#_Q_parms(if_any);Qunits_are_ln(q)
#_LO	HI	INIT	PRIOR	PR_SD	PR_type	PHASE	env-var	use_dev	dev_mnyr	dev_mxyr	dev_PH	Block	Blk_Fxn  #  parm_name
-25	25	-5.44532	0	1	0	-1	0	0	0	0	0	0	0	#_LnQ     
  0	 1	0.190777	0	1	0	-1	0	0	0	0	0	0	0	#_extra se
#_no timevary Q parameters
#
#_size_selex_patterns
#_Pattern	Discard	Male	Special
24	0	0	0	#_1 F01_LL.JP.1   
24	0	0	0	#_2 F02_LL.JP.2   
24	0	0	0	#_3 F03_LL.JP.3   
24	0	0	0	#_4 F04_LL.JP.4   
 5	0	0	4	#_5 F05_LL.TW.4   
24	0	0	0	#_6 F06_LL.AU.2   
24	0	0	0	#_7 F07_LL.AU.3   
 5	0	0	7	#_8 F08_LL.NZ.3   
24	0	0	0	#_9 F09.REC.AU.3  
 1	0	0	0	#_10 F10_REC.NZ.3 
 5	0	0	1	#_11 F11_LL.ALL.1 
24	0	0	0	#_12 F12_LL.ALL.2 
 5	0	0	3	#_13 F13_LL.ALL.3 
 5	0	0	4	#_14 F14_LL.ALL.4 
24	0	0	0	#_15 S01_INDEX.1-4
#
#_age_selex_patterns
#_Pattern	Discard	Male	Special
11	0	0	0	#_1 F01_LL.JP.1   
11	0	0	0	#_2 F02_LL.JP.2   
11	0	0	0	#_3 F03_LL.JP.3   
11	0	0	0	#_4 F04_LL.JP.4   
11	0	0	0	#_5 F05_LL.TW.4   
11	0	0	0	#_6 F06_LL.AU.2   
11	0	0	0	#_7 F07_LL.AU.3   
11	0	0	0	#_8 F08_LL.NZ.3   
11	0	0	0	#_9 F09.REC.AU.3  
11	0	0	0	#_10 F10_REC.NZ.3 
11	0	0	0	#_11 F11_LL.ALL.1 
11	0	0	0	#_12 F12_LL.ALL.2 
11	0	0	0	#_13 F13_LL.ALL.3 
11	0	0	0	#_14 F14_LL.ALL.4 
11	0	0	0	#_15 S01_INDEX.1-4
#
#_SizeSelex
#_LO	HI	INIT	PRIOR	PR_SD	PR_type	PHASE	env-var	use_dev	dev_mnyr	dev_mxyr	dev_PH	Block	Blk_Fxn  #  parm_name
  45	325	   175.466	 100	99	0	  2	0	0	0	0	0	0	0	#_SizeSel_P_1_F01_LL.JP.1(1)   
 -10	 10	       -10	 -10	99	0	 -3	0	0	0	0	0	0	0	#_SizeSel_P_2_F01_LL.JP.1(1)   
 -10	 10	   7.25927	   0	99	0	  4	0	0	0	0	0	0	0	#_SizeSel_P_3_F01_LL.JP.1(1)   
 -10	 10	   5.14734	   0	99	0	  5	0	0	0	0	0	0	0	#_SizeSel_P_4_F01_LL.JP.1(1)   
-999	  9	      -999	-999	99	0	 -4	0	0	0	0	0	0	0	#_SizeSel_P_5_F01_LL.JP.1(1)   
-999	  9	  -1.00126	   0	99	0	  3	0	0	0	0	0	0	0	#_SizeSel_P_6_F01_LL.JP.1(1)   
  45	325	   196.432	 100	99	0	  2	0	0	0	0	0	0	0	#_SizeSel_P_1_F02_LL.JP.2(2)   
 -10	 10	  0.108229	   0	99	0	  3	0	0	0	0	0	0	0	#_SizeSel_P_2_F02_LL.JP.2(2)   
 -10	 10	   6.93494	   0	99	0	  4	0	0	0	0	0	0	0	#_SizeSel_P_3_F02_LL.JP.2(2)   
 -10	 10	0.00147782	   0	99	0	  5	0	0	0	0	0	0	0	#_SizeSel_P_4_F02_LL.JP.2(2)   
-999	  9	      -999	-999	99	0	 -4	0	0	0	0	0	0	0	#_SizeSel_P_5_F02_LL.JP.2(2)   
-999	  9	      -999	-999	99	0	 -3	0	0	0	0	0	0	0	#_SizeSel_P_6_F02_LL.JP.2(2)   
  45	325	   192.025	 100	99	0	  2	0	0	0	0	0	0	0	#_SizeSel_P_1_F03_LL.JP.3(3)   
 -10	 10	       -10	 -10	99	0	 -3	0	0	0	0	0	0	0	#_SizeSel_P_2_F03_LL.JP.3(3)   
 -10	 10	   6.14855	   0	99	0	  4	0	0	0	0	0	0	0	#_SizeSel_P_3_F03_LL.JP.3(3)   
 -10	 10	   3.95241	   0	99	0	  5	0	0	0	0	0	0	0	#_SizeSel_P_4_F03_LL.JP.3(3)   
-999	  9	      -999	-999	99	0	 -4	0	0	0	0	0	0	0	#_SizeSel_P_5_F03_LL.JP.3(3)   
-999	  9	 -0.203926	   0	99	0	  3	0	0	0	0	0	0	0	#_SizeSel_P_6_F03_LL.JP.3(3)   
  45	325	   179.581	 100	99	0	  2	0	0	0	0	0	0	0	#_SizeSel_P_1_F04_LL.JP.4(4)   
 -10	 10	       -10	 -10	99	0	 -3	0	0	0	0	0	0	0	#_SizeSel_P_2_F04_LL.JP.4(4)   
 -10	 10	   7.22933	   0	99	0	  4	0	0	0	0	0	0	0	#_SizeSel_P_3_F04_LL.JP.4(4)   
 -10	 10	   4.96125	   0	99	0	  5	0	0	0	0	0	0	0	#_SizeSel_P_4_F04_LL.JP.4(4)   
-999	  9	      -999	-999	99	0	 -4	0	0	0	0	0	0	0	#_SizeSel_P_5_F04_LL.JP.4(4)   
-999	  9	  -1.46735	   0	99	0	  3	0	0	0	0	0	0	0	#_SizeSel_P_6_F04_LL.JP.4(4)   
 -99	 10	        -1	   1	99	0	-99	0	0	0	0	0	0	0	#_SizeSel_P_1_F05_LL.TW.4(5)   
 -99	 10	        -1	   1	99	0	-99	0	0	0	0	0	0	0	#_SizeSel_P_2_F05_LL.TW.4(5)   
  45	325	   182.632	 100	99	0	  2	0	0	0	0	0	0	0	#_SizeSel_P_1_F06_LL.AU.2(6)   
 -10	 10	       -10	 -10	99	0	 -3	0	0	0	0	0	0	0	#_SizeSel_P_2_F06_LL.AU.2(6)   
 -10	 10	   5.66444	   0	99	0	  4	0	0	0	0	0	0	0	#_SizeSel_P_3_F06_LL.AU.2(6)   
 -10	 10	   5.33581	   0	99	0	  5	0	0	0	0	0	0	0	#_SizeSel_P_4_F06_LL.AU.2(6)   
-999	  9	      -999	-999	99	0	 -4	0	0	0	0	0	0	0	#_SizeSel_P_5_F06_LL.AU.2(6)   
-999	  9	  -2.14536	   0	99	0	  3	0	0	0	0	0	0	0	#_SizeSel_P_6_F06_LL.AU.2(6)   
  45	325	   187.246	 100	99	0	  2	0	0	0	0	0	0	0	#_SizeSel_P_1_F07_LL.AU.3(7)   
 -10	 10	       -10	 -10	99	0	 -3	0	0	0	0	0	0	0	#_SizeSel_P_2_F07_LL.AU.3(7)   
 -10	 10	   6.06471	   0	99	0	  4	0	0	0	0	0	0	0	#_SizeSel_P_3_F07_LL.AU.3(7)   
 -10	 10	   5.14204	   0	99	0	  5	0	0	0	0	0	0	0	#_SizeSel_P_4_F07_LL.AU.3(7)   
-999	  9	      -999	-999	99	0	 -4	0	0	0	0	0	0	0	#_SizeSel_P_5_F07_LL.AU.3(7)   
-999	  9	  -1.55445	   0	99	0	  3	0	0	0	0	0	0	0	#_SizeSel_P_6_F07_LL.AU.3(7)   
 -99	 10	        -1	   1	99	0	-99	0	0	0	0	0	0	0	#_SizeSel_P_1_F08_LL.NZ.3(8)   
 -99	 10	        -1	   1	99	0	-99	0	0	0	0	0	0	0	#_SizeSel_P_2_F08_LL.NZ.3(8)   
  45	325	   183.489	 100	99	0	  2	0	0	0	0	0	0	0	#_SizeSel_P_1_F09.REC.AU.3(9)  
 -10	 10	 -0.686316	   0	99	0	  3	0	0	0	0	0	0	0	#_SizeSel_P_2_F09.REC.AU.3(9)  
 -10	 10	   4.61448	   0	99	0	  4	0	0	0	0	0	0	0	#_SizeSel_P_3_F09.REC.AU.3(9)  
 -10	 10	  -3.09942	   0	99	0	  5	0	0	0	0	0	0	0	#_SizeSel_P_4_F09.REC.AU.3(9)  
-999	  9	      -999	-999	99	0	 -4	0	0	0	0	0	0	0	#_SizeSel_P_5_F09.REC.AU.3(9)  
-999	  9	      -999	-999	99	0	 -3	0	0	0	0	0	0	0	#_SizeSel_P_6_F09.REC.AU.3(9)  
  45	325	   187.853	 100	99	0	  2	0	0	0	0	0	0	0	#_SizeSel_P_1_F10_REC.NZ.3(10) 
   0	500	   19.5876	  10	99	0	  3	0	0	0	0	0	0	0	#_SizeSel_P_2_F10_REC.NZ.3(10) 
 -99	 10	        -1	   1	99	0	-99	0	0	0	0	0	0	0	#_SizeSel_P_1_F11_LL.ALL.1(11) 
 -99	 10	        -1	   1	99	0	-99	0	0	0	0	0	0	0	#_SizeSel_P_2_F11_LL.ALL.1(11) 
  45	325	   177.721	 100	99	0	  2	0	0	0	0	0	0	0	#_SizeSel_P_1_F12_LL.ALL.2(12) 
 -10	 10	  -9.77624	   0	99	0	  3	0	0	0	0	0	0	0	#_SizeSel_P_2_F12_LL.ALL.2(12) 
 -10	 10	   7.29868	   0	99	0	  4	0	0	0	0	0	0	0	#_SizeSel_P_3_F12_LL.ALL.2(12) 
 -10	 10	   6.68782	   0	99	0	  5	0	0	0	0	0	0	0	#_SizeSel_P_4_F12_LL.ALL.2(12) 
-999	  9	      -999	-999	99	0	 -4	0	0	0	0	0	0	0	#_SizeSel_P_5_F12_LL.ALL.2(12) 
-999	  9	      -999	-999	99	0	 -3	0	0	0	0	0	0	0	#_SizeSel_P_6_F12_LL.ALL.2(12) 
 -99	 10	        -1	   1	99	0	-99	0	0	0	0	0	0	0	#_SizeSel_P_1_F13_LL.ALL.3(13) 
 -99	 10	        -1	   1	99	0	-99	0	0	0	0	0	0	0	#_SizeSel_P_2_F13_LL.ALL.3(13) 
 -99	 10	        -1	   1	99	0	-99	0	0	0	0	0	0	0	#_SizeSel_P_1_F14_LL.ALL.4(14) 
 -99	 10	        -1	   1	99	0	-99	0	0	0	0	0	0	0	#_SizeSel_P_2_F14_LL.ALL.4(14) 
  45	325	   166.437	 100	99	0	  2	0	0	0	0	0	0	0	#_SizeSel_P_1_S01_INDEX.1-4(15)
 -10	 10	  0.455676	   0	99	0	  3	0	0	0	0	0	0	0	#_SizeSel_P_2_S01_INDEX.1-4(15)
 -10	 10	   6.25712	   0	99	0	  4	0	0	0	0	0	0	0	#_SizeSel_P_3_S01_INDEX.1-4(15)
 -10	 10	0.00392247	   0	99	0	  5	0	0	0	0	0	0	0	#_SizeSel_P_4_S01_INDEX.1-4(15)
-999	  9	      -999	-999	99	0	 -4	0	0	0	0	0	0	0	#_SizeSel_P_5_S01_INDEX.1-4(15)
-999	  9	      -999	-999	99	0	 -3	0	0	0	0	0	0	0	#_SizeSel_P_6_S01_INDEX.1-4(15)
#_AgeSelex
 0	 10	 0	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_1_F01_LL.JP.1(1)   
10	100	10	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_2_F01_LL.JP.1(1)   
 0	 10	 0	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_1_F02_LL.JP.2(2)   
10	100	10	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_2_F02_LL.JP.2(2)   
 0	 10	 0	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_1_F03_LL.JP.3(3)   
10	100	10	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_2_F03_LL.JP.3(3)   
 0	 10	 0	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_1_F04_LL.JP.4(4)   
10	100	10	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_2_F04_LL.JP.4(4)   
 0	 10	 0	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_1_F05_LL.TW.4(5)   
10	100	10	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_2_F05_LL.TW.4(5)   
 0	 10	 0	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_1_F06_LL.AU.2(6)   
10	100	10	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_2_F06_LL.AU.2(6)   
 0	 10	 0	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_1_F07_LL.AU.3(7)   
10	100	10	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_2_F07_LL.AU.3(7)   
 0	 10	 0	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_1_F08_LL.NZ.3(8)   
10	100	10	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_2_F08_LL.NZ.3(8)   
 0	 10	 0	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_1_F09.REC.AU.3(9)  
10	100	10	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_2_F09.REC.AU.3(9)  
 0	 10	 0	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_1_F10_REC.NZ.3(10) 
10	100	10	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_2_F10_REC.NZ.3(10) 
 0	 10	 0	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_1_F11_LL.ALL.1(11) 
10	100	10	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_2_F11_LL.ALL.1(11) 
 0	 10	 0	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_1_F12_LL.ALL.2(12) 
10	100	10	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_2_F12_LL.ALL.2(12) 
 0	 10	 0	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_1_F13_LL.ALL.3(13) 
10	100	10	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_2_F13_LL.ALL.3(13) 
 0	 10	 0	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_1_F14_LL.ALL.4(14) 
10	100	10	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_2_F14_LL.ALL.4(14) 
 0	 10	 0	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_1_S01_INDEX.1-4(15)
10	100	10	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_2_S01_INDEX.1-4(15)
#_no timevary selex parameters
#
0 #  use 2D_AR1 selectivity(0/1):  experimental feature
#_no 2D_AR1 selex offset used
# Tag loss and Tag reporting parameters go next
0 # TG_custom:  0=no read; 1=read if tags exist
#_Cond -6 6 1 1 2 0.01 -4 0 0 0 0 0 0 0  #_placeholder if no parameters
#
# Input variance adjustments factors: 
#_factor	fleet	value
    1	 1	        0	#_Variance_adjustment_list1  
    2	 1	        0	#_Variance_adjustment_list2  
    3	 1	        0	#_Variance_adjustment_list3  
    4	 1	      0.1	#_Variance_adjustment_list4  
    5	 1	        0	#_Variance_adjustment_list5  
    6	 1	        0	#_Variance_adjustment_list6  
    7	 1	      0.1	#_Variance_adjustment_list7  
    1	 2	        0	#_Variance_adjustment_list8  
    2	 2	        0	#_Variance_adjustment_list9  
    3	 2	        0	#_Variance_adjustment_list10 
    4	 2	0.0833333	#_Variance_adjustment_list11 
    5	 2	        0	#_Variance_adjustment_list12 
    6	 2	        0	#_Variance_adjustment_list13 
    7	 2	      0.1	#_Variance_adjustment_list14 
    1	 3	        0	#_Variance_adjustment_list15 
    2	 3	        0	#_Variance_adjustment_list16 
    3	 3	        0	#_Variance_adjustment_list17 
    4	 3	      0.1	#_Variance_adjustment_list18 
    5	 3	        0	#_Variance_adjustment_list19 
    6	 3	        0	#_Variance_adjustment_list20 
    7	 3	      0.1	#_Variance_adjustment_list21 
    1	 4	        0	#_Variance_adjustment_list22 
    2	 4	        0	#_Variance_adjustment_list23 
    3	 4	        0	#_Variance_adjustment_list24 
    4	 4	      0.1	#_Variance_adjustment_list25 
    5	 4	        0	#_Variance_adjustment_list26 
    6	 4	        0	#_Variance_adjustment_list27 
    7	 4	      0.1	#_Variance_adjustment_list28 
    1	 5	        0	#_Variance_adjustment_list29 
    2	 5	        0	#_Variance_adjustment_list30 
    3	 5	        0	#_Variance_adjustment_list31 
    4	 5	      0.1	#_Variance_adjustment_list32 
    5	 5	        0	#_Variance_adjustment_list33 
    6	 5	        0	#_Variance_adjustment_list34 
    7	 5	      0.1	#_Variance_adjustment_list35 
    1	 6	        0	#_Variance_adjustment_list36 
    2	 6	        0	#_Variance_adjustment_list37 
    3	 6	        0	#_Variance_adjustment_list38 
    4	 6	0.0555556	#_Variance_adjustment_list39 
    5	 6	        0	#_Variance_adjustment_list40 
    6	 6	        0	#_Variance_adjustment_list41 
    7	 6	      0.1	#_Variance_adjustment_list42 
    1	 7	        0	#_Variance_adjustment_list43 
    2	 7	        0	#_Variance_adjustment_list44 
    3	 7	        0	#_Variance_adjustment_list45 
    4	 7	      0.1	#_Variance_adjustment_list46 
    5	 7	        0	#_Variance_adjustment_list47 
    6	 7	        0	#_Variance_adjustment_list48 
    7	 7	      0.1	#_Variance_adjustment_list49 
    1	 8	        0	#_Variance_adjustment_list50 
    2	 8	        0	#_Variance_adjustment_list51 
    3	 8	        0	#_Variance_adjustment_list52 
    4	 8	      0.1	#_Variance_adjustment_list53 
    5	 8	        0	#_Variance_adjustment_list54 
    6	 8	        0	#_Variance_adjustment_list55 
    7	 8	      0.1	#_Variance_adjustment_list56 
    1	 9	        0	#_Variance_adjustment_list57 
    2	 9	        0	#_Variance_adjustment_list58 
    3	 9	        0	#_Variance_adjustment_list59 
    4	 9	      0.1	#_Variance_adjustment_list60 
    5	 9	        0	#_Variance_adjustment_list61 
    6	 9	        0	#_Variance_adjustment_list62 
    7	 9	0.0833333	#_Variance_adjustment_list63 
    1	10	        0	#_Variance_adjustment_list64 
    2	10	        0	#_Variance_adjustment_list65 
    3	10	        0	#_Variance_adjustment_list66 
    4	10	      0.1	#_Variance_adjustment_list67 
    5	10	        0	#_Variance_adjustment_list68 
    6	10	        0	#_Variance_adjustment_list69 
    7	10	0.0434783	#_Variance_adjustment_list70 
    1	11	        0	#_Variance_adjustment_list71 
    2	11	        0	#_Variance_adjustment_list72 
    3	11	        0	#_Variance_adjustment_list73 
    4	11	   0.0625	#_Variance_adjustment_list74 
    5	11	        0	#_Variance_adjustment_list75 
    6	11	        0	#_Variance_adjustment_list76 
    7	11	      0.1	#_Variance_adjustment_list77 
    1	12	        0	#_Variance_adjustment_list78 
    2	12	        0	#_Variance_adjustment_list79 
    3	12	        0	#_Variance_adjustment_list80 
    4	12	0.0294118	#_Variance_adjustment_list81 
    5	12	        0	#_Variance_adjustment_list82 
    6	12	        0	#_Variance_adjustment_list83 
    7	12	      0.1	#_Variance_adjustment_list84 
    1	13	        0	#_Variance_adjustment_list85 
    2	13	        0	#_Variance_adjustment_list86 
    3	13	        0	#_Variance_adjustment_list87 
    4	13	      0.1	#_Variance_adjustment_list88 
    5	13	        0	#_Variance_adjustment_list89 
    6	13	        0	#_Variance_adjustment_list90 
    7	13	      0.1	#_Variance_adjustment_list91 
    1	14	        0	#_Variance_adjustment_list92 
    2	14	        0	#_Variance_adjustment_list93 
    3	14	        0	#_Variance_adjustment_list94 
    4	14	      0.1	#_Variance_adjustment_list95 
    5	14	        0	#_Variance_adjustment_list96 
    6	14	        0	#_Variance_adjustment_list97 
    7	14	      0.1	#_Variance_adjustment_list98 
    1	15	 0.190777	#_Variance_adjustment_list99 
    2	15	        0	#_Variance_adjustment_list100
    3	15	        0	#_Variance_adjustment_list101
    4	15	0.0666667	#_Variance_adjustment_list102
    5	15	        0	#_Variance_adjustment_list103
    6	15	        0	#_Variance_adjustment_list104
    7	15	0.0666667	#_Variance_adjustment_list105
-9999	 0	        0	#_terminator                 
#
1 #_maxlambdaphase
1 #_sd_offset; must be 1 if any growthCV, sigmaR, or survey extraSD is an estimated parameter
# read 38 changes to default Lambdas (default value is 1.0)
#_like_comp	fleet	phase	value	sizefreq_method
    1	 1	1	0	1	#_Surv_F01_LL.JP.1_Phz1                        
    1	 2	1	0	1	#_Surv_F02_LL.JP.2_Phz1                        
    1	 3	1	0	1	#_Surv_F03_LL.JP.3_Phz1                        
    1	 4	1	0	1	#_Surv_F04_LL.JP.4_Phz1                        
    1	 5	1	0	1	#_Surv_F05_LL.TW.4_Phz1                        
    1	 6	1	0	1	#_Surv_F06_LL.AU.2_Phz1                        
    1	 7	1	0	1	#_Surv_F07_LL.AU.3_Phz1                        
    1	 8	1	0	1	#_Surv_F08_LL.NZ.3_Phz1                        
    1	 9	1	0	1	#_Surv_F09.REC.AU.3_Phz1                       
    1	10	1	0	1	#_Surv_F10_REC.NZ.3_Phz1                       
    1	11	1	0	1	#_Surv_F11_LL.ALL.1_Phz1                       
    1	12	1	0	1	#_Surv_F12_LL.ALL.2_Phz1                       
    1	13	1	0	1	#_Surv_F13_LL.ALL.3_Phz1                       
    1	14	1	0	1	#_Surv_F14_LL.ALL.4_Phz1                       
    1	15	1	1	1	#_Surv_S01_INDEX.1-4_Phz1                      
    4	 1	1	1	0	#_length_F01_LL.JP.1_sizefreq_method_0_Phz1    
    4	 2	1	1	0	#_length_F02_LL.JP.2_sizefreq_method_0_Phz1    
    4	 3	1	1	0	#_length_F03_LL.JP.3_sizefreq_method_0_Phz1    
    4	 4	1	0	0	#_length_F04_LL.JP.4_sizefreq_method_0_Phz1    
    4	 5	1	1	0	#_length_F05_LL.TW.4_sizefreq_method_0_Phz1    
    4	 6	1	1	0	#_length_F06_LL.AU.2_sizefreq_method_0_Phz1    
    4	 7	1	1	0	#_length_F07_LL.AU.3_sizefreq_method_0_Phz1    
    4	 8	1	1	0	#_length_F08_LL.NZ.3_sizefreq_method_0_Phz1    
    4	 9	1	0	0	#_length_F09.REC.AU.3_sizefreq_method_0_Phz1   
    4	10	1	1	0	#_length_F10_REC.NZ.3_sizefreq_method_0_Phz1   
    4	11	1	1	0	#_length_F11_LL.ALL.1_sizefreq_method_0_Phz1   
    4	12	1	1	0	#_length_F12_LL.ALL.2_sizefreq_method_0_Phz1   
    4	13	1	1	0	#_length_F13_LL.ALL.3_sizefreq_method_0_Phz1   
    4	14	1	1	0	#_length_F14_LL.ALL.4_sizefreq_method_0_Phz1   
    4	15	1	1	0	#_length_S01_INDEX.1-4_sizefreq_method_0_Phz1  
    6	 2	1	1	1	#_SizeFreq_F02_LL.JP.2_sizefreq_method_1_Phz1  
    6	10	1	1	1	#_SizeFreq_F10_REC.NZ.3_sizefreq_method_1_Phz1 
    6	 1	1	1	1	#_SizeFreq_F01_LL.JP.1_sizefreq_method_1_Phz1  
    6	 3	1	1	1	#_SizeFreq_F03_LL.JP.3_sizefreq_method_1_Phz1  
    6	15	1	1	1	#_SizeFreq_S01_INDEX.1-4_sizefreq_method_1_Phz1
    6	 9	1	1	1	#_SizeFreq_F09.REC.AU.3_sizefreq_method_1_Phz1 
    6	 6	1	1	1	#_SizeFreq_F06_LL.AU.2_sizefreq_method_1_Phz1  
    6	 7	1	1	1	#_SizeFreq_F07_LL.AU.3_sizefreq_method_1_Phz1  
-9999	 0	0	0	0	#_terminator                                   
#
0 # 0/1 read specs for more stddev reporting
#
999
