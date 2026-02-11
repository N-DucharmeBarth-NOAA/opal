#C No comments
#C file created using an r4ss function
#C file write time: 2026-02-10  19:33:27
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
0 0
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
#_Age_0	Age_1	Age_2	Age_3	Age_4	Age_5	Age_6	Age_7	Age_8	Age_9	Age_10	Age_11	Age_12	Age_13	Age_14	Age_15	Age_16	Age_17	Age_18	Age_19	Age_20	Age_21	Age_22	Age_23	Age_24	Age_25	Age_26	Age_27	Age_28	Age_29	Age_30	Age_31	Age_32	Age_33	Age_34	Age_35	Age_36	Age_37	Age_38	Age_39	Age_40
0.2001	0.2001	0.1663	0.1337	0.1005	0.1	0.1	0.1	0.1	0.1001	0.1004	0.101	0.1026	0.1053	0.1086	0.1113	0.113	0.1138	0.1141	0.1141	0.1141	0.1139	0.1137	0.1135	0.1133	0.1131	0.1128	0.1126	0.1124	0.1122	0.1119	0.1117	0.1115	0.1113	0.111	0.1108	0.1106	0.1104	0.1102	0.1099	0.1097	#_natM1
1 # GrowthModel: 1=vonBert with L1&L2; 2=Richards with L1&L2; 3=age_specific_K_incr; 4=age_specific_K_decr;5=age_specific_K_each; 6=NA; 7=NA; 8=growth cessation
1 #_Age(post-settlement)_for_L1;linear growth below this
40 #_Growth_Age_for_L2 (999 to use as Linf)
-999 #_exponential decay for growth above maxage (value should approx initial Z; -999 replicates 3.24; -998 to not allow growth above maxage)
0 #_placeholder for future growth feature
#
0 #_SD_add_to_LAA (set to 0.1 for SS2 V1.x compatibility)
2 #_CV_Growth_Pattern:  0 CV=f(LAA); 1 CV=F(A); 2 SD=F(LAA); 3 SD=F(A); 4 logSD=F(A)
6 #_maturity_option:  1=length logistic; 2=age logistic; 3=read age-maturity matrix by growth_pattern; 4=read age-fecundity; 5=disabled; 6=read length-maturity
# Length Maturity: 
#_Len_1	Len_2	Len_3	Len_4	Len_5	Len_6	Len_7	Len_8	Len_9	Len_10	Len_11	Len_12	Len_13	Len_14	Len_15	Len_16	Len_17	Len_18	Len_19	Len_20	Len_21	Len_22	Len_23	Len_24	Len_25	Len_26	Len_27	Len_28	Len_29	Len_30	Len_31	Len_32	Len_33	Len_34	Len_35	Len_36	Len_37	Len_38	Len_39	Len_40	Len_41	Len_42	Len_43	Len_44	Len_45	Len_46	Len_47	Len_48	Len_49	Len_50	Len_51	Len_52	Len_53	Len_54	Len_55	Len_56	Len_57	Len_58	Len_59	Len_60	Len_61	Len_62	Len_63	Len_64	Len_65	Len_66	Len_67	Len_68	Len_69	Len_70	Len_71	Len_72	Len_73	Len_74	Len_75	Len_76	Len_77	Len_78	Len_79	Len_80	Len_81	Len_82	Len_83	Len_84	Len_85	Len_86	Len_87	Len_88	Len_89	Len_90	Len_91	Len_92	Len_93	Len_94	Len_95	Len_96	Len_97	Len_98	Len_99	Len_100	Len_101	Len_102	Len_103	Len_104	Len_105	Len_106	Len_107	Len_108	Len_109	Len_110	Len_111	Len_112	Len_113	Len_114	Len_115	Len_116	Len_117	Len_118	Len_119	Len_120	Len_121	Len_122	Len_123	Len_124	Len_125	Len_126	Len_127	Len_128	Len_129	Len_130	Len_131	Len_132	Len_133	Len_134	Len_135	Len_136	Len_137	Len_138	Len_139	Len_140	Len_141	Len_142	Len_143	Len_144	Len_145	Len_146	Len_147	Len_148	Len_149	Len_150	Len_151	Len_152	Len_153	Len_154	Len_155	Len_156	Len_157	Len_158	Len_159	Len_160	Len_161	Len_162	Len_163	Len_164	Len_165	Len_166	Len_167	Len_168	Len_169	Len_170	Len_171	Len_172	Len_173	Len_174	Len_175	Len_176	Len_177	Len_178	Len_179	Len_180	Len_181	Len_182	Len_183	Len_184	Len_185	Len_186	Len_187	Len_188	Len_189	Len_190	Len_191	Len_192	Len_193	Len_194	Len_195	Len_196	Len_197	Len_198	Len_199	Len_200	Len_201	Len_202	Len_203	Len_204	Len_205	Len_206	Len_207	Len_208	Len_209	Len_210	Len_211	Len_212	Len_213	Len_214	Len_215	Len_216
1.96075e-09	1.96075e-09	1.96075e-09	1.96075e-09	1.96075e-09	1.96075e-09	1.96075e-09	2.79532e-09	3.62989e-09	5.04564e-09	6.46139e-09	8.81499e-09	1.11686e-08	1.5019e-08	1.88693e-08	2.50862e-08	3.13032e-08	4.12321e-08	5.11609e-08	6.68718e-08	8.25828e-08	1.07246e-07	1.31909e-07	1.70357e-07	2.08806e-07	2.68382e-07	3.27959e-07	4.1978e-07	5.116e-07	6.52437e-07	7.93273e-07	1.00836e-06	1.22346e-06	1.55068e-06	1.87789e-06	2.37394e-06	2.86999e-06	3.61955e-06	4.36911e-06	5.49841e-06	6.62771e-06	8.32451e-06	1.00213e-05	1.25644e-05	1.51075e-05	1.89101e-05	2.27127e-05	2.83864e-05	3.40601e-05	4.25084e-05	5.09566e-05	6.35128e-05	7.60689e-05	9.46973e-05	0.000113326	0.000140917	0.000168508	0.000209306	0.000250104	0.000310322	0.00037054	0.000459265	0.00054799	0.000678489	0.000808988	0.00100059	0.00119219	0.00147293	0.00175368	0.0021641	0.00257453	0.0031729	0.00377127	0.00464074	0.0055102	0.00676822	0.00802623	0.00983643	0.0116466	0.0142326	0.0168186	0.0204777	0.0241367	0.0292486	0.0343606	0.0413835	0.0484065	0.0578466	0.0672866	0.0796269	0.0919671	0.10755	0.123132	0.14201	0.160887	0.182696	0.204505	0.228433	0.252361	0.277265	0.302168	0.32681	0.351451	0.374761	0.39807	0.419319	0.440568	0.459409	0.478251	0.494656	0.511062	0.525212	0.539362	0.551541	0.56372	0.574241	0.584761	0.593922	0.603082	0.611143	0.619204	0.626386	0.633567	0.640049	0.64653	0.652457	0.658383	0.663872	0.669362	0.674508	0.679654	0.684535	0.689415	0.694093	0.69877	0.703297	0.707823	0.712241	0.716659	0.721005	0.725351	0.729654	0.733957	0.738243	0.742529	0.746819	0.751108	0.75542	0.759731	0.764078	0.768426	0.772822	0.777218	0.781673	0.786128	0.790649	0.795171	0.799767	0.804362	0.809037	0.813711	0.818469	0.823227	0.828072	0.832917	0.83785	0.842784	0.847808	0.852833	0.857949	0.863065	0.868273	0.873482	0.878783	0.884084	0.889478	0.894871	0.900356	0.905842	0.911418	0.916995	0.922661	0.928328	0.934085	0.939841	0.945685	0.95153	0.957461	0.963393	0.96941	0.975427	0.981529	0.98763	0.993815	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	#_Length_Maturity1
1 #_First_Mature_Age
3 #_fecundity option:(1)eggs=Wt*(a+b*Wt);(2)eggs=a*L^b;(3)eggs=a*Wt^b; (4)eggs=a+b*L; (5)eggs=a+b*W
0 #_hermaphroditism option:  0=none; 1=female-to-male age-specific fxn; -1=male-to-female age-specific fxn
1 #_parameter_offset_approach (1=none, 2= M, G, CV_G as offset from female-GP1, 3=like SS2 V1.x)
#
#_growth_parms
#_LO	HI	INIT	PRIOR	PR_SD	PR_type	PHASE	env_var&link	dev_link	dev_minyr	dev_maxyr	dev_PH	Block	Block_Fxn
    5	     100	    30.9192	    30.9192	 99	0	 -3	0	0	0	0	  0	0	0	#_L_at_Amin_Fem_GP_1 
   50	     600	    153.443	    153.443	 99	0	 -4	0	0	0	0	  0	0	0	#_L_at_Amax_Fem_GP_1 
 0.01	    0.99	    0.09825	    0.09825	 99	0	 -5	0	0	0	0	  0	0	0	#_VonBert_K_Fem_GP_1 
 0.01	      20	      6.483	      6.483	 99	0	 -2	0	0	0	0	  0	0	0	#_CV_young_Fem_GP_1  
 0.01	      20	     16.497	     16.497	 99	0	 -3	0	0	0	0	  0	0	0	#_CV_old_Fem_GP_1    
   -3	       3	6.48085e-05	6.48085e-05	 99	0	 -3	0	0	0	0	  0	0	0	#_Wtlen_1_Fem_GP_1   
   -3	       5	    2.78104	    2.78104	 99	0	 -3	0	0	0	0	  0	0	0	#_Wtlen_2_Fem_GP_1   
    1	     300	    123.167	    123.167	 99	0	 -3	0	0	0	0	  0	0	0	#_Mat50%_Fem_GP_1    
 -200	       3	 -0.0521336	 -0.0521336	 99	0	 -3	0	0	0	0	  0	0	0	#_Mat_slope_Fem_GP_1 
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
3 #_Spawner-Recruitment; 2=Ricker (2 parms); 3=std_B-H(2); 4=SCAA(2);5=Hockey(3); 6=B-H_flattop(2); 7=Survival(3);8=Shepard(3);9=Ricker_Power(3);10=B-H_a,b(4)
0 # 0/1 to use steepness in initial equ recruitment calculation
0 # future feature: 0/1 to make realized sigmaR a function of SR curvature
#_LO	HI	INIT	PRIOR	PR_SD	PR_type	PHASE	env-var	use_dev	dev_mnyr	dev_mxyr	dev_PH	Block	Blk_Fxn # parm_name
   2	  20	 15	 15	  99	0	 1	0	0	0	0	0	0	0	#_SR_LN(R0)  
 0.2	0.99	0.8	0.8	1000	6	-2	0	0	0	0	0	0	0	#_SR_BH_steep
0.05	 1.9	0.6	0.6	1000	6	-4	0	0	0	0	0	0	0	#_SR_sigmaR  
  -4	   4	  0	  0	  99	0	-1	0	0	0	0	0	0	0	#_SR_regime  
   0	   0	  0	  0	  99	0	-1	0	0	0	0	0	0	0	#_SR_autocorr
#_no timevary SR parameters
1 #do_recdev:  0=none; 1=devvector (R=F(SSB)+dev); 2=deviations (R=F(SSB)+dev); 3=deviations (R=R0*dev; dev2=R-f(SSB)); 4=like 3 with sum(dev2) adding penalty
81 # first year of main recr_devs; early devs can preceed this era
267 # last year of main recr_devs; forecast devs start in following year
1 #_recdev phase
1 # (0/1) to read 13 advanced options
-40 #_recdev_early_start (0=none; neg value makes relative to recdev_start)
2 #_recdev_early_phase
0 #_forecast_recruitment phase (incl. late recr) (0 value resets to maxphase+1)
1 #_lambda for Fcast_recr_like occurring before endyr+1
41 #_last_yr_nobias_adj_in_MPD; begin of ramp
252 #_first_yr_fullbias_adj_in_MPD; begin of plateau
267 #_last_yr_fullbias_adj_in_MPD
267 #_end_yr_for_ramp_in_MPD (can be in forecast to shape ramp, but SS sets bias_adj to 0.0 for fcast yrs)
0.6 #_max_bias_adj_in_MPD (-1 to override ramp and set biasadj=1.0 for all estimated recdevs)
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
3 # N iterations for tuning F in hybrid method (recommend 3 to 7)
#
#_initial_F_parms; count = 0
#
#_Q_setup for fleets with cpue or survey data
#_fleet	link	link_info	extra_se	biasadj	float  #  fleetname
   15	1	0	0	0	1	#_S01_INDEX 
-9999	0	0	0	0	0	#_terminator
#_Q_parms(if_any);Qunits_are_ln(q)
#_LO	HI	INIT	PRIOR	PR_SD	PR_type	PHASE	env-var	use_dev	dev_mnyr	dev_mxyr	dev_PH	Block	Blk_Fxn  #  parm_name
-25	25	-10.265	0	1	0	-1	0	0	0	0	0	0	0	#_LnQ_base_S01_INDEX(15)
#_no timevary Q parameters
#
#_size_selex_patterns
#_Pattern	Discard	Male	Special
24	0	0	0	#_1 F01_LL.NORTH    
24	0	0	0	#_2 F02_LL.US       
24	0	0	0	#_3 F03_LL.OFFSHORE 
24	0	0	0	#_4 F04_LL.EQUAT    
24	0	0	0	#_5 F05_LL.WEST     
24	0	0	0	#_6 F06_LL.SOUTH    
24	0	0	0	#_7 F07_LL.AU       
24	0	0	0	#_8 F08_PS.ASS      
24	0	0	0	#_9 F09_PS.UNA      
24	0	0	0	#_10 F10_DOM.MISC   
 1	0	0	0	#_11 F11_DOM.HL     
24	0	0	0	#_12 F12_PS.JP.NORTH
24	0	0	0	#_13 F13_PL.JP.NORTH
24	0	0	0	#_14 F14_PL.EQUAT   
 1	0	0	0	#_15 S01_INDEX      
#
#_age_selex_patterns
#_Pattern	Discard	Male	Special
11	0	0	0	#_1 F01_LL.NORTH    
11	0	0	0	#_2 F02_LL.US       
11	0	0	0	#_3 F03_LL.OFFSHORE 
11	0	0	0	#_4 F04_LL.EQUAT    
11	0	0	0	#_5 F05_LL.WEST     
11	0	0	0	#_6 F06_LL.SOUTH    
11	0	0	0	#_7 F07_LL.AU       
11	0	0	0	#_8 F08_PS.ASS      
11	0	0	0	#_9 F09_PS.UNA      
11	0	0	0	#_10 F10_DOM.MISC   
11	0	0	0	#_11 F11_DOM.HL     
11	0	0	0	#_12 F12_PS.JP.NORTH
11	0	0	0	#_13 F13_PL.JP.NORTH
11	0	0	0	#_14 F14_PL.EQUAT   
11	0	0	0	#_15 S01_INDEX      
#
#_SizeSelex
#_LO	HI	INIT	PRIOR	PR_SD	PR_type	PHASE	env-var	use_dev	dev_mnyr	dev_mxyr	dev_PH	Block	Blk_Fxn  #  parm_name
10.1	200	  100.539	100	99	0	 2	0	0	0	0	0	0	0	#_SizeSel_P_1_F01_LL.NORTH(1)    
  -7	  7	       -5	 -5	99	0	-3	0	0	0	0	0	0	0	#_SizeSel_P_2_F01_LL.NORTH(1)    
  -7	  7	  6.21843	  0	99	0	-4	0	0	0	0	0	0	0	#_SizeSel_P_3_F01_LL.NORTH(1)    
  -7	  7	 -3.30342	  0	99	0	 5	0	0	0	0	0	0	0	#_SizeSel_P_4_F01_LL.NORTH(1)    
-999	  9	       -9	 -9	99	0	-4	0	0	0	0	0	0	0	#_SizeSel_P_5_F01_LL.NORTH(1)    
  -9	  9	   0.9987	  9	99	0	 6	0	0	0	0	0	0	0	#_SizeSel_P_6_F01_LL.NORTH(1)    
10.1	200	  105.512	100	99	0	 2	0	0	0	0	0	0	0	#_SizeSel_P_1_F02_LL.US(2)       
  -7	  7	       -5	 -5	99	0	-3	0	0	0	0	0	0	0	#_SizeSel_P_2_F02_LL.US(2)       
  -7	  7	  6.22456	  0	99	0	-4	0	0	0	0	0	0	0	#_SizeSel_P_3_F02_LL.US(2)       
  -7	  7	 -4.62114	  0	99	0	 5	0	0	0	0	0	0	0	#_SizeSel_P_4_F02_LL.US(2)       
-999	  9	       -9	 -9	99	0	-4	0	0	0	0	0	0	0	#_SizeSel_P_5_F02_LL.US(2)       
  -9	  9	  1.45676	  9	99	0	 6	0	0	0	0	0	0	0	#_SizeSel_P_6_F02_LL.US(2)       
10.1	200	  137.652	100	99	0	 2	0	0	0	0	0	0	0	#_SizeSel_P_1_F03_LL.OFFSHORE(3) 
  -7	  7	       -5	 -5	99	0	-3	0	0	0	0	0	0	0	#_SizeSel_P_2_F03_LL.OFFSHORE(3) 
  -7	  7	  6.51135	  0	99	0	-4	0	0	0	0	0	0	0	#_SizeSel_P_3_F03_LL.OFFSHORE(3) 
  -7	  7	  6.84653	  0	99	0	-5	0	0	0	0	0	0	0	#_SizeSel_P_4_F03_LL.OFFSHORE(3) 
-999	  9	       -9	 -9	99	0	-4	0	0	0	0	0	0	0	#_SizeSel_P_5_F03_LL.OFFSHORE(3) 
-999	  9	     -495	  9	99	0	 6	0	0	0	0	0	0	0	#_SizeSel_P_6_F03_LL.OFFSHORE(3) 
10.1	200	  127.789	100	99	0	 2	0	0	0	0	0	0	0	#_SizeSel_P_1_F04_LL.EQUAT(4)    
  -7	  7	       -5	 -5	99	0	-3	0	0	0	0	0	0	0	#_SizeSel_P_2_F04_LL.EQUAT(4)    
  -7	  7	  6.87501	  0	99	0	-4	0	0	0	0	0	0	0	#_SizeSel_P_3_F04_LL.EQUAT(4)    
  -7	  7	  6.87243	  0	99	0	-5	0	0	0	0	0	0	0	#_SizeSel_P_4_F04_LL.EQUAT(4)    
-999	  9	       -9	 -9	99	0	-4	0	0	0	0	0	0	0	#_SizeSel_P_5_F04_LL.EQUAT(4)    
-999	  9	     -495	  9	99	0	 6	0	0	0	0	0	0	0	#_SizeSel_P_6_F04_LL.EQUAT(4)    
10.1	200	  63.4502	100	99	0	 2	0	0	0	0	0	0	0	#_SizeSel_P_1_F05_LL.WEST(5)     
  -7	  7	       -5	 -5	99	0	-3	0	0	0	0	0	0	0	#_SizeSel_P_2_F05_LL.WEST(5)     
  -7	  7	 -6.66193	  0	99	0	-4	0	0	0	0	0	0	0	#_SizeSel_P_3_F05_LL.WEST(5)     
  -7	  7	-0.702849	  0	99	0	 5	0	0	0	0	0	0	0	#_SizeSel_P_4_F05_LL.WEST(5)     
-999	  9	       -9	 -9	99	0	-4	0	0	0	0	0	0	0	#_SizeSel_P_5_F05_LL.WEST(5)     
  -9	  9	   5.8491	  9	99	0	 6	0	0	0	0	0	0	0	#_SizeSel_P_6_F05_LL.WEST(5)     
10.1	200	  116.594	100	99	0	 2	0	0	0	0	0	0	0	#_SizeSel_P_1_F06_LL.SOUTH(6)    
  -7	  7	       -5	 -5	99	0	-3	0	0	0	0	0	0	0	#_SizeSel_P_2_F06_LL.SOUTH(6)    
  -7	  7	  6.70807	  0	99	0	-4	0	0	0	0	0	0	0	#_SizeSel_P_3_F06_LL.SOUTH(6)    
  -7	  7	  6.55203	  0	99	0	-5	0	0	0	0	0	0	0	#_SizeSel_P_4_F06_LL.SOUTH(6)    
-999	  9	       -9	 -9	99	0	-4	0	0	0	0	0	0	0	#_SizeSel_P_5_F06_LL.SOUTH(6)    
  -9	  9	 0.889611	  9	99	0	 6	0	0	0	0	0	0	0	#_SizeSel_P_6_F06_LL.SOUTH(6)    
10.1	200	  110.759	100	99	0	 2	0	0	0	0	0	0	0	#_SizeSel_P_1_F07_LL.AU(7)       
  -7	  7	       -5	 -5	99	0	-3	0	0	0	0	0	0	0	#_SizeSel_P_2_F07_LL.AU(7)       
  -7	  7	  6.24774	  0	99	0	-4	0	0	0	0	0	0	0	#_SizeSel_P_3_F07_LL.AU(7)       
  -7	  7	  6.93819	  0	99	0	-5	0	0	0	0	0	0	0	#_SizeSel_P_4_F07_LL.AU(7)       
-999	  9	       -9	 -9	99	0	-4	0	0	0	0	0	0	0	#_SizeSel_P_5_F07_LL.AU(7)       
-999	  9	 -1.29248	  9	99	0	 6	0	0	0	0	0	0	0	#_SizeSel_P_6_F07_LL.AU(7)       
10.1	200	  50.7393	 40	99	0	 2	0	0	0	0	0	0	0	#_SizeSel_P_1_F08_PS.ASS(8)      
  -7	  7	       -5	 -5	99	0	-3	0	0	0	0	0	0	0	#_SizeSel_P_2_F08_PS.ASS(8)      
  -7	  7	  4.05358	  0	99	0	 4	0	0	0	0	0	0	0	#_SizeSel_P_3_F08_PS.ASS(8)      
  -7	  7	   4.4569	  0	99	0	 5	0	0	0	0	0	0	0	#_SizeSel_P_4_F08_PS.ASS(8)      
-999	  9	       -9	 -9	99	0	-4	0	0	0	0	0	0	0	#_SizeSel_P_5_F08_PS.ASS(8)      
-999	  9	       -9	 -9	99	0	-6	0	0	0	0	0	0	0	#_SizeSel_P_6_F08_PS.ASS(8)      
10.1	200	  61.9319	 40	99	0	 2	0	0	0	0	0	0	0	#_SizeSel_P_1_F09_PS.UNA(9)      
  -7	  7	       -5	 -5	99	0	-3	0	0	0	0	0	0	0	#_SizeSel_P_2_F09_PS.UNA(9)      
  -7	  7	  5.63087	  0	99	0	 4	0	0	0	0	0	0	0	#_SizeSel_P_3_F09_PS.UNA(9)      
  -7	  7	  6.99613	  0	99	0	-5	0	0	0	0	0	0	0	#_SizeSel_P_4_F09_PS.UNA(9)      
-999	  9	       -9	 -9	99	0	-4	0	0	0	0	0	0	0	#_SizeSel_P_5_F09_PS.UNA(9)      
-999	  9	       -9	 -9	99	0	-6	0	0	0	0	0	0	0	#_SizeSel_P_6_F09_PS.UNA(9)      
10.1	200	  31.4346	 40	99	0	 2	0	0	0	0	0	0	0	#_SizeSel_P_1_F10_DOM.MISC(10)   
  -7	  7	       -5	 -5	99	0	-3	0	0	0	0	0	0	0	#_SizeSel_P_2_F10_DOM.MISC(10)   
  -7	  7	  4.27727	  0	99	0	 4	0	0	0	0	0	0	0	#_SizeSel_P_3_F10_DOM.MISC(10)   
  -7	  7	  5.67679	  0	99	0	 5	0	0	0	0	0	0	0	#_SizeSel_P_4_F10_DOM.MISC(10)   
-999	  9	     -495	  0	99	0	 6	0	0	0	0	0	0	0	#_SizeSel_P_5_F10_DOM.MISC(10)   
-999	  9	       -9	 -9	99	0	-6	0	0	0	0	0	0	0	#_SizeSel_P_6_F10_DOM.MISC(10)   
   5	200	  140.008	100	99	0	 2	0	0	0	0	0	0	0	#_SizeSel_P_1_F11_DOM.HL(11)     
   0	500	  56.9346	 10	99	0	 3	0	0	0	0	0	0	0	#_SizeSel_P_2_F11_DOM.HL(11)     
10.1	200	  50.6613	 40	99	0	 2	0	0	0	0	0	0	0	#_SizeSel_P_1_F12_PS.JP.NORTH(12)
  -7	  7	       -5	 -5	99	0	-3	0	0	0	0	0	0	0	#_SizeSel_P_2_F12_PS.JP.NORTH(12)
  -7	  7	  2.20661	  0	99	0	 4	0	0	0	0	0	0	0	#_SizeSel_P_3_F12_PS.JP.NORTH(12)
  -7	  7	  6.97642	  0	99	0	-5	0	0	0	0	0	0	0	#_SizeSel_P_4_F12_PS.JP.NORTH(12)
-999	  9	       -9	 -9	99	0	-4	0	0	0	0	0	0	0	#_SizeSel_P_5_F12_PS.JP.NORTH(12)
-999	  9	       -9	 -9	99	0	-6	0	0	0	0	0	0	0	#_SizeSel_P_6_F12_PS.JP.NORTH(12)
10.1	200	  14.5691	 40	99	0	 2	0	0	0	0	0	0	0	#_SizeSel_P_1_F13_PL.JP.NORTH(13)
  -7	  7	       -5	 -5	99	0	-3	0	0	0	0	0	0	0	#_SizeSel_P_2_F13_PL.JP.NORTH(13)
  -7	  7	  6.61726	  0	99	0	-4	0	0	0	0	0	0	0	#_SizeSel_P_3_F13_PL.JP.NORTH(13)
  -7	  7	  3.99651	  0	99	0	 5	0	0	0	0	0	0	0	#_SizeSel_P_4_F13_PL.JP.NORTH(13)
-999	  9	       -9	 -9	99	0	-4	0	0	0	0	0	0	0	#_SizeSel_P_5_F13_PL.JP.NORTH(13)
-999	  9	       -9	 -9	99	0	-6	0	0	0	0	0	0	0	#_SizeSel_P_6_F13_PL.JP.NORTH(13)
10.1	200	  34.7299	 40	99	0	 2	0	0	0	0	0	0	0	#_SizeSel_P_1_F14_PL.EQUAT(14)   
  -7	  7	       -5	 -5	99	0	-3	0	0	0	0	0	0	0	#_SizeSel_P_2_F14_PL.EQUAT(14)   
  -7	  7	  4.71981	  0	99	0	 4	0	0	0	0	0	0	0	#_SizeSel_P_3_F14_PL.EQUAT(14)   
  -7	  7	  4.51849	  0	99	0	 5	0	0	0	0	0	0	0	#_SizeSel_P_4_F14_PL.EQUAT(14)   
-999	  9	       -9	 -9	99	0	-4	0	0	0	0	0	0	0	#_SizeSel_P_5_F14_PL.EQUAT(14)   
-999	  9	       -9	 -9	99	0	-6	0	0	0	0	0	0	0	#_SizeSel_P_6_F14_PL.EQUAT(14)   
   5	200	  89.7616	100	99	0	 2	0	0	0	0	0	0	0	#_SizeSel_P_1_S01_INDEX(15)      
   0	500	  18.2857	 10	99	0	 3	0	0	0	0	0	0	0	#_SizeSel_P_2_S01_INDEX(15)      
#_AgeSelex
 0	 10	 0	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_1_F01_LL.NORTH(1)    
10	100	40	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_2_F01_LL.NORTH(1)    
 0	 10	 0	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_1_F02_LL.US(2)       
10	100	40	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_2_F02_LL.US(2)       
 0	 10	 0	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_1_F03_LL.OFFSHORE(3) 
10	100	40	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_2_F03_LL.OFFSHORE(3) 
 0	 10	 0	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_1_F04_LL.EQUAT(4)    
10	100	40	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_2_F04_LL.EQUAT(4)    
 0	 10	 0	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_1_F05_LL.WEST(5)     
10	100	40	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_2_F05_LL.WEST(5)     
 0	 10	 0	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_1_F06_LL.SOUTH(6)    
10	100	40	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_2_F06_LL.SOUTH(6)    
 0	 10	 0	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_1_F07_LL.AU(7)       
10	100	40	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_2_F07_LL.AU(7)       
 0	 10	 0	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_1_F08_PS.ASS(8)      
10	100	40	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_2_F08_PS.ASS(8)      
 0	 10	 0	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_1_F09_PS.UNA(9)      
10	100	40	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_2_F09_PS.UNA(9)      
 0	 10	 0	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_1_F10_DOM.MISC(10)   
10	100	40	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_2_F10_DOM.MISC(10)   
 0	 10	 0	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_1_F11_DOM.HL(11)     
10	100	40	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_2_F11_DOM.HL(11)     
 0	 10	 0	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_1_F12_PS.JP.NORTH(12)
10	100	40	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_2_F12_PS.JP.NORTH(12)
 0	 10	 0	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_1_F13_PL.JP.NORTH(13)
10	100	40	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_2_F13_PL.JP.NORTH(13)
 0	 10	 0	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_1_F14_PL.EQUAT(14)   
10	100	40	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_2_F14_PL.EQUAT(14)   
 0	 10	 0	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_1_S01_INDEX(15)      
10	100	40	0	99	0	-99	0	0	0	0	0	0	0	#_AgeSel_P_2_S01_INDEX(15)      
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
    1	 1	      0	#_Variance_adjustment_list1  
    2	 1	      0	#_Variance_adjustment_list2  
    3	 1	      0	#_Variance_adjustment_list3  
    4	 1	2.5e-05	#_Variance_adjustment_list4  
    5	 1	      0	#_Variance_adjustment_list5  
    6	 1	      0	#_Variance_adjustment_list6  
    7	 1	2.5e-05	#_Variance_adjustment_list7  
    1	 2	      0	#_Variance_adjustment_list8  
    2	 2	      0	#_Variance_adjustment_list9  
    3	 2	      0	#_Variance_adjustment_list10 
    4	 2	  5e-05	#_Variance_adjustment_list11 
    5	 2	      0	#_Variance_adjustment_list12 
    6	 2	      0	#_Variance_adjustment_list13 
    7	 2	  5e-05	#_Variance_adjustment_list14 
    1	 3	      0	#_Variance_adjustment_list15 
    2	 3	      0	#_Variance_adjustment_list16 
    3	 3	      0	#_Variance_adjustment_list17 
    4	 3	  5e-05	#_Variance_adjustment_list18 
    5	 3	      0	#_Variance_adjustment_list19 
    6	 3	      0	#_Variance_adjustment_list20 
    7	 3	  5e-05	#_Variance_adjustment_list21 
    1	 4	      0	#_Variance_adjustment_list22 
    2	 4	      0	#_Variance_adjustment_list23 
    3	 4	      0	#_Variance_adjustment_list24 
    4	 4	2.5e-05	#_Variance_adjustment_list25 
    5	 4	      0	#_Variance_adjustment_list26 
    6	 4	      0	#_Variance_adjustment_list27 
    7	 4	2.5e-05	#_Variance_adjustment_list28 
    1	 5	      0	#_Variance_adjustment_list29 
    2	 5	      0	#_Variance_adjustment_list30 
    3	 5	      0	#_Variance_adjustment_list31 
    4	 5	2.5e-05	#_Variance_adjustment_list32 
    5	 5	      0	#_Variance_adjustment_list33 
    6	 5	      0	#_Variance_adjustment_list34 
    7	 5	2.5e-05	#_Variance_adjustment_list35 
    1	 6	      0	#_Variance_adjustment_list36 
    2	 6	      0	#_Variance_adjustment_list37 
    3	 6	      0	#_Variance_adjustment_list38 
    4	 6	2.5e-05	#_Variance_adjustment_list39 
    5	 6	      0	#_Variance_adjustment_list40 
    6	 6	      0	#_Variance_adjustment_list41 
    7	 6	2.5e-05	#_Variance_adjustment_list42 
    1	 7	      0	#_Variance_adjustment_list43 
    2	 7	      0	#_Variance_adjustment_list44 
    3	 7	      0	#_Variance_adjustment_list45 
    4	 7	  5e-05	#_Variance_adjustment_list46 
    5	 7	      0	#_Variance_adjustment_list47 
    6	 7	      0	#_Variance_adjustment_list48 
    7	 7	  5e-05	#_Variance_adjustment_list49 
    1	 8	      0	#_Variance_adjustment_list50 
    2	 8	      0	#_Variance_adjustment_list51 
    3	 8	      0	#_Variance_adjustment_list52 
    4	 8	  5e-05	#_Variance_adjustment_list53 
    5	 8	      0	#_Variance_adjustment_list54 
    6	 8	      0	#_Variance_adjustment_list55 
    7	 8	  5e-05	#_Variance_adjustment_list56 
    1	 9	      0	#_Variance_adjustment_list57 
    2	 9	      0	#_Variance_adjustment_list58 
    3	 9	      0	#_Variance_adjustment_list59 
    4	 9	  5e-05	#_Variance_adjustment_list60 
    5	 9	      0	#_Variance_adjustment_list61 
    6	 9	      0	#_Variance_adjustment_list62 
    7	 9	  5e-05	#_Variance_adjustment_list63 
    1	10	      0	#_Variance_adjustment_list64 
    2	10	      0	#_Variance_adjustment_list65 
    3	10	      0	#_Variance_adjustment_list66 
    4	10	  5e-05	#_Variance_adjustment_list67 
    5	10	      0	#_Variance_adjustment_list68 
    6	10	      0	#_Variance_adjustment_list69 
    7	10	  5e-05	#_Variance_adjustment_list70 
    1	11	      0	#_Variance_adjustment_list71 
    2	11	      0	#_Variance_adjustment_list72 
    3	11	      0	#_Variance_adjustment_list73 
    4	11	  5e-05	#_Variance_adjustment_list74 
    5	11	      0	#_Variance_adjustment_list75 
    6	11	      0	#_Variance_adjustment_list76 
    7	11	  5e-05	#_Variance_adjustment_list77 
    1	12	      0	#_Variance_adjustment_list78 
    2	12	      0	#_Variance_adjustment_list79 
    3	12	      0	#_Variance_adjustment_list80 
    4	12	  5e-05	#_Variance_adjustment_list81 
    5	12	      0	#_Variance_adjustment_list82 
    6	12	      0	#_Variance_adjustment_list83 
    7	12	  5e-05	#_Variance_adjustment_list84 
    1	13	      0	#_Variance_adjustment_list85 
    2	13	      0	#_Variance_adjustment_list86 
    3	13	      0	#_Variance_adjustment_list87 
    4	13	  5e-05	#_Variance_adjustment_list88 
    5	13	      0	#_Variance_adjustment_list89 
    6	13	      0	#_Variance_adjustment_list90 
    7	13	  5e-05	#_Variance_adjustment_list91 
    1	14	      0	#_Variance_adjustment_list92 
    2	14	      0	#_Variance_adjustment_list93 
    3	14	      0	#_Variance_adjustment_list94 
    4	14	  5e-05	#_Variance_adjustment_list95 
    5	14	      0	#_Variance_adjustment_list96 
    6	14	      0	#_Variance_adjustment_list97 
    7	14	  5e-05	#_Variance_adjustment_list98 
    1	15	      0	#_Variance_adjustment_list99 
    2	15	      0	#_Variance_adjustment_list100
    3	15	      0	#_Variance_adjustment_list101
    4	15	2.5e-05	#_Variance_adjustment_list102
    5	15	      0	#_Variance_adjustment_list103
    6	15	      0	#_Variance_adjustment_list104
    7	15	2.5e-05	#_Variance_adjustment_list105
-9999	 0	      0	#_terminator                 
#
1 #_maxlambdaphase
1 #_sd_offset; must be 1 if any growthCV, sigmaR, or survey extraSD is an estimated parameter
# read 38 changes to default Lambdas (default value is 1.0)
#_like_comp	fleet	phase	value	sizefreq_method
    1	 1	1	0	0	#_1         
    1	 2	1	0	0	#_2         
    1	 3	1	0	0	#_3         
    1	 4	1	0	0	#_4         
    1	 5	1	0	0	#_5         
    1	 6	1	0	0	#_6         
    1	 7	1	0	0	#_7         
    1	 8	1	0	0	#_8         
    1	 9	1	0	0	#_9         
    1	10	1	0	0	#_10        
    1	11	1	0	0	#_11        
    1	12	1	0	0	#_12        
    1	13	1	0	0	#_13        
    1	14	1	0	0	#_14        
    1	15	1	1	0	#_15        
    4	 1	1	0	0	#_16        
    4	 2	1	0	0	#_17        
    4	 3	1	0	0	#_18        
    4	 4	1	0	0	#_19        
    4	 5	1	0	0	#_20        
    4	 6	1	0	0	#_21        
    4	 7	1	0	0	#_22        
    4	 8	1	1	0	#_23        
    4	 9	1	1	0	#_24        
    4	10	1	1	0	#_25        
    4	11	1	1	0	#_26        
    4	12	1	1	0	#_27        
    4	13	1	1	0	#_28        
    4	14	1	1	0	#_29        
    4	15	1	0	0	#_30        
    6	 1	1	1	1	#_31        
    6	 2	1	1	1	#_32        
    6	 3	1	1	1	#_33        
    6	 4	1	1	1	#_34        
    6	 5	1	1	1	#_35        
    6	 6	1	1	1	#_36        
    6	 7	1	1	1	#_37        
    6	15	1	1	1	#_38        
-9999	 0	0	0	0	#_terminator
#
2 # 0/1 read specs for more stddev reporting
0 0 0 0 # selex_fleet, 1=len/2=age/3=both, year, N selex bins
0 0       # 0 or Growth pattern, N growth ages
0 0 0    # 0 or NatAge_area(-1 for sum), NatAge_yr, N Natages
0 0 1 0     # Mortality, Dyn B0 (>3.30.16), SmryBio (>3.30.16) 
#
999
