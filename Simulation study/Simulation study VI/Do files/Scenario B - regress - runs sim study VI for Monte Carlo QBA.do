/*************************************************************************************************************************************************
                                                       RUNS SIMULATION STUDY I FOR MONTE CARLO QBA 
											SCENARIO B: Y CONTINUOUS, X CONTINUOUS, MULTIPLE U=(U1,U2) CONTINUOUS

ARGUMENTS:
		samplesize       - SAMPLE SIZE OF EACH SIMULATED DATASET
		numsimdatasets   - NUMBER OF SIMULATED DATASETS
		numMCreps        - NUMBER OF MONTE CARLO REPLICATIONS
		priorlevel		 - LEVEL OF PRECISION OF BIAS PARAMETERS' PRIORS
		uniform			 - PRIOR FOR eta IS UNIFORM DISTRIBUTION FOR (uniform=1); OTHERWISE GAMMA DISTRIBUTION ON 1/etasq
		
RETURNS: 
	DATASET CONTAINING RESULTS STORED BY POSTFILE; STORED UNDER FILENAME OF MACRO resultsfile
	STORED IN FOLDER "Results"
*************************************************************************************************************************************************/
args samplesize numsimdatasets numMCreps priorlevel uniform

noisily di "Analysing `numsimdatasets' simulated datasets; each of `samplesize' observations"
noisily di "Applying Monte Carlo QBA to generate a frequency distribution of `numMCreps' estimates"

* DROP ALL SCALARS
scalar drop _all

/****************************************************
   MACROS AND TEMPORARY NAMES, FILES, AND VARIABLES
***************************************************/
tempfile mcqba_replications
tempname memhold alpha2x_hyperp eta2_hyperp alpha1x_hyperp eta1_hyperp betau1_hyperp betau2_hyperp   

	* MACROLISTS
local exposurelist "xcon"
local confounderlist `"ccon cbin _Icnom_1 _Icnom_2"'
local parameterlist `"alpha2star_xcon eta2star alpha1star_xcon eta1star betastar_ucon1 betastar_ucon2 "'

* HYPERPARAMETERS OF THE PRIOR DISTRIBUTIONS FOR THE BIAS PARAMETERS
	* UNIFORM DISTRIBUTION FOR eta1 AND eta2, OR GAMMA DISTRIBUTION FOR 1/eta1sq AND 1/eta2sq
	* NORMAL DISTRIBUTIONS FOR COEFFCIENTS 
		* FOR PRIOR LEVELS 2 AND 3: MEANS SET TO TRUE VALUES; SDs SET TO DIFFERENT VALUES
		* FOR PRIOR LEVELS 4 AND 5: MEANS SET TO VALUES SHIFTED BY 30%; SDs SET SET TO DIFFERENT VALUES 
if `priorlevel'==2 {
	noisily di "Applying Monte Carlo QBA using accurate, informative prior distributions"
	local priortxt "medium"			
	* alpha1_xcon
	local mean_alpha1_xcon 1.260
	local sd_alpha1_xcon 0.1929
	* alpha2_xcon
	local mean_alpha2_xcon 1.281
	local sd_alpha2_xcon 0.1961
	* beta_ucon1
	local mean_beta_ucon1 0.1299
	local sd_beta_ucon1 0.01988 
	* beta_ucon2
	local mean_beta_ucon2 0.06747
	local sd_beta_ucon2 0.01033
	* eta1
	local luni_eta1 6.427
	local uuni_eta1 12.36
	local shape_inveta1sq 13.19
	local scale_inveta1sq = 1/1164
	* eta2
	local luni_eta2 7.047
	local uuni_eta2 13.55
	local shape_inveta2sq 13.19
	local scale_inveta2sq = 1/1399
}
else if `priorlevel'==3 {
	noisily di "Applying Monte Carlo QBA using accurate, very informative prior distributions"
	local priortxt "most"
	local mean_alpha2_xcon 1.281
	local mean_alpha1_xcon 1.260
	local mean_beta_ucon1 0.1299
	local mean_beta_ucon2 0.06747

	local sd_alpha2_xcon 0.06536
	local luni_eta2 9.216
	local uuni_eta2 11.38
	local shape_inveta2sq .
	local scale_inveta2sq .
	local sd_alpha1_xcon 0.06429
	local luni_eta1 8.405
	local uuni_eta1 10.38
	local shape_inveta1sq .
	local scale_inveta1sq .
	local sd_beta_ucon1 0.006628
	local sd_beta_ucon2 0.003442
}  
else if `priorlevel'==4 {
	noisily di "Applying Monte Carlo QBA using inaccurate (30% below truth), informative prior distributions"
	local priortxt "inaccmedium"
	* alpha1_xcon
	local mean_alpha1_xcon = 0.7*1.260
	local sd_alpha1_xcon 0.1929
	* alpha2_xcon
	local mean_alpha2_xcon = 0.7*1.281
	local sd_alpha2_xcon 0.1961
	* beta_ucon1
	local mean_beta_ucon1 = 0.7*0.1299
	local sd_beta_ucon1 0.01988 
	* beta_ucon2
	local mean_beta_ucon2 = 0.7*0.06747
	local sd_beta_ucon2 0.01033
	* eta1
	local acc_luni_eta1 6.427
	local acc_uuni_eta1 12.36
	local mean_eta1 = 0.7*9.394
	local luni_eta1 = `mean_eta1' - 0.5*(`acc_uuni_eta1'-`acc_luni_eta1')
	local uuni_eta1 = `mean_eta1' + 0.5*(`acc_uuni_eta1'-`acc_luni_eta1')
	local shape_inveta1sq .
	local scale_inveta1sq .
	* eta2
	local acc_luni_eta2 7.047
	local acc_uuni_eta2 13.55
	local mean_eta2 = 0.7*10.30
	local luni_eta2 = `mean_eta2' - 0.5*(`acc_uuni_eta2'-`acc_luni_eta2')
	local uuni_eta2 = `mean_eta2' + 0.5*(`acc_uuni_eta2'-`acc_luni_eta2')
	local shape_inveta2sq .
	local scale_inveta2sq .
}
else if `priorlevel'==5 {
	noisily di "Applying Monte Carlo QBA using inaccurate (30% below truth), less informative (sd double size of informative) priors"
	local priortxt "inaccweak"
	* alpha1_xcon
	local mean_alpha1_xcon = 0.7*1.260
	local sd_alpha1_xcon = 2*0.1929
	* alpha2_xcon
	local mean_alpha2_xcon = 0.7*1.281
	local sd_alpha2_xcon = 2*0.1961
	* beta_ucon1
	local mean_beta_ucon1 = 0.7*0.1299
	local sd_beta_ucon1 = 2*0.01988 
	* beta_ucon2
	local mean_beta_ucon2 = 0.7*0.06747
	local sd_beta_ucon2 = 2*0.01033
	* eta1
	local acc_luni_eta1 6.427
	local acc_uuni_eta1 12.36
	local mean_eta1 = 0.7*9.394
	local luni_eta1 = `mean_eta1' - (`acc_uuni_eta1'-`acc_luni_eta1')
	local uuni_eta1 = `mean_eta1' + (`acc_uuni_eta1'-`acc_luni_eta1')
	local shape_inveta1sq .
	local scale_inveta1sq .
	* eta2
	local acc_luni_eta2 7.047
	local acc_uuni_eta2 13.55
	local mean_eta2 = 0.7*10.30
	local luni_eta2 = `mean_eta2' - (`acc_uuni_eta2'-`acc_luni_eta2')
	local uuni_eta2 = `mean_eta2' + (`acc_uuni_eta2'-`acc_luni_eta2')
	local shape_inveta2sq .
	local scale_inveta2sq .
}
matrix `alpha2x_hyperp' = (`mean_alpha2_xcon', `sd_alpha2_xcon')
matrix `eta2_hyperp' = (`luni_eta2',`uuni_eta2',`shape_inveta2sq',`scale_inveta2sq')
matrix `alpha1x_hyperp' = (`mean_alpha1_xcon', `sd_alpha1_xcon')
matrix `eta1_hyperp' = (`luni_eta1',`uuni_eta1',`shape_inveta1sq',`scale_inveta1sq')
matrix `betau1_hyperp' = (`mean_beta_ucon1', `sd_beta_ucon1')
matrix `betau2_hyperp' = (`mean_beta_ucon2', `sd_beta_ucon2')

/*********************************
  MACROLISTS FOR FRAME STATEMENTS
*********************************/
* BIAS ANALYSIS
local vars_bapd1 `"b_bapd_sys_betaxcon sd_bapd_sys_betaxcon lci_bapd_sys_betaxcon uci_bapd_sys_betaxcon b_bapd_all_betaxcon sd_bapd_all_betaxcon lci_bapd_all_betaxcon uci_bapd_all_betaxcon bapd_min_alpha2star_xcon bapd_p025_alpha2star_xcon bapd_p25_alpha2star_xcon bapd_p50_alpha2star_xcon bapd_p75_alpha2star_xcon bapd_p975_alpha2star_xcon bapd_max_alpha2star_xcon bapd_min_eta2star bapd_p025_eta2star bapd_p25_eta2star bapd_p50_eta2star bapd_p75_eta2star bapd_p975_eta2star bapd_max_eta2star bapd_min_alpha1star_xcon bapd_p025_alpha1star_xcon bapd_p25_alpha1star_xcon bapd_p50_alpha1star_xcon bapd_p75_alpha1star_xcon bapd_p975_alpha1star_xcon bapd_max_alpha1star_xcon"' 

local vars_bapd2 `"bapd_min_eta1star bapd_p025_eta1star bapd_p25_eta1star bapd_p50_eta1star bapd_p75_eta1star bapd_p975_eta1star bapd_max_eta1star bapd_min_betastar_ucon1 bapd_p025_betastar_ucon1 bapd_p25_betastar_ucon1 bapd_p50_betastar_ucon1 bapd_p75_betastar_ucon1 bapd_p975_betastar_ucon1 bapd_max_betastar_ucon1 bapd_min_betastar_ucon2 bapd_p025_betastar_ucon2 bapd_p25_betastar_ucon2 bapd_p50_betastar_ucon2 bapd_p75_betastar_ucon2 bapd_p975_betastar_ucon2 bapd_max_betastar_ucon2 bapd_timer bapd_num_analysisrc0 bapd_min_analysis_eN bapd_p25_analysis_eN bapd_p50_analysis_eN bapd_p75_analysis_eN bapd_max_analysis_eN"'

* LOCAL MACROLISTS FOR RESULTS TO BE POSTED TO THE FRAME
local list_roots `"bapd1 bapd2"' 	
local list_vars ""
local list_posts ""
foreach root of local list_roots {
	local post_`root' ""
	foreach item of local vars_`root' {
		local posting `"(`item')"'
		local post_`root': list post_`root' | posting
	}
	local list_vars: list list_vars | vars_`root'
	local list_posts: list list_posts | post_`root'
}
*noisily di "`list_vars'"
*noisily di "`list_posts'"

/*************************************
   CODE TO RUN THE SIMULATION STUDY
*************************************/
local resultsfile ///
`"Simulation study V\Results\ScenarioB_regress_`samplesize'obs_`priortxt'inf_uniform`uniform'_`numMCreps'reps_strengthUC1.dta"'

frame create `memhold' dataset strL(state) mean_alpha2_xcon sd_alpha2_xcon luni_eta2 uuni_eta2 shape_inveta2sq scale_inveta2sq ///
                                           mean_alpha1_xcon sd_alpha1_xcon luni_eta1 uuni_eta1 shape_inveta1sq scale_inveta1sq ///
										   mean_beta_ucon1 sd_beta_ucon1 mean_beta_ucon2 sd_beta_ucon2 ///
										   bapd_MCQBA_rc `list_vars' 
forvalues dataset=1(1)`numsimdatasets' {
	noisily di ///
    "YconXconUcon; regress; `priortxt' informative prior; Uniform=`uniform'; `numMCreps' MCQBA reps: Processing dataset `dataset' of sample size `samplesize'"	
	
quietly {   
	* RECORD THE rngstate FOR THIS LOOP
	local state = c(rngstate)
	
	* IMPORT DATASET
	import delimited "Simulation study I\Scenario B\Data\Dataset_`dataset'.csv", clear
	
	* GENERATE DUMMY VARIABLES OF cnom FOR ANALYSIS
	xi i.cnom
	
	/***********************************************
	 APPLY MONTE CARLO QBA TO SIMULATED DATASET 
    ***********************************************/
	timer clear 1
	timer on 1
	capture run "Fits MCQBA\Scenario B - fits MCQBA.do" "`exposurelist'" "`confounderlist'" ///
																		  `alpha2x_hyperp' `eta2_hyperp' ///
																		  `alpha1x_hyperp' `eta1_hyperp' ///
																		  `betau1_hyperp' `betau2_hyperp' ///
																		  `numMCreps' `mcqba_replications' `uniform'
	local bapd_MCQBA_rc = _rc
	timer off 1
	timer list 1
	scalar bapd_timer = r(t1)

	if `bapd_MCQBA_rc'==0 {
		
		* PROCESS MCQBA RESULTS
		use `mcqba_replications', clear
		
		* SUMMARY STATISTICS ON THE SUCCESSFUL BIAS PARAMETER DRAWS
		foreach parameter of local parameterlist {
			centile `parameter', centile(2.5)
			scalar bapd_p025_`parameter' = r(c_1)
			centile `parameter', centile(97.5)
			scalar bapd_p975_`parameter' = r(c_1)
			
			summarize `parameter', detail
			scalar bapd_min_`parameter' = r(min)
			scalar bapd_p25_`parameter' = r(p25)
			scalar bapd_p50_`parameter' = r(p50)
			scalar bapd_p75_`parameter' = r(p75)
			scalar bapd_max_`parameter' = r(max)
		}
		
		* CONVERGENCE STATISTICS OF THE ANALYSIS MODEL
			* NUMBER THAT CONVERGED
		count if analysis_rc == 0
		scalar bapd_num_analysisrc0 = r(N)
		
			* KEEP THE CONVERGED RESULTS
		keep if analysis_rc == 0
		
		summ analysis_eN, detail
		scalar bapd_min_analysis_eN = r(min)
		scalar bapd_p25_analysis_eN = r(p25)
		scalar bapd_p50_analysis_eN = r(p50)
		scalar bapd_p75_analysis_eN = r(p75)
		scalar bapd_max_analysis_eN = r(max)

			* RESULTS ACCOUNTING FOR UNCERTAINTY ABOUT SYSTEMATIC ERROR ONLY
		centile betahat_xcon, centile(50)
		scalar b_bapd_sys_betaxcon = r(c_1)
		centile betahat_xcon, centile(2.5)
		scalar lci_bapd_sys_betaxcon = r(c_1)
		centile betahat_xcon, centile(97.5)
		scalar uci_bapd_sys_betaxcon = r(c_1)
		summarize betahat_xcon
		scalar sd_bapd_sys_betaxcon	= r(sd)
		
			* RESULTS ACCOUNTING FOR UNCERTAINTY ABOUT SYSTEMATIC ERROR AND RANDOM SAMPLING
		centile betastar_xcon, centile(50)
		scalar b_bapd_all_betaxcon = r(c_1)
		centile betastar_xcon, centile(2.5)
		scalar lci_bapd_all_betaxcon = r(c_1)
		centile betastar_xcon, centile(97.5)
		scalar uci_bapd_all_betaxcon = r(c_1)
		summarize betastar_xcon 
		scalar sd_bapd_all_betaxcon = r(sd)	
	}
	else {
		foreach item of local vars_bapd {
			scalar `item' = .
		}			
	}	
	
	frame post `memhold' (`dataset') ("`state'") (`mean_alpha2_xcon') (`sd_alpha2_xcon') (`luni_eta2') (`uuni_eta2') ///
												 (`shape_inveta2sq') (`scale_inveta2sq') ///
	                                             (`mean_alpha1_xcon') (`sd_alpha1_xcon') (`luni_eta1') (`uuni_eta1') ///
												 (`shape_inveta1sq') (`scale_inveta1sq') ///
										         (`mean_beta_ucon1') (`sd_beta_ucon1') (`mean_beta_ucon2') (`sd_beta_ucon2') ///
	                                             (`bapd_MCQBA_rc') `list_posts'
} // END OF QUIETLY STATEMENT	
} // END OF dataset FOR-LOOP
frame `memhold': save "`resultsfile'", replace		

* END OF DO-FILE