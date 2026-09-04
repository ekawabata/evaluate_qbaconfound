/***********************************************************************************************************************************
                                                 RUNS SIMULATION STUDY VI FOR MONTE CARLO QBA 
											SCENARIO A: Y BINARY, X BINARY, SINGLE U=(U1) CONTINUOUS

ARGUMENTS:
		samplesize       - SAMPLE SIZE OF EACH SIMULATED DATASET
		numsimdatasets   - NUMBER OF SIMULATED DATASETS
		numMCreps        - NUMBER OF MONTE CARLO REPLICATIONS
		priorlevel		 - LEVEL OF PRECISION OF BIAS PARAMETERS' PRIORS
		uniform			 - PRIOR FOR eta IS UNIFORM DISTRIBUTION FOR (uniform=1); OTHERWISE GAMMA DISTRIBUTION ON 1/etasq
		
RETURNS: 
	DATASET CONTAINING RESULTS STORED BY POSTFILE; STORED UNDER FILENAME OF MACRO resultsfile
	STORED IN FOLDER "Results"
***********************************************************************************************************************************/
args samplesize numsimdatasets numMCreps priorlevel uniform

noisily di "Analysing `numsimdatasets' simulated datasets; each of `samplesize' observations"
noisily di "Applying Monte Carlo QBA to generate a frequency distribution of `numMCreps' estimates"

* DROP ALL SCALARS
scalar drop _all

/****************************************************
  MACROS AND TEMPORARY NAMES, FILES, AND VARIABLES
***************************************************/
tempfile mcqba_replications
tempname memhold alphax_hyperp betau_hyperp eta_hyperp

* MACROLISTS
local exposurelist "xbin"
local confounderlist `"ccon cbin _Icnom_1 _Icnom_2"'
local parameterlist `"alphastar_xbin etastar betastar_ucon"'

* HYPERPARAMETERS OF THE PRIOR DISTRIBUTIONS FOR THE BIAS PARAMETERS
	* UNIFORM DISTRIBUTION FOR eta1 AND eta2, OR GAMMA DISTRIBUTION FOR 1/eta1sq AND 1/eta2sq
	* NORMAL DISTRIBUTIONS FOR COEFFCIENTS 
		* FOR PRIOR LEVELS 2 AND 3: MEANS SET TO TRUE VALUES; SDs SET TO DIFFERENT VALUES
		* FOR PRIOR LEVELS 4 AND 5: MEANS SET TO VALUES SHIFTED BY 30%; SDs SET SET TO DIFFERENT VALUES 
if `priorlevel'==2 {
	noisily di "Applying Monte Carlo QBA using accurate, informative prior distributions"
	local priortxt "medium"	
	* alpha_xbin
	local mean_alpha_xbin 6.636	
	local sd_alpha_xbin 1.016
	* beta_ucon
	local mean_beta_ucon 0.06588
	local sd_beta_ucon 0.01008 
	* eta
	local luni_eta 6.516
	local uuni_eta 12.53
	local shape_invetasq 13.19
	local scale_invetasq = 1/1196	
}
else if `priorlevel'==3 {
	noisily di "Applying Monte Carlo QBA using accurate, very informative prior distributions"
	local priortxt "most"	   
	* alpha_xbin
	local mean_alpha_xbin 6.636	
	local sd_alpha_xbin 0.3386
	* beta_ucon
	local mean_beta_ucon 0.06588
	local sd_beta_ucon 0.003361 
	* eta
	local luni_eta 8.521
	local uuni_eta 10.53
	local shape_invetasq .
	local scale_invetasq .
}
else if `priorlevel'==4 {
	noisily di "Applying Monte Carlo QBA using inaccurate (30% below truth), informative prior distributions"
	local priortxt "inaccmedium"
	* alpha_xbin
	local mean_alpha_xbin = 0.7*6.636	
	local sd_alpha_xbin 1.016
	* beta_ucon
	local mean_beta_ucon = 0.7*0.06588
	local sd_beta_ucon 0.01008 
	* eta
	local acc_luni_eta 6.516
	local acc_uuni_eta 12.53
	local mean_eta = 0.7*9.523
	local luni_eta = `mean_eta' - 0.5*(`acc_uuni_eta'-`acc_luni_eta')
	local uuni_eta = `mean_eta' + 0.5*(`acc_uuni_eta'-`acc_luni_eta')
	local shape_invetasq .
	local scale_invetasq .	
}
else if `priorlevel'==5 {
	noisily di "Applying Monte Carlo QBA using inaccurate (30% below truth), less informative (sd double size of informative) priors"
	local priortxt "inaccweak"
	* alpha_xbin
	local mean_alpha_xbin = 0.7*6.636	
	local sd_alpha_xbin = 2*1.016
	* beta_ucon
	local mean_beta_ucon = 0.7*0.06588
	local sd_beta_ucon = 2*0.01008 
	* eta
	local acc_luni_eta 6.516
	local acc_uuni_eta 12.53
	local mean_eta = 0.7*9.523
	local luni_eta = `mean_eta' - (`acc_uuni_eta'-`acc_luni_eta')
	local uuni_eta = `mean_eta' + (`acc_uuni_eta'-`acc_luni_eta')
	local shape_invetasq .
	local scale_invetasq .	
}
matrix `alphax_hyperp' = (`mean_alpha_xbin', `sd_alpha_xbin')
matrix `betau_hyperp' = (`mean_beta_ucon', `sd_beta_ucon')
matrix `eta_hyperp' = (`luni_eta',`uuni_eta',`shape_invetasq',`scale_invetasq')

/**********************************
  MACROLISTS FOR FRAME STATEMENTS
*********************************/
* LOCAL MACROLISTS FOR THE VARIABLE NAMES
local vars_bapd `"b_bapd_sys_betaxbin sd_bapd_sys_betaxbin lci_bapd_sys_betaxbin uci_bapd_sys_betaxbin b_bapd_all_betaxbin sd_bapd_all_betaxbin lci_bapd_all_betaxbin uci_bapd_all_betaxbin bapd_min_alphastar_xbin bapd_p025_alphastar_xbin bapd_p25_alphastar_xbin bapd_p50_alphastar_xbin bapd_p75_alphastar_xbin bapd_p975_alphastar_xbin bapd_max_alphastar_xbin bapd_min_etastar bapd_p025_etastar bapd_p25_etastar bapd_p50_etastar bapd_p75_etastar bapd_p975_etastar bapd_max_etastar bapd_min_betastar_ucon bapd_p025_betastar_ucon bapd_p25_betastar_ucon bapd_p50_betastar_ucon bapd_p75_betastar_ucon bapd_p975_betastar_ucon bapd_max_betastar_ucon bapd_timer bapd_num_analysisrc0 bapd_min_analysis_eN bapd_p25_analysis_eN bapd_p50_analysis_eN bapd_p75_analysis_eN bapd_max_analysis_eN"'

* LOCAL MACROLISTS FOR RESULTS TO BE POSTED TO THE FRAME
local list_roots `"bapd"' 	
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

/**************************************
   CODE TO RUN THE SIMULATION STUDY
*************************************/
local resultsfile ///
`"Simulation study VI\Results\ScenarioA_logit_`samplesize'obs_`priortxt'inf_uniform`uniform'_`numMCreps'reps_strengthUC1.dta"'

frame create `memhold' dataset strL(state) mean_alpha_xbin sd_alpha_xbin luni_eta uuni_eta mean_beta_ucon sd_beta_ucon ///
                                           bapd_MCQBA_rc `list_vars' 
							
forvalues dataset=1(1)`numsimdatasets' {
	noisily di ///
    "YbinXbinUcon; logit; `priortxt' informative prior; Uniform=`uniform'; `numMCreps' MCQBA reps: Processing dataset `dataset' of sample size `samplesize'"
	
quietly {   
	* RECORD THE rngstate FOR THIS LOOP
	local state = c(rngstate)
	
	* IMPORT DATASET
	import delimited "Simulation study I\Scenario A\Data\Dataset_`dataset'.csv", clear
	
	* GENERATE DUMMY VARIABLES OF cnom FOR ANALYSIS
	xi i.cnom	
	
	/***********************************************
	 APPLY MONTE CARLO QBA TO SIMULATED DATASET 
    ***********************************************/
	timer clear 1
	timer on 1
	capture run "Fits MCQBA\Scenario A - fits MCQBA.do" "`exposurelist'" "`confounderlist'" "`alphax_hyperp'" ///
																		"`eta_hyperp'" "`betau_hyperp'" ///
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
		centile betahat_xbin, centile(50)
		scalar b_bapd_sys_betaxbin = r(c_1)
		centile betahat_xbin, centile(2.5)
		scalar lci_bapd_sys_betaxbin = r(c_1)
		centile betahat_xbin, centile(97.5)
		scalar uci_bapd_sys_betaxbin = r(c_1)
		summarize betahat_xbin
		scalar sd_bapd_sys_betaxbin	= r(sd)
		
			* RESULTS ACCOUNTING FOR UNCERTAINTY ABOUT SYSTEMATIC ERROR AND RANDOM SAMPLING
		centile betastar_xbin, centile(50)
		scalar b_bapd_all_betaxbin = r(c_1)
		centile betastar_xbin, centile(2.5)
		scalar lci_bapd_all_betaxbin = r(c_1)
		centile betastar_xbin, centile(97.5)
		scalar uci_bapd_all_betaxbin = r(c_1)
		summarize betastar_xbin 
		scalar sd_bapd_all_betaxbin = r(sd)	
	}
	else {
		foreach item of local vars_bapd {
			scalar `item' = .
		}			
	}	
	
	frame post `memhold' (`dataset') ("`state'") (`mean_alpha_xbin') (`sd_alpha_xbin') (`luni_eta') (`uuni_eta') ///
	                     (`mean_beta_ucon') (`sd_beta_ucon') (`bapd_MCQBA_rc') `list_posts'
						 
} // END OF QUIETLY STATEMENT	
} // END OF dataset FOR-LOOP
frame `memhold': save "`resultsfile'", replace		

* END OF DO-FILE