/***********************************************************************************************************************************
                                               RUNS SIMULATION STUDY IV FOR MONTE CARLO QBA 
											SCENARIO F: Y NOMINAL, X CONTINUOUS, SINGLE U=(CONTINUOUS U1)

ARGUMENTS:
		samplesize       - SAMPLE SIZE OF EACH SIMULATED DATASET
		numsimdatasets   - NUMBER OF SIMULATED DATASETS
		numMCreps        - NUMBER OF MONTE CARLO REPLICATIONS
		priorlevel		 - LEVEL OF PRECISION OF BIAS PARAMETERS' PRIORS
		uniform			 - PRIOR FOR eta IS UNIFORM DISTRIBUTION FOR (uniform=1); OTHERWISE GAMMA DISTRIBUTION ON 1/etasq
		strength_UC      - MULTIPLIER OF ASSOCIATIONS BETWEEN U AND C (0 INDEPENDENCE BETWEEN U AND C, 1 ASSOCIATIONS OBSERVED IN REAL DATASET,
		                   2 DOUBLE THE OBSERVED ASSOCIATIONS)
		
RETURNS: 
	DATASET CONTAINING RESULTS STORED BY POSTFILE; STORED UNDER FILENAME OF MACRO resultsfile
	STORED IN FOLDER "Results"
***********************************************************************************************************************************/
args samplesize numsimdatasets numMCreps priorlevel uniform strength_UC

noisily di "Analysing `numsimdatasets' simulated datasets; each of `samplesize' observations"
noisily di "Applying Monte Carlo QBA to generate a frequency distribution of `numMCreps' estimates"
noisily di "Multiply the observed U-C associations by `strength_UC'"

* DROP ALL SCALARS
scalar drop _all

/****************************************************
  MACROS AND TEMPORARY NAMES, FILES, AND VARIABLES
***************************************************/
tempfile mcqba_replications
tempname memhold alphax_hyperp betau_hyperp eta_hyperp

	* MACROLISTS
local exposurelist "xcon"
local confounderlist `"ccon cbin cnom1 cnom2"'
local parameterlist `"alphastar_xcon etastar beta0star_ucon beta1star_ucon"'

* HYPERPARAMETERS OF THE PRIOR DISTRIBUTIONS FOR THE BIAS PARAMETERS (VALUES SUPPLIED BY EMILY)
	* UNIFORM DISTRIBUTION FOR eta
	* NORMAL DISTRIBUTIONS FOR COEFFCIENTS alpha_xcon AND beta_ucon
		* MEANS SET TO TRUE VALUES; SDs SET TO DIFFERENT VALUES	 
local mean_alpha_xcon 1.281	
local mean_beta0_ucon 0.1
local mean_beta1_ucon 0.1

if `priorlevel'==2 {
	noisily di "Applying Monte Carlo QBA using very informative prior distributions"
	local priortxt "medium"
	local sd_alpha_xcon .
	local sd_beta0_ucon .
	local sd_beta1_ucon .
	local luni_eta .
	local uuni_eta .
	local shape_invetasq .
	local scale_invetasq .	
}
else if `priorlevel'==3 {
	noisily di "Applying Monte Carlo QBA using very informative prior distributions"
	local priortxt "most"
	local sd_alpha_xcon 0.06536
	local sd_beta0_ucon 0.005102
	local sd_beta1_ucon 0.005102 
	local luni_eta  9.216
	local uuni_eta 11.38
	local shape_invetasq .
	local scale_invetasq .
}  
matrix `alphax_hyperp' = (`mean_alpha_xcon', `sd_alpha_xcon')
matrix `betau_hyperp' = (`mean_beta0_ucon', `sd_beta0_ucon' \ `mean_beta1_ucon', `sd_beta1_ucon')
matrix `eta_hyperp' = (`luni_eta',`uuni_eta', `shape_invetasq',`scale_invetasq')

/**********************************
  MACROLISTS FOR FRAME STATEMENTS
*********************************/
* FULL ANALYSIS MODEL
local vars_ynom `"b_ynom0_xcon se_ynom0_xcon b_ynom1_xcon se_ynom1_xcon ynom_eN ynom_econverged ynom_rc"'

* NAIVE MODEL
local vars_naive `"b_naive0_xcon se_naive0_xcon b_naive1_xcon se_naive1_xcon naive_eN naive_econverged naive_rc"'

* BIAS ANALYSIS
local vars_bapd_beta0 `"b_bapd_sys_beta0xcon sd_bapd_sys_beta0xcon lci_bapd_sys_beta0xcon uci_bapd_sys_beta0xcon b_bapd_all_beta0xcon sd_bapd_all_beta0xcon lci_bapd_all_beta0xcon uci_bapd_all_beta0xcon b_bapd_sys_beta1xcon sd_bapd_sys_beta1xcon lci_bapd_sys_beta1xcon uci_bapd_sys_beta1xcon b_bapd_all_beta1xcon sd_bapd_all_beta1xcon lci_bapd_all_beta1xcon uci_bapd_all_beta1xcon"'
local vars_bapd_other `"bapd_min_alphastar_xcon bapd_p025_alphastar_xcon bapd_p25_alphastar_xcon bapd_p50_alphastar_xcon bapd_p75_alphastar_xcon bapd_p975_alphastar_xcon bapd_max_alphastar_xcon bapd_min_etastar bapd_p025_etastar bapd_p25_etastar bapd_p50_etastar bapd_p75_etastar bapd_p975_etastar bapd_max_etastar bapd_min_beta0star_ucon bapd_p025_beta0star_ucon bapd_p25_beta0star_ucon bapd_p50_beta0star_ucon bapd_p75_beta0star_ucon bapd_p975_beta0star_ucon bapd_max_beta0star_ucon bapd_min_beta1star_ucon bapd_p025_beta1star_ucon bapd_p25_beta1star_ucon bapd_p50_beta1star_ucon bapd_p75_beta1star_ucon bapd_p975_beta1star_ucon bapd_max_beta1star_ucon bapd_timer bapd_num_analysisrc0 bapd_min_analysis_eN bapd_p25_analysis_eN bapd_p50_analysis_eN bapd_p75_analysis_eN bapd_max_analysis_eN"'

* LOCAL MACROLISTS FOR RESULTS TO BE POSTED TO THE FRAME
local list_roots `"ynom naive bapd_beta0 bapd_other"' 	
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
local resultsfile `"Simulation study IV\Results\ScenarioF_mlogit_`samplesize'obs_`priortxt'inf_uniform`uniform'_`numMCreps'reps_strengthUC`strength_UC'.dta"'

frame create `memhold' dataset strL(state) bapd_MCQBA_rc `list_vars' 
							
forvalues dataset=1(1)`numsimdatasets' {
	noisily di ///
	"Scenario F; mlogit; strengthUC=`strength_UC'; Uniform=`uniform'; `numMCreps' MCQBA reps: Processing dataset `dataset' of sample size `samplesize'"

quietly {   
	* RECORD THE rngstate FOR THIS LOOP
	local state = c(rngstate)
	
	* SIMULATE DATASET
	run "Data simulation\Scenario F - simulates a dataset.do" `samplesize' `strength_UC'
	
	* FULL ANALYSIS MODEL
	capture mlogit ynom `exposurelist' `confounderlist' ucon, baseoutcome(2)
	scalar ynom_rc = _rc
	if ynom_rc==0 {
		scalar b_ynom0_xcon = _b[0: xcon]
		scalar se_ynom0_xcon = _se[0: xcon]
		scalar b_ynom1_xcon = _b[1: xcon]
		scalar se_ynom1_xcon = _se[1: xcon]
				
		* OTHER
		scalar ynom_eN = e(N)
		scalar ynom_econverged = e(converged)
	}
	else {			
		foreach item of local vars_ynom {
			scalar `item' = .
		}		
	}

	* NAIVE MODEL
	capture mlogit ynom `exposurelist' `confounderlist', baseoutcome(2)
	scalar naive_rc = _rc
	if naive_rc==0 {
		scalar b_naive0_xcon = _b[0: xcon]
		scalar se_naive0_xcon = _se[0: xcon]
		scalar b_naive1_xcon = _b[1: xcon]
		scalar se_naive1_xcon = _se[1: xcon]
		
		* OTHER
		scalar naive_eN = e(N)
		scalar naive_econverged = e(converged)
	}
	else {			
		foreach item of local vars_naive {
			scalar `item' = .
		}		
	}
	
	/*********************************************
	 APPLY MONTE CARLO QBA TO SIMULATED DATASET 
    *********************************************/
	timer clear 1
	timer on 1
	capture run "Fits MCQBA\Scenario F - fits MCQBA.do" "`exposurelist'" "`confounderlist'" "`alphax_hyperp'" ///
	                                                    "`eta_hyperp'" "`betau_hyperp'" `numMCreps' `mcqba_replications' `uniform'
	local bapd_MCQBA_rc = _rc
	timer off 1
	timer list 1
	scalar bapd_timer = r(t1)

	if `bapd_MCQBA_rc'==0 {
		
		* PROCESS MCSA RESULTS
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

		forvalues k=0(1)1 {
			* RESULTS ACCOUNTING FOR UNCERTAINTY ABOUT SYSTEMATIC ERROR ONLY
			centile beta`k'hat_xcon, centile(50)
			scalar b_bapd_sys_beta`k'xcon = r(c_1)
			centile beta`k'hat_xcon, centile(2.5)
			scalar lci_bapd_sys_beta`k'xcon = r(c_1)
			centile beta`k'hat_xcon, centile(97.5)
			scalar uci_bapd_sys_beta`k'xcon = r(c_1)
			summarize beta`k'hat_xcon
			scalar sd_bapd_sys_beta`k'xcon	= r(sd)
			
			* RESULTS ACCOUNTING FOR UNCERTAINTY ABOUT SYSTEMATIC ERROR AND RANDOM SAMPLING
			centile beta`k'star_xcon, centile(50)
			scalar b_bapd_all_beta`k'xcon = r(c_1)
			centile beta`k'star_xcon, centile(2.5)
			scalar lci_bapd_all_beta`k'xcon = r(c_1)
			centile beta`k'star_xcon, centile(97.5)
			scalar uci_bapd_all_beta`k'xcon = r(c_1)
			summarize beta`k'star_xcon 
			scalar sd_bapd_all_beta`k'xcon = r(sd)	
		}
	}
	else {
		foreach item of local vars_bapd {
			scalar `item' = .
		}			
	}	
	
	frame post `memhold' (`dataset') ("`state'") (`bapd_MCQBA_rc') `list_posts'
						 
} // END OF QUIETLY STATEMENT	
} // END OF dataset FOR-LOOP
frame `memhold': save "`resultsfile'", replace		

* END OF DO-FILE