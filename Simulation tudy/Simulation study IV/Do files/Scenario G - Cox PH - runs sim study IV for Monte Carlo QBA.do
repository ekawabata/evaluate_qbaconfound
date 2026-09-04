/***********************************************************************************************************************************
                                               RUNS SIMULATION STUDY IV FOR MONTE CARLO QBA 
											SCENARIO G: Y SURVIVAL, X CONTINUOUS, SINGLE U=(CONTINUOUS U1)

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
local confounderlist `"ccon cbin _Icnom_2 _Icnom_3 _Icnom_4"'
local parameterlist `"alphastar_xcon etastar betastar_ucon"'

* HYPERPARAMETERS OF THE PRIOR DISTRIBUTIONS FOR THE BIAS PARAMETERS (VALUES SUPPLIED BY EMILY)
	* UNIFORM DISTRIBUTION FOR eta
	* NORMAL DISTRIBUTIONS FOR COEFFCIENTS alpha_xcon AND beta_ucon
		* MEANS SET TO TRUE VALUES; SDs SET TO DIFFERENT VALUES	 
local mean_alpha_xcon -7.8	
local mean_beta_ucon 0.0400
if `priorlevel'==2 {
	noisily di "Applying Monte Carlo QBA using less informative prior distributions"
	local priortxt "medium"
	local sd_alpha_xcon 1.194
	local sd_beta_ucon 0.006122
	local luni_eta 6.842
	local uuni_eta 13.16
	local shape_invetasq .
	local scale_invetasq .	
}
else if `priorlevel'==3 {
	noisily di "Applying Monte Carlo QBA using very informative prior distributions"
	local priortxt "most"
	local sd_alpha_xcon 0.3980
	local sd_beta_ucon 0.002041 
	local luni_eta 8.947
	local uuni_eta 11.05
	local shape_invetasq .
	local scale_invetasq .
}  
matrix `alphax_hyperp' = (`mean_alpha_xcon', `sd_alpha_xcon')
matrix `betau_hyperp' = (`mean_beta_ucon', `sd_beta_ucon')
matrix `eta_hyperp' = (`luni_eta',`uuni_eta',`shape_invetasq',`scale_invetasq')

/**********************************
  MACROLISTS FOR FRAME STATEMENTS
*********************************/
* LOCAL MACROLISTS FOR THE VARIABLE NAMES
local vars_bapd `"b_bapd_sys_betaxcon sd_bapd_sys_betaxcon lci_bapd_sys_betaxcon uci_bapd_sys_betaxcon b_bapd_all_betaxcon sd_bapd_all_betaxcon lci_bapd_all_betaxcon uci_bapd_all_betaxcon bapd_min_alphastar_xcon bapd_p025_alphastar_xcon bapd_p25_alphastar_xcon bapd_p50_alphastar_xcon bapd_p75_alphastar_xcon bapd_p975_alphastar_xcon bapd_max_alphastar_xcon bapd_min_etastar bapd_p025_etastar bapd_p25_etastar bapd_p50_etastar bapd_p75_etastar bapd_p975_etastar bapd_max_etastar bapd_min_betastar_ucon bapd_p025_betastar_ucon bapd_p25_betastar_ucon bapd_p50_betastar_ucon bapd_p75_betastar_ucon bapd_p975_betastar_ucon bapd_max_betastar_ucon bapd_timer bapd_num_analysisrc0 bapd_min_analysis_eN bapd_p25_analysis_eN bapd_p50_analysis_eN bapd_p75_analysis_eN bapd_max_analysis_eN"'

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
local resultsfile `"Simulation study IV\Results\ScenarioG_CoxPH_`samplesize'obs_`priortxt'inf_uniform`uniform'_`numMCreps'reps.dta"'

frame create `memhold' dataset strL(state) bapd_MCQBA_rc `list_vars' 
							
forvalues dataset=1(1)`numsimdatasets' {
		noisily di "Scenario G; Cox PH; `priortxt' informative prior; Uniform=`uniform'; `numMCreps' MCQBA reps: Processing dataset `dataset' of sample size `samplesize'"
	
quietly {   
	* RECORD THE rngstate FOR THIS LOOP
	local state = c(rngstate)
	
	* IMPORT DATASET
	import delimited "Simulation study IV\Data\Dataset_`dataset'.csv", clear
	
	* GENERATE DUMMY VARIABLES OF cnom FOR ANALYSIS
	xi i.cnom
	
	/*********************************************
	 APPLY MONTE CARLO QBA TO SIMULATED DATASET 
    *********************************************/
	timer clear 1
	timer on 1
	capture run "Fits MCQBA\Scenario G - fits MCQBA.do" /// 
									"`exposurelist'" "`confounderlist'" "`alphax_hyperp'" ///
	                                "`eta_hyperp'" "`betau_hyperp'" `numMCreps' `mcqba_replications' `uniform'
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
	
	frame post `memhold' (`dataset') ("`state'") (`bapd_MCQBA_rc') `list_posts'
						 
} // END OF QUIETLY STATEMENT	
} // END OF dataset FOR-LOOP
frame `memhold': save "`resultsfile'", replace		

* END OF DO-FILE