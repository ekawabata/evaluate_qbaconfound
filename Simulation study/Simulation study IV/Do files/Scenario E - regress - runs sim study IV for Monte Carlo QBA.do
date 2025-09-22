/***********************************************************************************************************************************
                                               RUNS SIMULATION STUDY IV FOR MONTE CARLO QBA 
											SCENARIO E: Y CONTINUOUS, X NOMINAL, MULTIPLE U=(CONTINUOUS U1, BINARY U2)

ARGUMENTS:
		samplesize       - SAMPLE SIZE OF EACH SIMULATED DATASET
		numsimdatasets   - NUMBER OF SIMULATED DATASETS
		numMCreps        - NUMBER OF MONTE CARLO REPLICATIONS
		priorlevel		 - LEVEL OF PRECISION OF BIAS PARAMETERS' PRIORS
		uniform			 - PRIOR FOR eta IS UNIFORM DISTRIBUTION FOR (uniform=1); OTHERWISE GAMMA DISTRIBUTION ON 1/etasq
		strength_UC      - MULTIPLIER OF ASSOCIATIONS BETWEEN U AND C (0 INDEPENDENCE BETWEEN U AND C, 1 ASSOCIATIONS OBSERVED IN REAL DATASET,
		                   2 DOUBLE THE OBSERVED ASSOCIATIONS)
		pilevel          - PREVALENCE OF U2 (1:Pr(U2=1)=59.77%)
		
RETURNS: 
	DATASET CONTAINING RESULTS STORED BY POSTFILE; STORED UNDER FILENAME OF MACRO resultsfile
	STORED IN FOLDER "Results"
***********************************************************************************************************************************/
args samplesize numsimdatasets numMCreps priorlevel uniform strength_UC pilevel

noisily di "Analysing `numsimdatasets' simulated datasets; each of `samplesize' observations"
noisily di "Applying Monte Carlo QBA to generate a frequency distribution of `numMCreps' estimates"
noisily di "Multiply the observed U-C associations by `strength_UC'"

* DROP ALL SCALARS
scalar drop _all

/****************************************************
   MACROS AND TEMPORARY NAMES, FILES, AND VARIABLES
***************************************************/
tempfile mcqba_replications
tempname memhold alphax_hyperp eta1_hyperp pi2_hyperp betau_hyperp 

	* MACROLISTS
local exposurelist "xnom0 xnom1"
local confounderlist `"ccon cbin cnom1 cnom2"'
local parameterlist `"alpha1star_xnom0 alpha1star_xnom1 eta1star betastar_ucon1 alpha2star_xnom0 alpha2star_xnom1 pi2star betastar_ubin2"' 

* HYPERPARAMETERS OF THE PRIOR DISTRIBUTIONS FOR THE BIAS PARAMETERS (VALUES SUPPLIED BY EMILY)
	* UNIFORM DISTRIBUTIONS FOR eta1 AND pi2
	* NORMAL DISTRIBUTIONS FOR ALL COEFFCIENTS
		* MEANS SET TO TRUE VALUES; SDs SET TO DIFFERENT VALUES	 
local mean_alpha1_xnom0 -7.492
local mean_alpha1_xnom1 -4.175
local mean_alpha2_xnom0 -1.172
local mean_alpha2_xnom1 -0.6502
local mean_beta_ucon1 0.1299
local mean_beta_ubin2 0.5000
if `priorlevel'==2 {
	noisily di "Applying Monte Carlo QBA using less informative prior distributions"
	local priortxt "medium"
	local sd_alpha1_xnom0 1.147
	local sd_alpha1_xnom1 0.6390
	local sd_alpha2_xnom0 0.1794
	local sd_alpha2_xnom1 0.09951 
	local sd_beta_ucon1 0.01988
	local sd_beta_ubin2 0.07653
	local luni_eta1 6.537
	local uuni_eta1 12.57
	local shape_inveta1sq .
	local scale_inveta1sq .
	if `pilevel'==1 {
		local luni_pi2 0.4089
		local uuni_pi2 0.7864       
		local a_pi2 .
		local b_pi2 .
	}
}
else if `priorlevel'==3 {
	noisily di "Applying Monte Carlo QBA using very informative prior distributions"
	local priortxt "most"
	local sd_alpha1_xnom0 0.3822
	local sd_alpha1_xnom1 0.2130
	local sd_alpha2_xnom0 0.05981 
	local sd_alpha2_xnom1 0.03317
	local sd_beta_ucon1 0.006628
	local sd_beta_ubin2 0.02551
	local luni_eta1 8.548
	local uuni_eta1 10.56
	local shape_inveta1sq .
	local scale_inveta1sq .
	if `pilevel'==1 {
		local luni_pi2 0.5348
		local uuni_pi2 0.6606
		local a_pi2 .
		local b_pi2 .
	}
}  
matrix `alphax_hyperp' = (`mean_alpha1_xnom0', `sd_alpha1_xnom0', `mean_alpha1_xnom1', `sd_alpha1_xnom1' \ ///
                          `mean_alpha2_xnom0', `sd_alpha2_xnom0', `mean_alpha2_xnom1', `sd_alpha2_xnom1')
matrix `eta1_hyperp' = (`luni_eta1',`uuni_eta1',`shape_inveta1sq',`scale_inveta1sq')
matrix `pi2_hyperp' = (`luni_pi2',`uuni_pi2',`a_pi2',`b_pi2')
matrix `betau_hyperp' = (`mean_beta_ucon1', `sd_beta_ucon1' \ ///
                         `mean_beta_ubin2', `sd_beta_ubin2')

/*********************************
  MACROLISTS FOR FRAME STATEMENTS
*********************************/
* FULL ANALYSIS MODEL
local vars_ycon `"b_ycon_xnom0 b_ycon_xnom1 se_ycon_xnom0 se_ycon_xnom1 ycon_eN ycon_rc"'

* NAIVE MODEL
local vars_naive `"b_naive_xnom0 b_naive_xnom1 se_naive_xnom0 se_naive_xnom1 naive_eN naive_rc"'

* BIAS ANALYSIS
local vars_bapd1 `"b_bapd_sys_betaxnom0 b_bapd_sys_betaxnom1 sd_bapd_sys_betaxnom0 sd_bapd_sys_betaxnom1 lci_bapd_sys_betaxnom0 lci_bapd_sys_betaxnom1 uci_bapd_sys_betaxnom0 uci_bapd_sys_betaxnom1 b_bapd_all_betaxnom0 b_bapd_all_betaxnom1 sd_bapd_all_betaxnom0 sd_bapd_all_betaxnom1 lci_bapd_all_betaxnom0 lci_bapd_all_betaxnom1 uci_bapd_all_betaxnom0 uci_bapd_all_betaxnom1 bapd_min_alpha1star_xnom0 bapd_p025_alpha1star_xnom0 bapd_p25_alpha1star_xnom0 bapd_p50_alpha1star_xnom0 bapd_p75_alpha1star_xnom0 bapd_p975_alpha1star_xnom0 bapd_max_alpha1star_xnom0 bapd_min_alpha1star_xnom1 bapd_p025_alpha1star_xnom1 bapd_p25_alpha1star_xnom1 bapd_p50_alpha1star_xnom1 bapd_p75_alpha1star_xnom1 bapd_p975_alpha1star_xnom1 bapd_max_alpha1star_xnom1 bapd_min_eta1star bapd_p025_eta1star bapd_p25_eta1star bapd_p50_eta1star bapd_p75_eta1star bapd_p975_eta1star bapd_max_eta1star"' 
local vars_bapd2 `"bapd_min_alpha2star_xnom0 bapd_p025_alpha2star_xnom0 bapd_p25_alpha2star_xnom0 bapd_p50_alpha2star_xnom0 bapd_p75_alpha2star_xnom0 bapd_p975_alpha2star_xnom0 bapd_max_alpha2star_xnom0 bapd_min_alpha2star_xnom1 bapd_p025_alpha2star_xnom1 bapd_p25_alpha2star_xnom1 bapd_p50_alpha2star_xnom1 bapd_p75_alpha2star_xnom1 bapd_p975_alpha2star_xnom1 bapd_max_alpha2star_xnom1 bapd_min_pi2star bapd_p025_pi2star bapd_p25_pi2star bapd_p50_pi2star bapd_p75_pi2star bapd_p975_pi2star bapd_max_pi2star  bapd_min_betastar_ucon1 bapd_p025_betastar_ucon1 bapd_p25_betastar_ucon1 bapd_p50_betastar_ucon1 bapd_p75_betastar_ucon1 bapd_p975_betastar_ucon1 bapd_max_betastar_ucon1 bapd_min_betastar_ubin2 bapd_p025_betastar_ubin2 bapd_p25_betastar_ubin2 bapd_p50_betastar_ubin2 bapd_p75_betastar_ubin2 bapd_p975_betastar_ubin2 bapd_max_betastar_ubin2 bapd_timer bapd_num_analysisrc0 bapd_min_analysis_eN bapd_p25_analysis_eN bapd_p50_analysis_eN bapd_p75_analysis_eN bapd_max_analysis_eN"'

* LOCAL MACROLISTS FOR RESULTS TO BE POSTED TO THE FRAME
local list_roots `"ycon naive bapd1 bapd2"' 	
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
* RESULTS STORED HERE
local resultsfile `"Simulation study IV\Results\ScenarioE_regress_`samplesize'obs_`priortxt'inf_uniform`uniform'_`numMCreps'reps_strengthUC`strength_UC'_pilevel`pilevel'.dta"'

frame create `memhold' dataset strL(state) bapd_MCQBA_rc `list_vars' 
forvalues dataset=1(1)`numsimdatasets' {
	noisily di ///
	"Scenario E; regress; strengthUC=`strength_UC'; pilevel=`pilevel'; Uniform=`uniform'; `numMCreps' MCQBA reps: Processing dataset `dataset' of sample size `samplesize'"
	
quietly {   
	* RECORD THE rngstate FOR THIS LOOP
	local state = c(rngstate)
	
	* SIMULATE A DATASET
	run "Data simulation\Scenario E - simulates a dataset.do" `samplesize' `strength_UC' `pilevel'
	
	* FULL ANALYSIS MODEL
	capture regress ycon `exposurelist' `confounderlist' ucon1 ubin2
	scalar ycon_rc = _rc
	if ycon_rc==0 {
		scalar b_ycon_xnom0 = _b[xnom0]
		scalar b_ycon_xnom1 = _b[xnom1]
		scalar se_ycon_xnom0 = _se[xnom0]
		scalar se_ycon_xnom1 = _se[xnom1]
		scalar ycon_eN = e(N)
	}
	else {			
		foreach item of local vars_ycon {
			scalar `item' = .
		}		
	}

	* NAIVE MODEL
	capture regress ycon `exposurelist' `confounderlist'
	scalar naive_rc = _rc
	if naive_rc==0 {
		scalar b_naive_xnom0 = _b[xnom0]
		scalar b_naive_xnom1 = _b[xnom1]
		scalar se_naive_xnom0 = _se[xnom0]
		scalar se_naive_xnom1 = _se[xnom1]
		scalar naive_eN = e(N)
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
	capture run "Fits MCQBA\Scenario E - fits MCQBA.do" ///
											"`exposurelist'" "`confounderlist'" ///
											`alphax_hyperp' `eta1_hyperp' `pi2_hyperp' `betau_hyperp' ///
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

		* SUMMARIZE FREQUENCY DISTRIBUTION OF MONTE CARLO ESTIMATES	
		foreach exposure of local exposurelist {
			* RESULTS ACCOUNTING FOR UNCERTAINTY ABOUT SYSTEMATIC ERROR ONLY
			centile betahat_`exposure', centile(50)
			scalar b_bapd_sys_beta`exposure' = r(c_1)
			centile betahat_`exposure', centile(2.5)
			scalar lci_bapd_sys_beta`exposure' = r(c_1)
			centile betahat_`exposure', centile(97.5)
			scalar uci_bapd_sys_beta`exposure' = r(c_1)
			summarize betahat_`exposure'
			scalar sd_bapd_sys_beta`exposure'	= r(sd)
			
			* RESULTS ACCOUNTING FOR UNCERTAINTY ABOUT SYSTEMATIC ERROR AND RANDOM SAMPLING
			centile betastar_`exposure', centile(50)
			scalar b_bapd_all_beta`exposure' = r(c_1)
			centile betastar_`exposure', centile(2.5)
			scalar lci_bapd_all_beta`exposure' = r(c_1)
			centile betastar_`exposure', centile(97.5)
			scalar uci_bapd_all_beta`exposure' = r(c_1)
			summarize betastar_`exposure' 
			scalar sd_bapd_all_beta`exposure' = r(sd)	
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