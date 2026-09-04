/***********************************************************************************************************************************
                                               RUNS SIMULATION STUDY III FOR MONTE CARLO QBA 
											SCENARIO C: Y BINARY, X BINARY, SINGLE U=(U1) BINARY

ARGUMENTS:
		samplesize       - SAMPLE SIZE OF EACH SIMULATED DATASET
		numsimdatasets   - NUMBER OF SIMULATED DATASETS
		numMCreps        - NUMBER OF MONTE CARLO REPLICATIONS
		priorlevel		 - LEVEL OF PRECISION OF BIAS PARAMETERS' PRIORS
		uniform			 - PRIOR FOR eta IS UNIFORM DISTRIBUTION FOR (uniform=1); OTHERWISE GAMMA DISTRIBUTION ON 1/etasq
		strength_UC      - MULTIPLIER OF ASSOCIATIONS BETWEEN U AND C (0 INDEPENDENCE BETWEEN U AND C, 1 ASSOCIATIONS OBSERVED IN REAL DATASET,
		                   2 DOUBLE THE OBSERVED ASSOCIATIONS)
		pilevel          - LEVEL OF THE PREVALENCE OF U (0=~6%, 1=~20%, 2=~40%)
		
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
tempname memhold alphax_hyperp betau_hyperp pi_hyperp

	* MACROLISTS
local exposurelist "xbin"
local confounderlist `"ccon cbin _Icnom_1 _Icnom_2"'
local parameterlist `"alphastar_xbin pistar betastar_ubin"'

* HYPERPARAMETERS OF THE PRIOR DISTRIBUTIONS FOR THE BIAS PARAMETERS (VALUES SUPPLIED BY EMILY)
	* UNIFORM DISTRIBUTION OR BETA DISTRIBUTION FOR pi
	* NORMAL DISTRIBUTIONS FOR COEFFCIENTS alpha_xbin AND beta_ubin
		* MEANS SET TO TRUE VALUES; SDs SET TO DIFFERENT VALUES	 
local mean_alpha_xbin 1.178
local mean_beta_ubin 1.5
if `priorlevel'==2 {
	noisily di "Applying Monte Carlo QBA using less informative prior distributions"
	local priortxt "medium"
}
else if `priorlevel'==3 {
	noisily di "Applying Monte Carlo QBA using very informative prior distributions"
	local priortxt "most"
	local sd_alpha_xbin 0.06010
	local sd_beta_ubin 0.07653
	
	if `pilevel'==0 {
		noisily di "Prevalence of U is 6.5%"

		if `strength_UC'==0 {
			* pi=6.54%
			local luni_pi 0.5852
			local uuni_pi 0.7228	
		}
		else if `strength_UC'==1 {
			* pi=6.60%
			local luni_pi 0.5905
			local uuni_pi 0.7295	
		}
		else if `strength_UC'==2 {
			* pi=6.45%
			local luni_pi 0.5771
			local uuni_pi 0.7129	
		}
	}
	else if `pilevel'==1 {
		noisily di "Prevalence of U is 20.5%"

		if `strength_UC'==0 {
			* pi=20.47%
			local luni_pi 0.1832
			local uuni_pi 0.2262	
		}
		else if `strength_UC'==1 {
			* pi=20.52%
			local luni_pi 0.1836
			local uuni_pi 0.2268	
		}
		else if `strength_UC'==2 {
			* pi=20.67%
			local luni_pi 0.1849
			local uuni_pi 0.2285	
		}
	}
	else if `pilevel'==2 {
		noisily di "Prevalence of U is 40%"

		if `strength_UC'==0 {
			* pi=39.89%
			local luni_pi 0.3569
			local uuni_pi 0.4409	
		}
		else if `strength_UC'==1 {
			* pi=40.13%
			local luni_pi 0.3591
			local uuni_pi 0.4435	
		}
		else if `strength_UC'==2 {
			* pi=39.59%
			local luni_pi 0.3542
			local uuni_pi 0.4376	
		}
		
	}
	local a_pi . 			// HYPERPARAMETERS OF A BETA DISTRIBUTION FOR pi
	local b_pi . 	
}  
matrix `alphax_hyperp' = (`mean_alpha_xbin', `sd_alpha_xbin')
matrix `betau_hyperp' = (`mean_beta_ubin', `sd_beta_ubin')
matrix `pi_hyperp' = (`luni_pi',`uuni_pi',`a_pi',`b_pi')

/**********************************
  MACROLISTS FOR FRAME STATEMENTS
*********************************/
* LOCAL MACROLISTS FOR THE VARIABLE NAMES
* FULL ANALYSIS MODEL
local vars_ybin `"b_ybin_xbin se_ybin_xbin ybin_eN ybin_erules ybin_econverged ybin_rc"'

* NAIVE MODEL
local vars_naive `"b_naive_xbin se_naive_xbin naive_eN naive_erules naive_econverged naive_rc"'

* MONTE CARLO BIAS ANALYSIS
local vars_bapd `"b_bapd_sys_betaxbin sd_bapd_sys_betaxbin lci_bapd_sys_betaxbin uci_bapd_sys_betaxbin b_bapd_all_betaxbin sd_bapd_all_betaxbin lci_bapd_all_betaxbin uci_bapd_all_betaxbin bapd_min_alphastar_xbin bapd_p025_alphastar_xbin bapd_p25_alphastar_xbin bapd_p50_alphastar_xbin bapd_p75_alphastar_xbin bapd_p975_alphastar_xbin bapd_max_alphastar_xbin bapd_min_pistar bapd_p025_pistar bapd_p25_pistar bapd_p50_pistar bapd_p75_pistar bapd_p975_pistar bapd_max_pistar bapd_min_betastar_ubin bapd_p025_betastar_ubin bapd_p25_betastar_ubin bapd_p50_betastar_ubin bapd_p75_betastar_ubin bapd_p975_betastar_ubin bapd_max_betastar_ubin bapd_timer bapd_num_analysisrc0 bapd_min_analysis_eN bapd_p25_analysis_eN bapd_p50_analysis_eN bapd_p75_analysis_eN bapd_max_analysis_eN"'

* LOCAL MACROLISTS FOR RESULTS TO BE POSTED TO THE FRAME
local list_roots `"ybin naive bapd"' 	
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
`"Simulation study III\Results\ScenarioC_logit_`samplesize'obs_`priortxt'inf_uniform`uniform'_`numMCreps'reps_strengthUC`strength_UC'_pilevel`pilevel'.dta"'

frame create `memhold' dataset strL(state) mean_alpha_xbin sd_alpha_xbin luni_pi uuni_pi mean_beta_ubin sd_beta_ubin ///
                                           bapd_MCQBA_rc `list_vars' 
							
forvalues dataset=1(1)`numsimdatasets' {
	noisily di ///
    "YbinXbinUbin; logit; strengthUC=`strength_UC'; pilevel=`pilevel'; `priortxt' informative prior; Uniform=`uniform'; `numMCreps' MCQBA reps: Processing dataset `dataset' of sample size `samplesize'"
	
quietly {   
	* RECORD THE rngstate FOR THIS LOOP
	local state = c(rngstate)
	
	* SIMULATE DATASET
	run "Data simulation\Scenario C - simulates a dataset.do" `samplesize' `strength_UC' `pilevel'	
	
	* GENERATE DUMMY VARIABLES OF cnom FOR ANALYSIS
	xi i.cnom
	
	* FULL ANALYSIS MODEL
	capture logit ybin xbin `confounderlist' ubin
	scalar ybin_rc = _rc
	if ybin_rc==0 {
		scalar b_ybin_xbin = _b[xbin]
		scalar se_ybin_xbin = _se[xbin]
		
		* OTHER
		matrix erules = e(rules)
		local numrows = rowsof(erules)
		local numcols = colsof(erules)
		scalar ybin_erules = 0						// INDICATOR IF PERFECT PREDICTION HAS OCCURRED
		forvalues row=1(1)`numrows' {
			forvalues col=1(1)`numcols' {
				di "row=`row'; col=`col'"
				if erules[`row',`col'] !=0 {
					scalar ybin_erules = 1
					continue, break
				}
			}
		}
		scalar ybin_eN = e(N)
		scalar ybin_econverged = e(converged)
	}
	else {			
		foreach item of local vars_ybin {
			scalar `item' = .
		}		
	}

	* NAIVE MODEL
	capture logit ybin xbin `confounderlist'
	scalar naive_rc = _rc
	if naive_rc==0 {
		scalar b_naive_xbin = _b[xbin]
		scalar se_naive_xbin = _se[xbin]
		
		* OTHER
		matrix erules = e(rules)
		local numrows = rowsof(erules)
		local numcols = colsof(erules)
		scalar naive_erules = 0					
		forvalues row=1(1)`numrows' {
			forvalues col=1(1)`numcols' {
				di "row=`row'; col=`col'"
				if erules[`row',`col'] !=0 {
					scalar naive_erules = 1
					continue, break
				}
			}
		}
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
	capture run "Fits MCQBA\Scenario C - fits MCQBA.do" "`exposurelist'" "`confounderlist'" ///
																	     "`alphax_hyperp'" "`pi_hyperp'" "`betau_hyperp'" ///
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
	
	frame post `memhold' (`dataset') ("`state'") (`mean_alpha_xbin') (`sd_alpha_xbin') (`luni_pi') (`uuni_pi') ///
	                     (`mean_beta_ubin') (`sd_beta_ubin') (`bapd_MCQBA_rc') `list_posts'
						 
} // END OF QUIETLY STATEMENT	
} // END OF dataset FOR-LOOP
frame `memhold': save "`resultsfile'", replace		

* END OF DO-FILE