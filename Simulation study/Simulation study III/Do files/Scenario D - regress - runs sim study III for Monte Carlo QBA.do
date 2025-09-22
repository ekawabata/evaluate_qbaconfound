/***********************************************************************************************************************************
                                               RUNS SIMULATION STUDY III FOR MONTE CARLO QBA 
											SCENARIO D: Y CONTINUOUS, X CONTINUOUS, MULTIPLE U=(U1,U2) BINARY

ARGUMENTS:
		samplesize       - SAMPLE SIZE OF EACH SIMULATED DATASET
		numsimdatasets   - NUMBER OF SIMULATED DATASETS
		numMCreps        - NUMBER OF MONTE CARLO REPLICATIONS
		priorlevel		 - LEVEL OF PRECISION OF BIAS PARAMETERS' PRIORS
		uniform			 - PRIOR FOR eta IS UNIFORM DISTRIBUTION FOR (uniform=1); OTHERWISE GAMMA DISTRIBUTION ON 1/etasq
		strength_UC      - MULTIPLIER OF ASSOCIATIONS BETWEEN U AND C (0 INDEPENDENCE BETWEEN U AND C, 1 ASSOCIATIONS OBSERVED IN REAL DATASET,
		                   2 DOUBLE THE OBSERVED ASSOCIATIONS)
		pilevel          - PREVALENCES OF U1 AND U2 (1:Pr(U1=1)~14.1% AND Pr(U2=1)=18%)
		
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
tempname memhold alphax_hyperp pi_hyperp betau_hyperp   

	* MACROLISTS
local exposurelist "xcon"
local confounderlist `"ccon cbin _Icnom_1 _Icnom_2"'
local parameterlist `"alpha1star_xcon pi1star alpha2star_xcon pi2star betastar_ubin1 betastar_ubin2"'

* HYPERPARAMETERS OF THE PRIOR DISTRIBUTIONS FOR THE BIAS PARAMETERS (VALUES SUPPLIED BY EMILY)
	* UNIFORM DISTRIBUTION FOR eta
	* NORMAL DISTRIBUTIONS FOR COEFFCIENTS
		* MEANS SET TO TRUE VALUES; SDs SET TO DIFFERENT VALUES	 
local mean_alpha1_xcon 0.2500
local mean_alpha2_xcon 0.2550	// ESTIMATED FROM FITTING MODEL U2|X,C TO DATASET OF 50,000,000 OBS. ROUNDED TO NEAREST 3 D.P.
local mean_beta_ubin1 2.5
local mean_beta_ubin2 1.5

if `priorlevel'==2 {
	noisily di "Applying Monte Carlo QBA using less informative prior distributions"
	local priortxt "medium"
}
else if `priorlevel'==3 {
	noisily di "Applying Monte Carlo QBA using very informative prior distributions"
	local priortxt "most"
	local sd_alpha1_xcon 0.01276
	local sd_beta_ubin1 0.1276
	local sd_alpha2_xcon 0.01301
	local sd_beta_ubin2 0.07653
	if `pilevel'==1 {
		local luni_pi1 0.1262
		local uuni_pi1 0.1558
		local a_pi1 .
		local b_pi1 .
		local luni_pi2 0.161
		local uuni_pi2 0.199
		local a_pi2 .
		local b_pi2 .
	}
}  
matrix `alphax_hyperp' = (`mean_alpha1_xcon', `sd_alpha1_xcon' \ `mean_alpha2_xcon', `sd_alpha2_xcon')
matrix `betau_hyperp' = (`mean_beta_ubin1', `sd_beta_ubin1' \ `mean_beta_ubin2', `sd_beta_ubin2')
matrix `pi_hyperp' = (`luni_pi1',`uuni_pi1',`a_pi1',`b_pi1' \ ///
                      `luni_pi2',`uuni_pi2',`a_pi2',`b_pi2')

/*********************************
  MACROLISTS FOR FRAME STATEMENTS
*********************************/
* FULL ANALYSIS MODEL
local vars_ycon `"b_ycon_xcon se_ycon_xcon ycon_eN ycon_rc"'

* NAIVE MODEL
local vars_naive `"b_naive_xcon se_naive_xcon naive_eN naive_rc"'

* MONTE CARLO BIAS ANALYSIS
local vars_bapd1 `"b_bapd_sys_betaxcon sd_bapd_sys_betaxcon lci_bapd_sys_betaxcon uci_bapd_sys_betaxcon b_bapd_all_betaxcon sd_bapd_all_betaxcon lci_bapd_all_betaxcon uci_bapd_all_betaxcon bapd_timer bapd_num_analysisrc0 bapd_min_analysis_eN bapd_p25_analysis_eN bapd_p50_analysis_eN bapd_p75_analysis_eN bapd_max_analysis_eN"'

local vars_bapd2 `"bapd_min_alpha2star_xcon bapd_p025_alpha2star_xcon bapd_p25_alpha2star_xcon bapd_p50_alpha2star_xcon bapd_p75_alpha2star_xcon bapd_p975_alpha2star_xcon bapd_max_alpha2star_xcon bapd_min_pi2star bapd_p025_pi2star bapd_p25_pi2star bapd_p50_pi2star bapd_p75_pi2star bapd_p975_pi2star bapd_max_pi2star bapd_min_alpha1star_xcon bapd_p025_alpha1star_xcon bapd_p25_alpha1star_xcon bapd_p50_alpha1star_xcon bapd_p75_alpha1star_xcon bapd_p975_alpha1star_xcon bapd_max_alpha1star_xcon bapd_min_pi1star bapd_p025_pi1star bapd_p25_pi1star bapd_p50_pi1star bapd_p75_pi1star bapd_p975_pi1star bapd_max_pi1star bapd_min_betastar_ubin1 bapd_p025_betastar_ubin1 bapd_p25_betastar_ubin1 bapd_p50_betastar_ubin1 bapd_p75_betastar_ubin1 bapd_p975_betastar_ubin1 bapd_max_betastar_ubin1 bapd_min_betastar_ubin2 bapd_p025_betastar_ubin2 bapd_p25_betastar_ubin2 bapd_p50_betastar_ubin2 bapd_p75_betastar_ubin2 bapd_p975_betastar_ubin2 bapd_max_betastar_ubin2"'

* LOCAL MACROLISTS FOR RESULTS TO BE POSTED TO THE FRAME
local list_roots `"bapd1 bapd2 ycon naive"' 	
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
local resultsfile `"Simulation study III\Results\ScenarioD_regress_`samplesize'obs_`priortxt'inf_uniform`uniform'_`numMCreps'reps_strengthUC`strength_UC'_pilevel`pilevel'.dta"'

frame create `memhold' dataset strL(state) mean_alpha2_xcon sd_alpha2_xcon luni_pi2 uuni_pi2 a_pi2 b_pi2 mean_beta_ubin2 sd_beta_ubin2 ///
                                           mean_alpha1_xcon sd_alpha1_xcon luni_pi1 uuni_pi1 a_pi1 b_pi1 mean_beta_ubin1 sd_beta_ubin1  ///
										   bapd_MCQBA_rc `list_vars' 
forvalues dataset=1(1)`numsimdatasets' {
	noisily di ///
	"YconXconUbin; regress; strengthUC=`strength_UC'; pilevel=`pilevel'; Uniform=`uniform'; `numMCreps' MCQBA reps: Processing dataset `dataset' of sample size `samplesize'"
	
quietly {   
	* RECORD THE rngstate FOR THIS LOOP
	local state = c(rngstate)
	
	* SIMULATE DATASET
	run "Data simulation\Scenario D - simulates a dataset.do" `samplesize' `strength_UC' `pilevel'
	
	* GENERATE DUMMY VARIABLES OF cnom FOR ANALYSIS
	xi i.cnom
	
	* FULL ANALYSIS MODEL
	capture regress ycon xcon `confounderlist' ubin1 ubin2
	scalar ycon_rc = _rc
	if ycon_rc==0 {
		scalar b_ycon_xcon = _b[xcon]
		scalar se_ycon_xcon = _se[xcon]
		scalar ycon_eN = e(N)
	}
	else {			
		foreach item of local vars_ycon {
			scalar `item' = .
		}		
	}

	* NAIVE MODEL
	capture regress ycon xcon `confounderlist'
	scalar naive_rc = _rc
	if naive_rc==0 {
		scalar b_naive_xcon = _b[xcon]
		scalar se_naive_xcon = _se[xcon]
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
	capture run "Fits MCQBA\Scenario D - fits MCQBA.do" ///
		        "`exposurelist'" "`confounderlist'" ///
				`alphax_hyperp' `pi_hyperp' `betau_hyperp' ///
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
	
	frame post `memhold' (`dataset') ("`state'") (`mean_alpha2_xcon') (`sd_alpha2_xcon') (`luni_pi2') (`uuni_pi2') ///
												 (`a_pi2') (`b_pi2') (`mean_beta_ubin2') (`sd_beta_ubin2') ///
	                                             (`mean_alpha1_xcon') (`sd_alpha1_xcon') (`luni_pi1') (`uuni_pi1') ///
												 (`a_pi1') (`b_pi1')  (`mean_beta_ubin1') (`sd_beta_ubin1') ///
	                                             (`bapd_MCQBA_rc') `list_posts'
} // END OF QUIETLY STATEMENT	
} // END OF dataset FOR-LOOP
frame `memhold': save "`resultsfile'", replace		

* END OF DO-FILE