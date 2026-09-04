/*************************************************************************************************************************************************
                                                         APPLIES MONTE CARLO QBA TO NHANES EXAMPLE 

TWO UNMEASURED CONFOUNDERS: uconE (CONTINUOUS; energy) AND ubinM (BINARY; martialstatus)

PRIOR DISTRIBUTION FOR THE MARGINAL PRAVELNCE Pr(ubinM=1)=piM
													
RETURNS: 
	DATASET CONTAINING RESULTS STORED BY POSTFILE; STORED UNDER FILENAME OF MACRO resultsfile
*************************************************************************************************************************************************/
* SET WORKING DIRECTORY
* e.g., cd "C:\Github\Applied example\"

version 18

* DROP ALL SCALARS
scalar drop _all

/**********************************************
*       IMPORT NHANES DATASET

import delimited "Data\Data.csv", clear 

/**********************************
  CONVERT INTO NUMERICAL VARIABLES
***********************************/
* EXPOSURE OF INTEREST
encode alcoholx, gen(alcohol) 
tab alcoholx alcohol, m

* GENERATE DUMMY VARIABLES FOR THE EXPOSURE
xi i.alcohol
rename _Ialcohol_2 x3to4
rename _Ialcohol_3 xatleast5
rename _Ialcohol_4 xformer 
rename _Ialcohol_5 xnever

* MEASURED CONFOUNDERS
local confounderlist `"age race education physicalactivitylevel year recall1 recall2 smoking income misreporting"'
foreach item of local confounderlist {
	encode `item'c, gen(`item')	
	tab `item'c `item', m
}

* DROP UNWANTED VARIABLES
drop alcoholx-misreportingc
*********************************************/

/****************************************************
   MACROS AND TEMPORARY NAMES, FILES, AND VARIABLES
***************************************************/
* NUMBER OF MONTE CARLO REPLICATIONS
local numMCreps 10000

tempname postname
tempvar xb_utildeE sim_utildeE xb_utildeM pr_sim_utildeM sim_utildeM

* TEMPORARY VARIABLES TO BE USED LATER IN THE DO-FILE
gen double `sim_utildeE' = .
gen double `xb_utildeE' = .
gen double `sim_utildeM' = .
gen double `pr_sim_utildeM' = .
gen double `xb_utildeM' = .

* MACROLISTS
local exposurelist `"x3to4 xatleast5 xformer xnever"'
local confounderlist `"age race education physicalactivitylevel year recall1 recall2 smoking income misreporting"'

* HYPERPARAMETERS OF THE PRIOR DISTRIBUTIONS FOR THE BIAS PARAMETERS 
	* COEFFICIENTS OF EXPOSURE IN UNMEASURED CONFOUNDER MODEL FOR uconE
local mean_alphaE_x3to4 -65.37   
local mean_alphaE_xatleast5 -94.71   
local mean_alphaE_xformer 69.72   
local mean_alphaE_xnever -43.47  
local sd_alphaE_x3to4 18.07   
local sd_alphaE_xatleast5 19.64  
local sd_alphaE_xformer 18.63
local sd_alphaE_xnever 27.21
	* etaE=544.5
local luni_etaE 372.6
local uuni_etaE 716.5
	* COEFFICIENTS OF EXPOSURE IN UNMEASURED CONFOUNDER MODEL FOR ubinM
local mean_alphaM_x3to4 0.3831   
local mean_alphaM_xatleast5 0.3785   
local mean_alphaM_xformer -0.07963  
local mean_alphaM_xnever 0.4524   
local sd_alphaM_x3to4 0.07405  
local sd_alphaM_xatleast5 0.07958 
local sd_alphaM_xformer 0.08107
local sd_alphaM_xnever 0.1081
	* UNIFORM DISTRIBUTION FOR PREVALENCE OF ubinM
local luni_piM 0.2900
local uuni_piM 0.3484
	* COEFFICIENT OF uconE IN SUBSTANTIVE ANALYSIS MODEL
local mean_beta_uconE 0.006931 
local sd_beta_uconE 0.0002805
	* COEFFICIENT OF ubinM IN SUBSTANTIVE ANALYSIS MODEL
local mean_beta_ubinM -1.312 
local sd_beta_ubinM 0.3453

/*************************************
  GENERATE LISTS FOR POST STATEMENTS
*************************************/
* LOCAL MACROLISTS FOR EXPOSURE ESTIMATES TO BE POSTED
local list_roots `"betahat se_betahat betastar alphaEstar alphaMstar"' 	
local list_vars ""
local list_posts ""
foreach root of local list_roots {
	local post_`root' ""
	local vars_`root' ""
	foreach exposure of local exposurelist {
		local variable "`root'_`exposure'"
		local vars_`root': list vars_`root' | variable
		local posting `"(`root'_`exposure')"'
		local post_`root': list post_`root' | posting		
	}
	local list_vars: list list_vars | vars_`root'
	local list_posts: list list_posts | post_`root'
}
local vars_otherbp `"etaEstar piMstar betastar_uconE betastar_ubinM"'
local list_vars: list list_vars | vars_otherbp
local post_otherbp `"(etaEstar) (piMstar) (betastar_uconE) (betastar_ubinM)"'
local list_posts: list list_posts | post_otherbp
* CHECK
*noisily di "`list_vars'"
*noisily di "`list_posts'"
				
/************************************************************
   CODE TO APPLY THE MONTE CARLO QUANTITATIVE BIAS ANALYSIS
*************************************************************/
* SET UP DATASET TO STORE MONTE CARLO REPLICATIONS
local resultsfile `"Monte Carlo\mcqba_`numMCreps'reps_piM"'
postfile `postname' MCQBArep analysis_rc analysis_eN `list_vars' using "`resultsfile'", replace

set seed 1800

quietly {
forvalues t=1(1)`numMCreps' {
	noisily di "Monte Carlo replication `t'"
	
	* INITIALISE 
	quietly replace `xb_utildeE' = 0 
	quietly replace `xb_utildeM' = 0
	
	/***************************************
	  DRAW NEW SET OF BIAS PARAMETER VALUES
	***************************************/
	* DRAW VALUES FOR alphaE AND alphaM
	foreach exposure of local exposurelist {
		* DRAW VALUES OF BIAS PARAMETERS alpha1star_`exposure' AND alpha2star_`exposure'
		scalar alphaEstar_`exposure' = rnormal(`mean_alphaE_`exposure'',`sd_alphaE_`exposure'')
		scalar alphaMstar_`exposure' = rnormal(`mean_alphaM_`exposure'',`sd_alphaM_`exposure'')
		
		* SUM OVER EXPOSURE VARIABLE(S) TO GENERATE LINEAR PREDICTOR alphaEstar*x AND alphaMstar*x
		quietly replace `xb_utildeE' = `xb_utildeE' + scalar(alphaEstar_`exposure')*`exposure'
		quietly replace `xb_utildeM' = `xb_utildeM' + scalar(alphaMstar_`exposure')*`exposure'
	}

	* DRAW VALUE FOR BIAS PARAMETER etaE
	scalar etaEstar = runiform(`luni_etaE',`uuni_etaE')

	* DRAW VALUE FOR THE INTERCEPT OF THE UNMEASURED CONFOUNDER MODEL FOR ubinM
	scalar piMstar = runiform(`luni_piM', `uuni_piM')

	* DRAW VALUES FOR BIAS PARAMETERS betauconE AND betaubinM
	scalar betastar_uconE = rnormal(`mean_beta_uconE',`sd_beta_uconE')	
	scalar betastar_ubinM = rnormal(`mean_beta_ubinM',`sd_beta_ubinM')
	
	/************************************************
	  DERIVE THE INTERCEPT OF logit{Pr(UtildeM=1|X)}
	*************************************************/
	* INITIAL VALUES FOR INTERVAL SEARCH OF INTERCEPT; 10 UNITS AROUND pistar ON THE LOG-ODDS SCALE
	local interval 10
	local logodds = ln(scalar(piMstar)/(1-scalar(piMstar)))
	local lower = `logodds' - `interval'/2              
	local upper = `logodds' + `interval'/2

	* SEARCH FOR AN ESTIMATE OF alphacons THAT RESULTS IN A PR(sim_utildeM=1) WITHIN 1 DECIMAL PLACE OF piMstar
		* SET TO A LARGE VALUE TO ENTER while LOOP
	local smallabsdiff 5
	local count_while 0
	while `smallabsdiff'>0.0001 & `count_while'<10 {
		local ++count_while
			
		* SET TO A IMPOSSIBLE LARGE VALUE FOR FIRST CYCLE OF FOR-LOOP
		local absdiff_nearest 500
			
		* SIMULATE sim_utildeM FOR DIFFERENT VALUES OF INTERCEPT (i.e., alphaMcons) WITHIN GIVEN INTERVAL
				* UPDATE INTERVAL TO SEARCH WITHIN [lower,upper]; ALSO USED TO CALCULATE LOWER AND UPPER LIMITS FOR THE NEXT CYCLE
		local interval = `interval'/10
		forvalues alphaMcons=`lower'(`interval')`upper' {
			quietly replace `pr_sim_utildeM' = invlogit(`alphaMcons' + `xb_utildeM')
			quietly replace `sim_utildeM' = rbinomial(1,`pr_sim_utildeM')
				
			* CALCULATE PR(sim_utildeM) 
			quietly count if `sim_utildeM' !=.
			local denom = r(N)
			quietly count if `sim_utildeM'==1
			local current = r(N)/`denom'
				
			* CALCUATE ABSOLUTE DIFFERENCE FROM pistar
			local absdiff_current = abs(`current'-scalar(piMstar))
		
			if `absdiff_current'<`absdiff_nearest' {
				local alphaMcons_nearest `alphaMcons'	
				local absdiff_nearest `absdiff_current'
			}		
		}
		
		* UPDATE THE COUNTER OF THE while LOOP 
		local smallabsdiff `absdiff_nearest'
			
		* UPDATE THE LOWER AND UPPER LIMITS FOR THE NEXT CYCLE OF THE WHILE LOOP
		local lower = `alphaMcons_nearest' - `interval'/2
		local upper = `alphaMcons_nearest' + `interval'/2
	}
	
	* SET alphaM0 TO THE CLOSEST VALUE FROM THE WHILE LOOP
	local alphaM0 `alphaMcons_nearest' 

	/**********************************************
	  SIMULATE PROXIES FOR UNMEASURED CONFOUNDERS
	***********************************************/
	* SIMULATE sim_utildeE
	quietly replace `sim_utildeE' = `xb_utildeE' + scalar(etaEstar)*rnormal() 		

	* SIMULATE sim_utildeM 
	quietly replace `xb_utildeM' = `alphaM0' + `xb_utildeM'
	quietly replace `pr_sim_utildeM' = invlogit(`xb_utildeM')
	quietly replace `sim_utildeM' = rbinomial(1,`pr_sim_utildeM')
	
	/***************************************
	  FIT CONSTRAINED SUBSTANTIVE ANALYSIS 
	***************************************/
	* DEFINE CONSTRAINTS
	capture constraint drop _all
	constraint 1 `sim_utildeE' = scalar(betastar_uconE) 
	constraint 2 `sim_utildeM' = scalar(betastar_ubinM) 
	
	* FIT ANALYSIS MODEL
	capture cnsreg waist `exposurelist' i.(`confounderlist') `sim_utildeE' `sim_utildeM', constraints(1 2) 
	local analysis_rc = _rc

	if `analysis_rc'==0  {				// SUCCESSFULLY FITTED ANALYSIS MODEL
		* EXTRACT ANALYSIS RESULTS
		noisily di "Extracting results from the analysis model"
		local analysis_eN = e(N)		
		foreach exposure of local exposurelist {
			scalar betahat_`exposure' = _b[`exposure'] 
			scalar se_betahat_`exposure' = _se[`exposure'] 
			
			* INCORPORATE RANDOM SAMPLING ERROR
			scalar betastar_`exposure' = rnormal(scalar(betahat_`exposure'), scalar(se_betahat_`exposure'))	
		}
	}
	else {								// PROBLEM ENCOUNTERED WHEN FITTING ANALYSIS MODEL	
		noisily di "Analysis model failed to fit"
		local analysis_eN . 		
		foreach exposure of local exposurelist {
			scalar betahat_`exposure' = .
			scalar se_betahat_`exposure' = .
			scalar betastar_`exposure' = .
		}
	}

	post `postname' (`t') (`analysis_rc') (`analysis_eN') `list_posts' 
} // END OF t FOR-LOOP
postclose `postname'			// CLOSES postfile STATEMENT AND SAVES DATASET 
} // END OF QUIETLY

/****************************************************
   SUMMARIZE MONTE CARLO EMPIRICAL DISTRIBUTION
***************************************************/
* OPEN DATASET CONTAINING MONTE CARLO REPLICATIONS
use "`resultsfile'", clear

* NO ERRORS REPORTED WHEN FITTING THE ANALYSIS MODEL
tab analysis_rc, m

* ANALYSIS MODEL FITTED TO FULL SAMPLE EVERYTIME
tab analysis_eN, m

* MATRIX TO HOLD THE RESULTS FOR THE POINT AND INTERVAL ESTIMATES
matrix results = J(7,4,.)
matrix rownames results = "3to4:all" "." "5+:all" "." "form:all" "." "nev:all"
matrix colnames results = "Coef" "SD" "lint" "uint"

local exposurelist `"x3to4 xatleast5 xformer xnever"'

* SUMMARIZE FREQUENCY DISTRIBUTION OF MONTE CARLO ESTIMATES	
quietly {
local row 0
foreach exposure of local exposurelist {
	noisily di "Processing results for `exposure'"
	
	* RESULTS ACCOUNTING FOR UNCERTAINTY ABOUT SYSTEMATIC ERROR AND RANDOM SAMPLING
	local ++row
		* COEF
	centile betastar_`exposure', centile(50)
	matrix results[`row',1] = r(c_1)
		* SD
	summarize betastar_`exposure'
	matrix results[`row',2] = r(sd)
		* LOWER LIMIT OF INTERVAL
	centile betastar_`exposure', centile(2.5)
	matrix results[`row',3] = r(c_1)
		* UPPER LIMIT OF INTERVAL
	centile betastar_`exposure', centile(97.5)
	matrix results[`row',4] = r(c_1)
	
	* MISSING ROW
	local ++row
}	// END OF FOREACH STATEMENT
}	// END OF QUIETLY

* PRINT RESULTS TO THE SCREEN
matrix list results
	
* END OF DO-FILE