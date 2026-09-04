/***********************************************************************************************************************************
                       APPLIES A MONTE CARLO BIAS ANALYSIS TO DATA WITH A NOMINAL OUTCOME AND CONTINUOUS UNMEASURED CONFOUNDER

APPLICABLE WHEN:
	ANALYSIS OF INTEREST IS A MULTINOMIAL LOGISTIC REGRESSION WHERE EXPOSURE X IS INCLUDED LINEARLY (E.G., NO ITERACTIONS WITH X)
	POINT MASS PRIORS OR PROBABILITY DISTRIBUTION PRIORS CAN BE PLACED ON THE BIAS PARAMETERS
		(SET SD OF BIAS PARAMETERS TO 0 FOR POINT MASS PRIORS)

ARGUMENTS: 
   exposurelist       - LIST OF NAMES FOR THE EXPOSURE VARIABLES (ORDERED AS `"x1, x2, ..., xp"')
   confounderlist     - LIST OF NAMES FOR THE MEASURED CONFOUNDERS
   alphax_hyperp      - HYPERPARAMETERS FOR alphax (ORDERED AS mean_alpha_x1, sd_alphax1, ..., mean_alpha_xp, sd_alpha_xp) 
   eta_hyperp         - HYPERPARAMETERS FOR eta AND 1/etasq ORDERED AS (`luni_eta',`uuni_eta',`shape_invetasq',`scale_invetasq') 
   betau_hyperp       - HYPERPARAMETERS FOR betau (ORDERED AS mean_betau, sd_betau) 
   numMCreps          - NUMBER OF MONTE CARLO QBA REPLICATIONS
   namepostdata       - NAME OF FILE STORING MONTE CARLO FREQUENCY DISTRIBUTION OF BIAS-ADJUSTED EXPOSURE EFFECT ESTIMATES
   uniform            - 1 UNIFORM PRIOR; !1 GAMMA PRIOR

RETURNS: 
	DATASET CONTAINING RESULTS STORED BY POSTFILE; STORED UNDER FILENAME namepostdata
***********************************************************************************************************************************/
args exposurelist confounderlist alphax_hyperp eta_hyperp betau_hyperp numMCreps namepostdata uniform
	 
noisily di "Exposure variable(s): `exposurelist'"
noisily di "Measured confounder variable(s): `confounderlist'"
noisily di "Conduct Monte Carlo QBA with `numMCreps' replications"
noisily di "Store temporary postfile results under filename `namepostdata'"

/************************************************************
  DEFINE MACROS AND TEMPORARY NAMES, FILES, AND VARIABLES	
***********************************************************/
tempname postname

* TEMPORARY VARIABLES TO BE USED LATER IN THE DO-FILE
tempvar xb_utilde sim_utilde 
gen double `sim_utilde' = .
gen double `xb_utilde' = .

local numexposures : list sizeof exposurelist
noisily di "Number of exposure variables is `numexposures'"


* EXTRACTING HYPERPARAMETER VALUES FOR PRIOR DISTRIBUTION OF alphax
local row 1
foreach exposure of local exposurelist {
	local mean_alpha_`exposure' = `alphax_hyperp'[`row',1]
	local sd_alpha_`exposure' = `alphax_hyperp'[`row',2]
	noisily di "Draw alpha_`exposure' from Normal distribution with mean `mean_alpha_`exposure'' and standard deviation `sd_alpha_`exposure''"

	* UPDATE row INDICATOR
	local ++row
} // END OF exposure FOREACH LOOP	

* EXTRACTING AND DEFINING PRIOR DISTRIBUTION FOR eta
if `uniform'==1 {
	* UNIFORM DISTRIBUTION FOR eta
	local luni_eta = `eta_hyperp'[1,1]
	local uuni_eta = `eta_hyperp'[1,2]
	local etaprior `"scalar etastar = runiform(`luni_eta',`uuni_eta')"'	
	noisily di "Draw eta from a Uniform distribution with lower limit `luni_eta' and upper limit `uuni_eta'"
}
else {
	* GAMMA DISTRIBUTION FOR 1/eta^2
	local shape_invetasq = `eta_hyperp'[1,3]
	local scale_invetasq = `eta_hyperp'[1,4]
	local etaprior `"scalar etastar = sqrt(1/rgamma(`shape_invetasq',`scale_invetasq'))"'	
	noisily di "Draw 1/etasq from a Gamma distribution with shape `shape_invetasq' and scale `scale_invetasq'"
}
		
* EXTRACTING HYPERPARAMETER VALUES FOR PRIOR DISTRIBUTION OF beta0ucon AND beta1ucon
forvalues k=0(1)1 {
	local mean_beta`k'_ucon = `betau_hyperp'[`k'+1,1]
	local sd_beta`k'_ucon = `betau_hyperp'[`k'+1,2]
	noisily di "Draw beta`k'_ucon from Normal distribution with mean `mean_beta`k'_ucon' and standard deviation `sd_beta`k'_ucon'"
}

* LOCAL MACROLISTS FOR EXPOSURE ESTIMATES TO BE POSTED
local list_roots `"beta0hat se_beta0hat beta0star beta1hat se_beta1hat beta1star alphastar"' 	
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
local vars_otherbp `"etastar beta0star_ucon beta1star_ucon"'
local list_vars: list list_vars | vars_otherbp
local post_otherbp `"(etastar) (beta0star_ucon) (beta1star_ucon)"'
local list_posts: list list_posts | post_otherbp
* CHECK
*noisily di "`list_vars'"
*noisily di "`list_posts'"

/************************************************
  APPLY MONTE CARLO QUANTITATIVE BIAS ANALYSIS	
***********************************************/
* SET UP DATASET TO STORE MONTE CARLO REPLICATIONS 
postfile `postname' MCQBArep analysis_rc analysis_eN `list_vars' using `namepostdata', replace

forvalues t=1(1)`numMCreps' {
	noisily di "Monte Carlo step `t'"
	
	* INITIALISE 
	quietly replace `xb_utilde' = 0 

	/*****************************************
	  DRAW NEW SET OF BIAS PARAMETER VALUES
	*****************************************/
	* DRAW VALUES FOR alphax
	foreach exposure of local exposurelist {
		* DRAW VALUES OF BIAS PARAMETERS alphastar_`exposure' 
		scalar alphastar_`exposure' = rnormal(`mean_alpha_`exposure'',`sd_alpha_`exposure'')
			
		* SUM OVER EXPOSURE VARIABLE(S) TO GENERATE LINEAR PREDICTOR alphastar*x 
		quietly replace `xb_utilde' = `xb_utilde' + scalar(alphastar_`exposure')*`exposure'		
	}
	
	* DRAW VALUE FOR eta OR 1/etasq
	`etaprior'
	
	* DRAW VALUES FOR BIAS PARAMETERS beta0ucon AND beta1ucon
	scalar beta0star_ucon = rnormal(`mean_beta0_ucon',`sd_beta0_ucon')	
	scalar beta1star_ucon = rnormal(`mean_beta1_ucon',`sd_beta1_ucon')	
		
	/*******************
	   SIMULATE utilde
	*******************/
	quietly replace `sim_utilde' = `xb_utilde' + scalar(etastar)*rnormal() 		

	/*************************************************************************************
	  FIT ANALYSIS OF INTEREST y|x,c,utilde - CONSTRAINED REGRESSION
	  CONSTRAIN COEFFICIENT OF utilde TO DRAWN VALUES beta0star_ucon AND beta1star_ucon
	*************************************************************************************/
	* SET UP CONSTRAINTS
	capture constraint drop _all
	constraint 1 [0]`sim_utilde' = scalar(beta0star_ucon) 
	constraint 2 [1]`sim_utilde' = scalar(beta1star_ucon) 
	
	* FIT ANALYSIS OF INTEREST
	capture noisily mlogit ynom `exposurelist' `confounderlist' `sim_utilde', constraints(1 2) baseoutcome(2)
	local analysis_rc = _rc

	if `analysis_rc'==0  {				// SUCCESSFULLY FITTED ANALYSIS MODEL
		* EXTRACT ANALYSIS RESULTS
		noisily di "Extracting results from the analysis model"
		local analysis_eN = e(N)		
		foreach exposure of local exposurelist {
			forvalues k=0(1)1 {
				scalar beta`k'hat_`exposure' = _b[`k':`exposure'] 
				scalar se_beta`k'hat_`exposure' = _se[`k':`exposure'] 
			
				* INCORPORATE RANDOM SAMPLING ERROR
				scalar beta`k'star_`exposure' = rnormal(scalar(beta`k'hat_`exposure'), scalar(se_beta`k'hat_`exposure'))					
			}
		}
	}
	else {								// PROBLEM ENCOUNTERED WHEN FITTING ANALYSIS MODEL	
		noisily di "Analysis model failed to fit"
		local analysis_eN . 		
		foreach exposure of local exposurelist {
			forvalues k=0(1)1 {
				scalar beta`k'hat_`exposure' = .
				scalar se_beta`k'hat_`exposure' = .
				scalar beta`k'star_`exposure' = .
			}
		}
	}
	
	post `postname' (`t') (`analysis_rc') (`analysis_eN') `list_posts' 	
	
} // END OF t FOR-LOOP
postclose `postname'			

* END OF DO-FILE