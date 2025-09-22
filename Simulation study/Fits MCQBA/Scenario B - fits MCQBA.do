/*************************************************************************************************************************************************
                       APPLIES A MONTE CARLO BIAS ANALYSIS TO DATA WITH A CONTINUOUS OUTCOME AND TWO CONTINUOUS UNMEASURED CONFOUNDERS

APPLICABLE WHEN:
	ANALYSIS OF INTEREST IS A LINEAR REGRESSION WHERE EXPOSURE X IS INCLUDED LINEARLY (E.G., NO ITERACTIONS WITH X)
	POINT MASS PRIORS OR PROBABILITY DISTRIBUTION PRIORS CAN BE PLACED ON THE BIAS PARAMETERS
		(SET SD OF BIAS PARAMETERS TO 0 FOR POINT MASS PRIORS)

ARGUMENTS: 
   exposurelist       - LIST OF NAMES FOR THE EXPOSURE VARIABLES (ORDERED AS `"x1, x2, ..., xp"') 
   confounderlist     - LIST OF NAMES FOR THE MEASURED CONFOUNDERS
   alpha2x_hyperp     - HYPERPARAMETERS FOR alpha2x (ORDERED AS mean_alpha2_x1, sd_alpha2_x1, ..., mean_alpha2_xp, sd_alpha2_xp) 
   eta2_hyperp        - HYPERPARAMETERS FOR eta2 AND 1/eta2sq ORDERED AS (`luni_eta2',`uuni_eta2',`shape_inveta2sq',`scale_inveta2sq') 
   alpha1x_hyperp     - HYPERPARAMETERS FOR alpha1x (ORDERED AS mean_alpha1_x1, sd_alpha1_x1, ..., mean_alpha1_xp, sd_alpha1_xp) 
   eta1_hyperp        - HYPERPARAMETERS FOR eta1 AND 1/eta1sq ORDERED AS (`luni_eta1',`uuni_eta1',`shape_inveta1sq',`scale_inveta1sq') 
   betau1_hyperp      - HYPERPARAMETERS FOR betau1 (ORDERED AS mean_beta_ucon1, sd_beta_ucon1) 
   betau2_hyperp      - HYPERPARAMETERS FOR betau2 (ORDERED AS mean_beta_ucon2, sd_beta_ucon2) 
   numMCreps          - NUMBER OF MONTE CARLO QBA REPLICATIONS
   namepostdata       - NAME OF FILE STORING MONTE CARLO FREQUENCY DISTRIBUTION OF BIAS-ADJUSTED EXPOSURE EFFECT ESTIMATES
   uniform            - 1 UNIFORM PRIOR; !1 GAMMA PRIOR

RETURNS: 
	DATASET CONTAINING RESULTS STORED BY POSTFILE; STORED UNDER FILENAME namepostdata

AUTHOR: RACHAEL HUGHES
*************************************************************************************************************************************************/
args exposurelist confounderlist alpha2x_hyperp eta2_hyperp alpha1x_hyperp eta1_hyperp betau1_hyperp betau2_hyperp numMCreps namepostdata uniform

noisily di "Exposure variable(s): `exposurelist'"
noisily di "Measured confounder variable(s): `confounderlist'"
noisily di "Conduct Monte Carlo QBA with `numMCreps' replications"
noisily di "Store temporary postfile results under filename `namepostdata'"

/***********************************************************
  DEFINE MACROS AND TEMPORARY NAMES, FILES, AND VARIABLES	
***********************************************************/
tempname postname

* TEMPORARY VARIABLES TO BE USED LATER IN THE DO-FILE
tempvar xb_u1tilde sim_u1tilde xb_u2tilde sim_u2tilde
gen double `sim_u1tilde' = .
gen double `xb_u1tilde' = .
gen double `sim_u2tilde' = .
gen double `xb_u2tilde' = .

local numexposures : list sizeof exposurelist
noisily di "Number of exposure variables is `numexposures'"

* EXTRACTING HYPERPARAMETER VALUES FOR PRIOR DISTRIBUTION OF alpha1x AND alpha2x
local column -1		// SET SUCH THAT -1+2=1
foreach exposure of local exposurelist {
	* UPDATE column FOR THE NEXT EXPOSURE
	local column = `column' + 2
	
	* EXTRACT HYPERPARAMETERS CORRESPONDING TO COEFFICIENTS OF exposurelist
	local mean_alpha2_`exposure' = `alpha2x_hyperp'[1,`column']
	local sd_alpha2_`exposure' = `alpha2x_hyperp'[1,`column'+1]
	noisily di "Draw alpha2_`exposure' from Normal distribution with mean `mean_alpha2_`exposure'' and standard deviation `sd_alpha2_`exposure''"
	
	local mean_alpha1_`exposure' = `alpha1x_hyperp'[1,`column']
	local sd_alpha1_`exposure' = `alpha1x_hyperp'[1,`column'+1]
	noisily di "Draw alpha1_`exposure' from Normal distribution with mean `mean_alpha1_`exposure'' and standard deviation `sd_alpha1_`exposure''"
}

* EXTRACTING AND DEFINING FOR PRIOR DISTRIBUTION FOR eta2 OR THEIR PRECISION
if `uniform'==1 {
	* UNIFORM DISTRIBUTION FOR eta2
	local luni_eta2 = `eta2_hyperp'[1,1]
	local uuni_eta2 = `eta2_hyperp'[1,2]
	local eta2prior `"scalar eta2star = runiform(`luni_eta2',`uuni_eta2')"'	
	noisily di "Draw eta2 from a Uniform distribution with lower limit `luni_eta2' and upper limit `uuni_eta2'"

	local luni_eta1 = `eta1_hyperp'[1,1]
	local uuni_eta1 = `eta1_hyperp'[1,2]
	local eta1prior `"scalar eta1star = runiform(`luni_eta1',`uuni_eta1')"'	
	noisily di "Draw eta1 from a Uniform distribution with lower limit `luni_eta1' and upper limit `uuni_eta1'"	
}
else {
	* GAMMA DISTRIBUTION FOR 1/eta2^2 AND 1/eta1^2
	local shape_inveta2sq = `eta2_hyperp'[1,3]
	local scale_inveta2sq = `eta2_hyperp'[1,4]
	local eta2prior `"scalar eta2star = sqrt(1/rgamma(`shape_inveta2sq',`scale_inveta2sq'))"'	
	noisily di "Draw 1/eta2sq from a Gamma distribution with shape `shape_inveta2sq' and scale `scale_inveta2sq'"

	local shape_inveta1sq = `eta1_hyperp'[1,3]
	local scale_inveta1sq = `eta1_hyperp'[1,4]
	local eta1prior `"scalar eta1star = sqrt(1/rgamma(`shape_inveta1sq',`scale_inveta1sq'))"'	
	noisily di "Draw 1/eta1sq from a Gamma distribution with shape `shape_inveta1sq' and scale `scale_inveta1sq'"
}
* EXTRACTING HYPERPARAMETER VALUES FOR PRIOR DISTRIBUTION OF betau2
local mean_beta_ucon2 = `betau2_hyperp'[1,1]
local sd_beta_ucon2 = `betau2_hyperp'[1,2]
noisily di "Draw beta_ucon2 from Normal distribution with mean `mean_beta_ucon2' and standard deviation `sd_beta_ucon2'"

* EXTRACTING HYPERPARAMETER VALUES FOR PRIOR DISTRIBUTION OF betau1
local mean_beta_ucon1 = `betau1_hyperp'[1,1]
local sd_beta_ucon1 = `betau1_hyperp'[1,2]
noisily di "Draw beta_ucon1 from Normal distribution with mean `mean_beta_ucon1' and standard deviation `sd_beta_ucon1'"

* LOCAL MACROLISTS FOR EXPOSURE ESTIMATES TO BE POSTED
	* CORRESPONDING TO COEFFICIENTS OF exposurelist
local list_roots `"betahat se_betahat betastar alpha2star alpha1star"' 	
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
local vars_otherbp `"eta2star eta1star betastar_ucon1 betastar_ucon2"'
local list_vars: list list_vars | vars_otherbp
local post_otherbp `"(eta2star) (eta1star) (betastar_ucon1) (betastar_ucon2)"'
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
	noisily di "Monte Carlo replication `t'"

	* INITIALISE
	quietly replace `xb_u1tilde' = 0
	quietly replace `xb_u2tilde' = 0
	
	/*****************************************
	   DRAW NEW SET OF BIAS PARAMETER VALUES
	*****************************************/
	* DRAW VALUE FOR eta1 AND eta2 (OR 1/eta1sq AND 1/eta2sq)
	`eta2prior'
	`eta1prior'
	
	* DRAW VALUE FOR betau1 AND betau2
	scalar betastar_ucon1 = rnormal(`mean_beta_ucon1',`sd_beta_ucon1')	
	scalar betastar_ucon2 = rnormal(`mean_beta_ucon2',`sd_beta_ucon2')	
	
	* DRAW VALUES FOR alpha1x AND alpha2x
	foreach exposure of local exposurelist {
		* DRAW alpha2star
		scalar alpha2star_`exposure' = rnormal(`mean_alpha2_`exposure'',`sd_alpha2_`exposure'')
		
		* DRAW alpha1star
		scalar alpha1star_`exposure' = rnormal(`mean_alpha1_`exposure'',`sd_alpha1_`exposure'')
				
		* SUM OVER EXPOSURE VARIABLE(S) TO GENERATE LINEAR PREDICTOR alpha2star*x
		quietly replace `xb_u2tilde' = `xb_u2tilde' + scalar(alpha2star_`exposure')*`exposure'
		
		* SUM OVER EXPOSURE VARIABLE(S) TO GENERATE LINEAR PREDICTOR alpha1star*x
		quietly replace `xb_u1tilde' = `xb_u1tilde' + scalar(alpha1star_`exposure')*`exposure'
	}	
	
	/******************************
	  SIMULATE u1tilde AND u2tilde
	*******************************/
	* SIMULATE u1tilde FROM CONDITIONAL DISTRIBUTION u1tilde | x, 0*c
	quietly replace `sim_u2tilde' = `xb_u2tilde' + scalar(eta2star)*rnormal() 		

	* SIMULATE u2tilde FROM CONDITIONAL DISTRIBUTION u2tilde | x, 0*c
	quietly replace `sim_u1tilde' = `xb_u1tilde' + scalar(eta1star)*rnormal() 	
	
	/**************************************************************************************************************
	  FIT ANALYSIS OF INTEREST y|x,c,u1tilde, u2tilde - CONSTRAINED REGRESSION
	  CONSTRAIN COEFFICIENT OF u1tilde AND u2tilde TO DRAWN VALUES betastar_ucon1 AND betastar_ucon2, RESPECTIVELY
	***************************************************************************************************************/
	* SET UP CONSTRAINT
	capture constraint drop _all
	constraint 1 `sim_u1tilde' = scalar(betastar_ucon1) 
	constraint 2 `sim_u2tilde' = scalar(betastar_ucon2) 
	
	* FIT ANALYSIS OF INTEREST
	capture cnsreg ycon `exposurelist' `confounderlist' `sim_u1tilde' `sim_u2tilde', constraints(1 2) 
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

* END OF DO-FILE