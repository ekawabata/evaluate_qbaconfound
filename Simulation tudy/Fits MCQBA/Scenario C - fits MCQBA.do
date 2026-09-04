/***********************************************************************************************************************************
                       APPLIES A MONTE CARLO BIAS ANALYSIS TO DATA WITH A BINARY OUTCOME AND BINARY UNMEASURED CONFOUNDER

APPLICABLE WHEN:
	ANALYSIS OF INTEREST IS A LOGISTIC REGRESSION WHERE EXPOSURE X IS INCLUDED LINEARLY (E.G., NO ITERACTIONS WITH X)
	POINT MASS PRIORS OR PROBABILITY DISTRIBUTION PRIORS CAN BE PLACED ON THE BIAS PARAMETERS
		(SET SD OF BIAS PARAMETERS TO 0 FOR POINT MASS PRIORS)

ARGUMENTS: 
   exposurelist       - LIST OF NAMES FOR THE EXPOSURE VARIABLES (ORDERED AS `"x1, x2, ..., xp"')
   confounderlist     - LIST OF NAMES FOR THE MEASURED CONFOUNDERS
   alphax_hyperp      - HYPERPARAMETERS FOR alphax (ORDERED AS mean_alpha_x1, sd_alphax1, ..., mean_alpha_xp, sd_alpha_xp) 
   pi_hyperp          - HYPERPARAMETERS FOR pi ORDERED AS (`luni_pi',`uuni_pi',`a_pi',`b_pi') 
   betau_hyperp       - HYPERPARAMETERS FOR betau (ORDERED AS mean_betau, sd_betau) 
   numMCreps          - NUMBER OF MONTE CARLO QBA REPLICATIONS
   namepostdata       - NAME OF FILE STORING MONTE CARLO FREQUENCY DISTRIBUTION OF BIAS-ADJUSTED EXPOSURE EFFECT ESTIMATES
   uniform            - 1 UNIFORM PRIOR; !1 BETA PRIOR

RETURNS: 
	DATASET CONTAINING RESULTS STORED BY POSTFILE; STORED UNDER FILENAME namepostdata
***********************************************************************************************************************************/
args exposurelist confounderlist alphax_hyperp pi_hyperp betau_hyperp numMCreps namepostdata uniform

noisily di "Exposure variable(s): `exposurelist'"
noisily di "Measured confounder variable(s): `confounderlist'"
noisily di "Conduct Monte Carlo QBA with `numMCreps' replications"
noisily di "Store temporary postfile results under filename `namepostdata'"

/************************************************************
  DEFINE MACROS AND TEMPORARY NAMES, FILES, AND VARIABLES	
***********************************************************/
tempname postname

* TEMPORARY VARIABLES TO BE USED LATER IN THE DO-FILE
tempvar xb_utilde sim_utilde pr_sim_utilde
gen double `sim_utilde' = .
gen double `xb_utilde' = .
gen double `pr_sim_utilde' = .

local numexposures : list sizeof exposurelist
noisily di "Number of exposure variables is `numexposures'"

* EXTRACTING HYPERPARAMETER VALUES FOR PRIOR DISTRIBUTION OF alphax
local column -1		// SET SUCH THAT -1+2=1
foreach exposure of local exposurelist {
	* UPDATE column FOR THE NEXT EXPOSURE
	local column = `column' + 2

	* EXTRACT HYPERPARAMETERS CORRESPONDING TO COEFFICIENTS OF exposurelist
	local mean_alpha_`exposure' = `alphax_hyperp'[1,`column']
	local sd_alpha_`exposure' = `alphax_hyperp'[1,`column'+1]
	noisily di "Draw alpha_`exposure' from Normal distribution with mean `mean_alpha_`exposure'' and standard deviation `sd_alpha_`exposure''"
}

* EXTRACTING AND DEFINING FOR PRIOR DISTRIBUTION FOR pi
if `uniform'==1 {
	* UNIFORM DISTRIBUTION FOR pi
	local luni_pi = `pi_hyperp'[1,1]
	local uuni_pi = `pi_hyperp'[1,2]
	local piprior `"scalar pistar = runiform(`luni_pi',`uuni_pi')"'	
	noisily di "Draw pi from a Uniform distribution with lower limit `luni_pi' and upper limit `uuni_pi'"
}
else {
	* BETA DISTRIBUTION FOR pi
	local a_pi = `pi_hyperp'[1,3]
	local b_pi = `pi_hyperp'[1,4]
	local piprior `"scalar pistar = rbeta(`a_pi',`b_pi')"'	
	noisily di "Draw pi from a Beta distribution with a=`a_pi' and b=`b_pi'"
}

* EXTRACTING HYPERPARAMETER VALUES FOR PRIOR DISTRIBUTION OF betau
local mean_beta_ubin = `betau_hyperp'[1,1]
local sd_beta_ubin = `betau_hyperp'[1,2]
noisily di "Draw beta_ubin from Normal distribution with mean `mean_beta_ubin' and standard deviation `sd_beta_ubin'"

* LOCAL MACROLISTS FOR EXPOSURE ESTIMATES TO BE POSTED
local list_roots `"betahat se_betahat betastar alphastar"' 	
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
local vars_otherbp `"pistar betastar_ubin"'
local list_vars: list list_vars | vars_otherbp
local post_otherbp `"(pistar) (betastar_ubin)"'
local list_posts: list list_posts | post_otherbp
* CHECK
*noisily di "`list_vars'"
*noisily di "`list_posts'"

/************************************************
  APPLY MONTE CARLO QUANTITATIVE BIAS ANALYSIS	
***********************************************/
* SET UP DATASET TO STORE RESULTS MONTE CARLO REPLICATIONS 
postfile `postname' MCQBArep analysis_rc analysis_eN `list_vars' using `namepostdata', replace

forvalues t=1(1)`numMCreps' {
	noisily di "Monte Carlo replication `t'"
	
	* INITIALISE
	quietly replace `xb_utilde' = 0
	
	/*****************************************
	   DRAW NEW SET OF BIAS PARAMETER VALUES
	*****************************************/
	* DRAW VALUE FOR pi
	`piprior'
	
	* DRAW VALUE FOR betau
	scalar betastar_ubin = rnormal(`mean_beta_ubin',`sd_beta_ubin')	
	
	* DRAW VALUES FOR alphax
	foreach exposure of local exposurelist {
		* DRAW alphastar
		scalar alphastar_`exposure' = rnormal(`mean_alpha_`exposure'',`sd_alpha_`exposure'')
		
		* SUM OVER EXPOSURE VARIABLE(S) TO GENERATE LINEAR PREDICTOR alphastar*x
		quietly replace `xb_utilde' = `xb_utilde' + scalar(alphastar_`exposure')*`exposure'
	}
	
	/****************************************************************************************************************************************
	 * DERIVE THE INTERCEPT OF logit{Pr(U`j'=1|X)} = alpha0 + alphax1*x1 + ... + alphaxp*xp BASED ON pistar AND alphax1*x1 + ... + alphaxp*xp
		* INITIAL VALUES FOR INTERVAL SEARCH OF alphacons (I.E., alpha0); 10 UNITS AROUND pistar ON THE LOG-ODDS SCALE
	*****************************************************************************************************************************************/
	local interval 10
	local logodds = ln(scalar(pistar)/(1-scalar(pistar)))
	local lower = `logodds' - `interval'/2
	local upper = `logodds' + `interval'/2

	* SEARCH FOR AN ESTIMATE OF alpha0 THAT RESULTS IN PR(sim_utilde=1) WITHIN 1 DECIMAL PLACE OF PROBABILITY pistar
		* SET TO A RELATIVELY LARGE VALUE TO ENTER while LOOP
	local smallabsdiff 5
	local count_while 0
	while `smallabsdiff'>0.001 & `count_while'<10 {
		local ++count_while
		
		* SET TO AN IMPOSSIBLY LARGE VALUE FOR FIRST CYCLE OF FOR-LOOP BELOW
		local absdiff_nearest 500
		
		* SIMULATE sim_utilde FOR DIFFERENT VALUES OF alpha0 (I.E., alphacons) WITHIN GIVEN INTERVAL
				* UPDATE INTERVAL TO SEARCH WITHIN [lower,upper]; ALSO USED TO CALCULATE LOWER AND UPPER LIMITS FOR THE NEXT CYCLE BELOW
		local interval = `interval'/10
		forvalues alphacons=`lower'(`interval')`upper' {
			quietly replace `pr_sim_utilde' = invlogit(`alphacons' + `xb_utilde')
			quietly replace `sim_utilde' = rbinomial(1,`pr_sim_utilde')
			
			* CALCULATE PR(sim_utilde=1) 
			quietly count if `sim_utilde' !=.
			local denom = r(N)
			quietly count if `sim_utilde'==1
			local current = r(N)/`denom'
			
			* CALCUATE ABSOLUTE DIFFERENCE FROM pistar
			local absdiff_current = abs(`current'-scalar(pistar))

			* UPDATE IF THE CURRENT VALUE OF alpha0 GIVES A PREVALENCE CLOSER TO pistar THAN AVAILABLE BEST ESTIMATE
			if `absdiff_current'<`absdiff_nearest' {
				local alphacons_nearest `alphacons'	
				local absdiff_nearest `absdiff_current'
			}		
		}	// END OF FORVALUES LOOP
		
		* UPDATE THE COUNTER OF THE while LOOP 
		local smallabsdiff `absdiff_nearest'
		
		* UPDATE THE LOWER AND UPPER LIMITS FOR THE NEXT CYCLE OF THE WHILE LOOP
		local lower = `alphacons_nearest' - `interval'/2
		local upper = `alphacons_nearest' + `interval'/2
	}	// END OF THE WHILTE LOOP
	
	* SET alpha0 TO THE CLOSEST VALUE FROM THE WHILE LOOP
	local alpha0 `alphacons_nearest'

	/*******************
	   SIMULATE utilde
	*******************/
	quietly replace `pr_sim_utilde' = invlogit(`alpha0' + `xb_utilde')
	quietly replace `sim_utilde' = rbinomial(1,`pr_sim_utilde')
	
	/*****************************************************************
	  FIT ANALYSIS OF INTEREST y|x,c,utilde - CONSTRAINED REGRESSION
	  CONSTRAIN COEFFICIENT OF utilde TO DRAWN VALUE betastar_ubin
	*****************************************************************/
	* SET UP CONSTRAINT
	capture constraint drop _all
	constraint 1 `sim_utilde' = scalar(betastar_ubin) 
	
	* FIT ANALYSIS OF INTEREST
	capture logit ybin `exposurelist' `confounderlist' `sim_utilde', constraints(1) 
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
postclose `postname'			

* END OF DO-FILE