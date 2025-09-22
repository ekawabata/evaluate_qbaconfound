/***********************************************************************************************************************************
                       APPLIES A MONTE CARLO BIAS ANALYSIS TO DATA WITH A CONTINUOUS OUTCOME AND BINARY UNMEASURED CONFOUNDERS

APPLICABLE WHEN:
	ANALYSIS OF INTEREST IS A LOGISTIC REGRESSION WHERE EXPOSURE X IS INCLUDED LINEARLY (E.G., NO ITERACTIONS WITH X)
	POINT MASS PRIORS OR PROBABILITY DISTRIBUTION PRIORS CAN BE PLACED ON THE BIAS PARAMETERS
		(SET SD OF BIAS PARAMETERS TO 0 FOR POINT MASS PRIORS)

ARGUMENTS: 
   exposurelist       - LIST OF NAMES FOR THE EXPOSURE VARIABLES (ORDERED AS `"x1, x2, ..., xp"')
   confounderlist     - LIST OF NAMES FOR THE MEASURED CONFOUNDERS
   alphax_hyperp      - HYPERPARAMETERS FOR alphax (ORDERED AS mean_alpha_x1, sd_alphax1, ..., mean_alpha_xp, sd_alpha_xp) 
   pi_hyperp          - HYPERPARAMETERS FOR pi ORDERED AS (`luni_pi1',`uuni_pi1',`a_pi1',`b_pi1'
                                                           `luni_pi2',`uuni_pi2',`a_pi2',`b_pi2') 
   betau_hyperp       - HYPERPARAMETERS FOR betau (ORDERED AS mean_betau1, sd_betau1
                                                              mean_betau2, sd_betau2) 
   numMCreps          - NUMBER OF MONTE CARLO QBA REPLICATIONS
   namepostdata       - NAME OF FILE STORING MONTE CARLO FREQUENCY DISTRIBUTION OF BIAS-ADJUSTED EXPOSURE EFFECT ESTIMATES
   uniform            - 1 UNIFORM PRIOR; !1 BETA PRIOR

RETURNS: 
	DATASET CONTAINING RESULTS STORED BY POSTFILE; STORED UNDER FILENAME namepostdata
*********************************************************************************************************************************/
args exposurelist confounderlist alphax_hyperp pi_hyperp betau_hyperp numMCreps namepostdata uniform

noisily di "Exposure variable(s): `exposurelist'"
noisily di "Measured confounder variable(s): `confounderlist'"
noisily di "Conduct Monte Carlo QBA with `numMCreps' replications"
noisily di "Store temporary postfile results under filename `namepostdata'"

/************************************************************
  DEFINE MACROS AND TEMPORARY NAMES, FILES, AND VARIABLES	
***********************************************************/
tempname postname
tempvar xb_utilde1 sim_utilde1 pr_sim_utilde1 xb_utilde2 sim_utilde2 pr_sim_utilde2

local numexposures : list sizeof exposurelist
noisily di "Number of exposure variables is `numexposures'"

forvalues j=1(1)2 {
	* TEMPORARY VARIABLES TO BE USED LATER IN THE DO-FILE
	gen double `sim_utilde`j'' = .
	gen double `xb_utilde`j'' = .
	gen double `pr_sim_utilde`j'' = .
	
	* EXTRACTING HYPERPARAMETER VALUES FOR PRIOR DISTRIBUTION OF alphax`j'
	local column -1		// SET SUCH THAT -1+2=1
	foreach exposure of local exposurelist {
		* UPDATE column FOR THE NEXT EXPOSURE
		local column = `column' + 2

		* EXTRACT HYPERPARAMETERS CORRESPONDING TO COEFFICIENTS OF exposurelist
		local mean_alpha`j'_`exposure' = `alphax_hyperp'[`j',`column']
		local sd_alpha`j'_`exposure' = `alphax_hyperp'[`j',`column'+1]
		noisily di "Draw alpha`j'_`exposure' from Normal distribution with mean `mean_alpha`j'_`exposure'' and standard deviation `sd_alpha`j'_`exposure''"
	}
	
	* EXTRACTING AND DEFINING FOR PRIOR DISTRIBUTION FOR pi`j'
	if `uniform'==1 {
		* UNIFORM DISTRIBUTION FOR pi`j'
		local luni_pi`j' = `pi_hyperp'[`j',1]
		local uuni_pi`j' = `pi_hyperp'[`j',2]
		local pi`j'prior `"scalar pi`j'star = runiform(`luni_pi`j'',`uuni_pi`j'')"'	
		noisily di "Draw pi`j' from a Uniform distribution with lower limit `luni_pi`j'' and upper limit `uuni_pi`j''"
	}
	else {
		* BETA DISTRIBUTION FOR pi`j'
		local a_pi`j' = `pi_hyperp'[`j',3]
		local b_pi`j' = `pi_hyperp'[`j',4]
		local pi`j'prior `"scalar pi`j'star = rbeta(`a_pi`j'',`b_pi`j'')"'	
		noisily di "Draw pi`j' from a Beta distribution with a=`a_pi`j'' and b=`b_pi`j''"
	}
	
	* EXTRACTING HYPERPARAMETER VALUES FOR PRIOR DISTRIBUTION OF betau`j'
	local mean_beta_ubin`j' = `betau_hyperp'[`j',1]
	local sd_beta_ubin`j' = `betau_hyperp'[`j',2]
	noisily di "Draw beta_ubin`j' from Normal distribution with mean `mean_beta_ubin`j'' and standard deviation `sd_beta_ubin`j''"
}

* LOCAL MACROLISTS FOR EXPOSURE ESTIMATES TO BE POSTED
local list_roots `"betahat se_betahat betastar alpha1star alpha2star"' 	
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
local vars_otherbp `"pi1star betastar_ubin1 pi2star betastar_ubin2"'
local list_vars: list list_vars | vars_otherbp
local post_otherbp `"(pi1star) (betastar_ubin1) (pi2star) (betastar_ubin2)"'
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
	
	forvalues j=1(1)2 {
		* INITIALISE 
		quietly replace `xb_utilde`j'' = 0
		
	   /***************************************
	    DRAW NEW SET OF BIAS PARAMETER VALUES
	   ***************************************/
		* DRAW VALUE FOR pi`j'
		`pi`j'prior'
		
		* DRAW VALUE FOR betau`j'
		scalar betastar_ubin`j' = rnormal(`mean_beta_ubin`j'',`sd_beta_ubin`j'')	
		
		* DRAW VALUES FOR alphax`j'
		foreach exposure of local exposurelist {
			* DRAW alphastar
			scalar alpha`j'star_`exposure' = rnormal(`mean_alpha`j'_`exposure'',`sd_alpha`j'_`exposure'')
			
			* SUM OVER EXPOSURE VARIABLE(S) TO GENERATE LINEAR PREDICTOR alpha`j'star*x
			quietly replace `xb_utilde`j'' = `xb_utilde`j'' + scalar(alpha`j'star_`exposure')*`exposure'
		}
		
		/****************************************************************************************************************************************
		DERIVE THE INTERCEPT OF logit{Pr(U`j'=1|X)} = alpha`j'0 + alpha`j'star_x*x 
			INITIAL VALUES FOR INTERVAL SEARCH OF alphacons (I.E., alpha`j'0); 10 UNITS AROUND pi`j'star ON THE LOG-ODDS SCALE
		*****************************************************************************************************************************************/
		local interval 10
		local logodds = ln(scalar(pi`j'star)/(1-scalar(pi`j'star)))
		local lower = `logodds' - `interval'/2
		local upper = `logodds' + `interval'/2

		* SEARCH FOR AN ESTIMATE OF alphacons THAT RESULTS IN A PR(sim_utilde`j'=1) WITHIN 1 DECIMAL PLACE OF pistar`j'
			* SET TO A RELATIVELY LARGE VALUE TO ENTER while LOOP
		local smallabsdiff 500
		local count_while 0
		while `smallabsdiff'>0.001 & `count_while'<10 {
			local ++count_while
			
			* SET TO A IMPOSSIBLE LARGE VALUE FOR FIRST CYCLE OF FOR-LOOP
			local absdiff_nearest 500
			
			* SIMULATE sim_utilde`j' FOR DIFFERENT VALUES OF alphacons WITHIN GIVEN INTERVAL
					* UPDATE INTERVAL TO SEARCH WITHIN [lower,upper]; ALSO USED TO CALCULATE LOWER AND UPPER LIMITS FOR THE NEXT CYCLE BELOW
			local interval = `interval'/10
			forvalues alphacons=`lower'(`interval')`upper' {
				quietly replace `pr_sim_utilde`j'' = invlogit(`alphacons' + `xb_utilde`j'')
				quietly replace `sim_utilde`j'' = rbinomial(1,`pr_sim_utilde`j'')
				
				* CALCULATE PR(sim_utilde`j'=1) 
				quietly count if `sim_utilde`j'' !=.
				local denom = r(N)
				quietly count if `sim_utilde`j''==1
				local current = r(N)/`denom'
				
				* CALCUATE ABSOLUTE DIFFERENCE FROM pistar`j'
				local absdiff_current = abs(`current'-scalar(pi`j'star))

				* UPDATE IF THE CURRENT VALUE OF alpha`j'0 GIVES A PREVALENCE CLOSER TO pi`j'star THAN AVAILABLE BEST ESTIMATE
				if `absdiff_current'<`absdiff_nearest' {
					local alphacons_nearest `alphacons'	
					local absdiff_nearest `absdiff_current'
				}		
			}	// END OF FORVALUES alphacons LOOP
						
			* UPDATE THE COUNTER OF THE while LOOP 
			local smallabsdiff `absdiff_nearest'
			
			* UPDATE THE LOWER AND UPPER LIMITS FOR THE NEXT CYCLE OF THE WHILE LOOP
			local lower = `alphacons_nearest' - `interval'/2
			local upper = `alphacons_nearest' + `interval'/2
		}	// END OF THE WHILE LOOP
	
		* SET alpha`j'0 TO THE CLOSEST VALUE FROM THE WHILE LOOP
		local alpha`j'0 `alphacons_nearest'
		
		/*******************
	     SIMULATE utilde`j'
	    *******************/
		quietly replace `pr_sim_utilde`j'' = invlogit(`alpha`j'0' + `xb_utilde`j'')
		quietly replace `sim_utilde`j'' = rbinomial(1,`pr_sim_utilde`j'')
	} // END OF FORVALUES j LOOP
	
	/**************************************************************************
	  FIT ANALYSIS OF INTEREST y|x,c,utilde1, utilde2 - CONSTRAINED REGRESSION
	  CONSTRAIN COEFFICIENT OF utilde TO DRAWN VALUE betastar_ubin
	***************************************************************************/
	* SET UP CONSTRAINTS
	capture constraint drop _all
	constraint 1 `sim_utilde1' = scalar(betastar_ubin1) 
	constraint 2 `sim_utilde2' = scalar(betastar_ubin2) 
	
	* FIT ANALYSIS OF INTEREST
	capture cnsreg ycon `exposurelist' `confounderlist' `sim_utilde1' `sim_utilde2', constraints(1 2) 
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