/************************************************************************************************************************************************
                   APPLIES A MONTE CARLO BIAS ANALYSIS TO DATA WITH A CONTINUOUS OUTCOME AND BINARY AND CONTINUOUS UNMEASURED CONFOUNDERS

APPLICABLE WHEN:
	ANALYSIS OF INTEREST IS A LINEAR REGRESSION WHERE EXPOSURE X IS INCLUDED LINEARLY (E.G., NO ITERACTIONS WITH X)
	POINT MASS PRIORS OR PROBABILITY DISTRIBUTION PRIORS CAN BE PLACED ON THE BIAS PARAMETERS
		(SET SD OF BIAS PARAMETERS TO 0 FOR POINT MASS PRIORS)

ARGUMENTS: 
   exposurelist       - LIST OF NAMES FOR THE EXPOSURE VARIABLES (ORDERED AS `"x1, x2, ..., xp"')
   confounderlist     - LIST OF NAMES FOR THE MEASURED CONFOUNDERS
   alpha1x_hyperp     - MATRIX OF HYPERPARAMETERS FOR alpha1x (mean_alpha1_x1, sd_alpha1_x1 \ 
                                                               ............................ \ 
															   mean_alpha1_xp, sd_alpha1_xp)
															   
   eta1_hyperp        - VECTOR OF HYPERPARAMETERS FOR eta1 (`luni_eta1',`uuni_eta1',`shape_inveta1sq',`scale_inveta1sq')
   alpha1x_hyperp     - MATRIX OF HYPERPARAMETERS FOR alpha2x (mean_alpha2_x1, sd_alpha2_x1 \ 
                                                               ............................ \ 
															   mean_alpha2_xp, sd_alpha2_xp)
															   
   pi_hyperp          - VECTOR OF HYPERPARAMETERS FOR pi2  (`luni_pi2',`uuni_pi2',`a_pi2',`b_pi2')
   betau1_hyperp      - VECTOR OF HYPERPARAMETERS FOR betau1 (mean_betau1, sd_betau1)
   betau2_hyperp      - VECTOR OF HYPERPARAMETERS FOR betau2 (mean_betau2, sd_betau2)
   numMCreps          - NUMBER OF MONTE CARLO QBA REPLICATIONS
   namepostdata       - NAME OF FILE STORING MONTE CARLO FREQUENCY DISTRIBUTION OF BIAS-ADJUSTED EXPOSURE EFFECT ESTIMATES
   uniform            - 1 UNIFORM PRIOR; !1 BETA PRIOR

RETURNS: 
	DATASET CONTAINING RESULTS STORED BY POSTFILE; STORED UNDER FILENAME namepostdata
*********************************************************************************************************************************/
args exposurelist confounderlist alphax_hyperp eta1_hyperp pi2_hyperp betau_hyperp numMCreps namepostdata uniform
	 
noisily di "Exposure variable(s): `exposurelist'"
noisily di "Measured confounder variable(s): `confounderlist'"
noisily di "Conduct Monte Carlo QBA with `numMCreps' replications"
noisily di "Store temporary postfile results under filename `namepostdata'"

/************************************************************
  DEFINE MACROS AND TEMPORARY NAMES, FILES, AND VARIABLES	
***********************************************************/
tempname postname
tempvar xb_utilde1 sim_utilde1 xb_utilde2 sim_utilde2 pr_sim_utilde2

local numexposures : list sizeof exposurelist
noisily di "Number of exposure variables is `numexposures'"

* TEMPORARY VARIABLES TO BE USED LATER IN THE DO-FILE
gen double `sim_utilde1' = .
gen double `xb_utilde1' = .
gen double `sim_utilde2' = .
gen double `xb_utilde2' = .
gen double `pr_sim_utilde2' = .

forvalues j=1(1)2 {
	* EXTRACTING HYPERPARAMETER VALUES FOR PRIOR DISTRIBUTION OF alpha`j'x
	local column 1
	foreach exposure of local exposurelist {
		local mean_alpha`j'_`exposure' = `alphax_hyperp'[`j',`column']
		local sd_alpha`j'_`exposure' = `alphax_hyperp'[`j',`column'+1]
		noisily di "Draw alpha`j'_`exposure' from Normal distribution with mean `mean_alpha`j'_`exposure'' and standard deviation `sd_alpha`j'_`exposure''"

		* UPDATE row INDICATOR
		local column = `column' + 2
	} // END OF exposure FOREACH LOOP	
} // END OF j FORVALUES LOOP			

* EXTRACTING AND DEFINING PRIOR DISTRIBUTION FOR eta1 and pi2
if `uniform'==1 {
	* UNIFORM DISTRIBUTION FOR eta1
	local luni_eta1 = `eta1_hyperp'[1,1]
	local uuni_eta1 = `eta1_hyperp'[1,2]
	local eta1prior `"scalar eta1star = runiform(`luni_eta1',`uuni_eta1')"'	
	noisily di "Draw eta1 from a Uniform distribution with lower limit `luni_eta1' and upper limit `uuni_eta1'"
	
	* UNIFORM DISTRIBUTION FOR pi2
	local luni_pi2 = `pi2_hyperp'[1,1]
	local uuni_pi2 = `pi2_hyperp'[1,2]
	local pi2prior `"scalar pi2star = runiform(`luni_pi2',`uuni_pi2')"'	
	noisily di "Draw pi2 from a Uniform distribution with lower limit `luni_pi2' and upper limit `uuni_pi2'"
}
else {
	* GAMMA DISTRIBUTION FOR 1/eta1^2
	local shape_inveta1sq = `eta1_hyperp'[1,3]
	local scale_inveta1sq = `eta1_hyperp'[1,4]
	local eta1prior `"scalar eta1star = sqrt(1/rgamma(`shape_inveta1sq',`scale_inveta1sq'))"'	
	noisily di "Draw 1/eta1sq from a Gamma distribution with shape `shape_inveta1sq' and scale `scale_inveta1sq'"

	* BETA DISTRIBUTION FOR pi
	local a_pi2 = `pi_hyperp'[`j',3]
	local b_pi2 = `pi_hyperp'[`j',4]
	local pi2prior `"scalar pi2star = rbeta(`a_pi2',`b_pi2')"'	
	noisily di "Draw pi2 from a Beta distribution with a=`a_pi2' and b=`b_pi2'"
}
	
* EXTRACTING HYPERPARAMETER VALUES FOR PRIOR DISTRIBUTION OF betaucon1
local mean_beta_ucon1 = `betau_hyperp'[1,1]
local sd_beta_ucon1 = `betau_hyperp'[1,2]
noisily di "Draw beta_ucon1 from Normal distribution with mean `mean_beta_ucon1' and standard deviation `sd_beta_ucon1'"
	
* EXTRACTING HYPERPARAMETER VALUES FOR PRIOR DISTRIBUTION OF betaubin2
local mean_beta_ubin2 = `betau_hyperp'[2,1]
local sd_beta_ubin2 = `betau_hyperp'[2,2]
noisily di "Draw beta_ubin2 from Normal distribution with mean `mean_beta_ubin2' and standard deviation `sd_beta_ubin2'"
		
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
local vars_otherbp `"eta1star betastar_ucon1 pi2star betastar_ubin2"'
local list_vars: list list_vars | vars_otherbp
local post_otherbp `"(eta1star) (betastar_ucon1) (pi2star) (betastar_ubin2)"'
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
	quietly replace `xb_utilde1' = 0 
	quietly replace `xb_utilde2' = 0

	/***************************************
	  DRAW NEW SET OF BIAS PARAMETER VALUES
	***************************************/
	* DRAW VALUES FOR alphax1 AND alphax2
	foreach exposure of local exposurelist {
		* DRAW VALUES OF BIAS PARAMETERS alpha1star_`exposure' AND alpha2star_`exposure'
		scalar alpha1star_`exposure' = rnormal(`mean_alpha1_`exposure'',`sd_alpha1_`exposure'')
		scalar alpha2star_`exposure' = rnormal(`mean_alpha2_`exposure'',`sd_alpha2_`exposure'')
		
		* SUM OVER EXPOSURE VARIABLE(S) TO GENERATE LINEAR PREDICTOR alpha1star*x AND alpha2star*x
		quietly replace `xb_utilde1' = `xb_utilde1' + scalar(alpha1star_`exposure')*`exposure'
		quietly replace `xb_utilde2' = `xb_utilde2' + scalar(alpha2star_`exposure')*`exposure'
	}
	
	* DRAW VALUE FOR BIAS PARAMETER eta1
	`eta1prior'
	
	* DRAW VALUE FOR BIAS PARAMETER pi2
	`pi2prior'
	
	* DRAW VALUES FOR BIAS PARAMETERS betaucon1 AND betaubin2
	scalar betastar_ucon1 = rnormal(`mean_beta_ucon1',`sd_beta_ucon1')	
	scalar betastar_ubin2 = rnormal(`mean_beta_ubin2',`sd_beta_ubin2')	
		
	/******************************************************************************************************	
	  DERIVE THE INTERCEPT OF logit{Pr(U2=1|X)} = alpha20 + alpha2x*X
		  INITIAL VALUES FOR INTERVAL SEARCH OF alpha20; 10 UNITS AROUND pi2star ON THE LOG-ODDS SCALE
	*******************************************************************************************************/
	local interval 10
	local logodds = ln(scalar(pi2star)/(1-scalar(pi2star)))
	local lower = `logodds' - `interval'/2
	local upper = `logodds' + `interval'/2

	* SEARCH FOR AN ESTIMATE OF alphacons THAT RESULTS IN A PR(sim_utilde2=1) WITHIN 1 DECIMAL PLACE OF pi2star
		* SET TO A RELATIVELY LARGE VALUE TO ENTER while LOOP
	local smallabsdiff 5
	local count_while 0
	while `smallabsdiff'>0.001 & `count_while'<10 {
		local ++count_while
		
		* SET TO A IMPOSSIBLE LARGE VALUE FOR FIRST CYCLE OF FOR-LOOP
		local absdiff_nearest 500
		
		* SIMULATE sim_utilde2 FOR DIFFERENT VALUES OF alphacons (I.E., alpha20) WITHIN GIVEN INTERVAL
				* UPDATE INTERVAL TO SEARCH WITHIN [lower,upper]; ALSO USED TO CALCULATE LOWER AND UPPER LIMITS FOR THE NEXT CYCLE BELOW
		local interval = `interval'/10
		forvalues alphacons=`lower'(`interval')`upper' {
			quietly replace `pr_sim_utilde2' = invlogit(`alphacons' + `xb_utilde2')
			quietly replace `sim_utilde2' = rbinomial(1,`pr_sim_utilde2')
			
			* CALCULATE PR(sim_utilde2=1) 
			quietly count if `sim_utilde2' !=.
			local denom = r(N)
			quietly count if `sim_utilde2'==1
			local current = r(N)/`denom'
			
			* CALCUATE ABSOLUTE DIFFERENCE FROM pistar
			local absdiff_current = abs(`current'-scalar(pi2star))

			* UPDATE IF THE CURRENT VALUE OF alpha20 GIVES A PREVALENCE CLOSER TO pi2star THAN AVAILABLE BEST ESTIMATE
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
	
	* SET alpha20 TO THE CLOSEST VALUE FROM THE WHILE LOOP
	local alpha20 `alphacons_nearest'

	/******************************
	  SIMULATE utilde1 AND utilde2
	*******************************/
	* SIMULATE sim_utilde2 
	quietly replace `pr_sim_utilde2' = invlogit(`alpha20' + `xb_utilde2')
	quietly replace `sim_utilde2' = rbinomial(1,`pr_sim_utilde2')

	* SIMULATE sim_utilde1
	quietly replace `sim_utilde1' = `xb_utilde1' + scalar(eta1star)*rnormal() 		

	/**************************************************************************************************************
	  FIT ANALYSIS OF INTEREST y|x,c,utilde1,utilde2 - CONSTRAINED REGRESSION
	  CONSTRAIN COEFFICIENT OF utilde1 AND utilde2 TO DRAWN VALUES betastar_ucon1 AND betastar_ubin2, RESPECTIVELY
	***************************************************************************************************************/
	* SET UP CONSTRAINTS
	capture constraint drop _all
	constraint 1 `sim_utilde1' = scalar(betastar_ucon1) 
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