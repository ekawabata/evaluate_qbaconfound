/*************************************************************************************************************************************************
                                           SIMULATES A SINGLE DATASET FOR SCENARIO F OF SIMULATION STUDY IV
										   
SCENARIO F: Y NOMINAL, X CONTINUOUS, SINGLE U=(CONTINUOUS U1, BINARY U2), C=(CONTINUOUS ccon, BINARY cbin, AND NOMINAL cnom)

ARGUMENTS:
		samplesize       - SAMPLE SIZE OF EACH SIMULATED DATASET
		strength_UC      - MULTIPLIER OF ASSOCIATIONS BETWEEN U AND C (0 INDEPENDENCE BETWEEN U AND C, 1 ASSOCIATIONS OBSERVED IN REAL DATASET,
						   						   
RETURNS: 
	A SINGLE SIMULATED DATASET CONTAINING Y, X, U, C							   
*************************************************************************************************************************************************/
args samplesize strength_UC

noisily di "Simulate a dataset of sample size `samplesize' with `strength_UC' x U-C associations"

if `strength_UC'==1 {
	local alpha0 58.11
}

clear
set obs `samplesize'

* MULTINOMIAL DISTRIBUTION OF NOMINAL CONFOUNDER cnom (0,1,2)
local pi0 0.1261
local pi1 0.7477
local pi2 0.1261
local pi0pi1 = `pi0' + `pi1'
local pi0pi1pi2 = `pi0' + `pi1' + `pi2'
capture drop random
capture drop cnom
gen random  = uniform()
gen cnom = cond(random < `pi0', 0, cond(random < `pi0pi1', 1, 2))
tab cnom, m
xi i.cnom   // CREATES DUMMY VARIABLES _Icnom_1 AND _Icnom_2
rename _Icnom_1 cnom1
rename _Icnom_2 cnom2

* BINARY CONFOUNDER cbin
capture drop pr_cbin cbin
local omega0 0.2624
local omegacnom1 -0.3651
local omegacnom2 -0.4075
gen pr_cbin = invlogit(`omega0' + `omegacnom1'*cnom1 + `omegacnom2'*cnom2)
gen cbin = rbinomial(1,pr_cbin)

* CONTINUOUS MEASURED CONFOUNDER ccon
capture drop ccon
local psi0 2.916
local psicnom1 0.5791
local psicnom2 0.6344
local psicbin -0.1057
local psisd 0.4495
gen ccon = `psi0' + `psicnom1'*cnom1 + `psicnom2'*cnom2 + `psicbin'*cbin + `psisd'*rnormal()

* CONTINUOUS EXPOSURE xcon
capture drop xcon
local theta0 13.41
local thetaccon 1.719
local thetacbin 0.1442
local thetacnom1 -0.1988
local thetacnom2 -0.3895
local thetasd 2.176
gen xcon = `theta0' + `thetacnom1'*cnom1 + `thetacnom2'*cnom2 + `thetacbin'*cbin + `thetaccon'*ccon + `thetasd'*rnormal()

* CONTINUOUS UNMEASURED CONFOUNDER ucon
capture drop ucon
local alphaxcon 1.281
local alphaccon = `strength_UC'*0.1617
local alphacbin = `strength_UC'*0.05031
local alphacnom1 = `strength_UC'*-0.4942
local alphacnom2 = `strength_UC'*-0.8886
local eta 10.30
gen ucon = `alpha0' + `alphacnom1'*cnom1 + `alphacnom2'*cnom2 + `alphacbin'*cbin + `alphaccon'*ccon + `alphaxcon'*xcon + `eta'*rnormal()

* NOMINAL OUTCOME ynom
capture drop exp_xb0 exp_xb1 pr_ynom0 pr_ynom1 pr_ynom01 denominator ynom
* COEFFICIENTS FOR CATEGORY ynom=0
local beta0xcon 0
local beta0cnom1 -0.1022
local beta0cnom2 -0.1217
local beta0cbin 0.04580
local beta0ccon -0.02349
local beta0ucon 0.1
local beta00 -7
* COEFFICIENTS FOR CATEGORY ynom=1
local beta1xcon 0
local beta1cnom1 -0.1959
local beta1cnom2 -0.2413  
local beta1cbin 0.08650
local beta1ccon -0.05253 
local beta1ucon 0.1 
local beta10 -7

gen double exp_xb0 = exp(`beta00' + `beta0cnom1'*cnom1 + `beta0cnom2'*cnom2 + `beta0cbin'*cbin + `beta0ccon'*ccon + `beta0xcon'*xcon + `beta0ucon'*ucon)
gen double exp_xb1 = exp(`beta10' + `beta1cnom1'*cnom1 + `beta1cnom2'*cnom2 + `beta1cbin'*cbin + `beta1ccon'*ccon + `beta1xcon'*xcon + `beta1ucon'*ucon)
gen double denominator = 1 + exp_xb0 + exp_xb1
gen double pr_ynom0 = exp_xb0/denominator
gen double pr_ynom1 = exp_xb1/denominator
gen double pr_ynom01 = pr_ynom0 + pr_ynom1
replace random  = uniform()
gen ynom = cond(random < pr_ynom0, 0, cond(random < pr_ynom01, 1, 2))

keep ynom xcon cnom1 cnom2 ccon cbin ucon	
	
* END OF DO FILE