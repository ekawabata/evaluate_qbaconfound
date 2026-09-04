/*************************************************************************************************************************************************
                                           SIMULATES A SINGLE DATASET FOR SCENARIO B OF SIMULATION STUDIES I AND II
										   
SCENARIO B: Y CONTINUOUS, X CONTINUOUS, MULTIPLE U=(U1,U2) CONTINUOUS, C=(CONTINUOUS ccon, BINARY cbin, AND NOMINAL cnom)

ARGUMENTS:
		samplesize       - SAMPLE SIZE OF EACH SIMULATED DATASET
		strength_UC      - MULTIPLIER OF ASSOCIATIONS BETWEEN U AND C (0 INDEPENDENCE BETWEEN U AND C, 1 ASSOCIATIONS OBSERVED IN REAL DATASET,
		                   2 DOUBLE THE OBSERVED ASSOCIATIONS)
RETURNS: 
	A SINGLE SIMULATED DATASET CONTAINING Y, X, U, C							   
*************************************************************************************************************************************************/
args samplesize strength_UC betaxcon

noisily di "Simulate a dataset of sample size `samplesize' with true exposure effect of `betaxcon' and `strength_UC' x U-C and `strength_UC' x U1-U2 associations"

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

* CONTINUOUS UNMEASURED CONFOUNDER ucon2
capture drop ucon2 mu2
local alpha2ccon = `strength_UC'*0.1617
local alpha2cbin = `strength_UC'*0.05031
local alpha2cnom1 = `strength_UC'*-0.4942
local alpha2cnom2 = `strength_UC'*-0.8886
local alpha2xcon 1.281
local alpha20 58.11
local eta2 10.30
gen mu2 = `alpha20' + `alpha2cnom1'*cnom1 + `alpha2cnom2'*cnom2 + `alpha2cbin'*cbin + `alpha2ccon'*ccon + `alpha2xcon'*xcon
gen ucon2 = mu2 + `eta2'*rnormal()

* CONTINUOUS UNMEASURED CONFOUNDER ucon1
capture drop mu1 mubar ucon1
local alpha1cnom1 = `strength_UC'*-0.8979
local alpha1cnom2 = `strength_UC'*-0.1473
local alpha1cbin = `strength_UC'*0.5159
local alpha1ccon = `strength_UC'*2.209
local alpha1xcon 1.260
local alpha10 37.91
local eta1 9.394
gen mu1 = `alpha10' + `alpha1cnom1'*cnom1 + `alpha1cnom2'*cnom2 + `alpha1cbin'*cbin + `alpha1ccon'*ccon + `alpha1xcon'*xcon
local rho12 = 4.502*`strength_UC'
local slope = `rho12'/(`eta2'^2)
gen mubar = mu1 + `slope'*(ucon2-mu2)
local epsilonsq = (`eta1'^2) - `slope'*`rho12'
gen ucon1 = mubar + sqrt(`epsilonsq')*rnormal()

* CONTINUOUS OUTCOME ycon
capture drop ycon
local betaccon -0.1783
local betacbin -0.2531		
local betacnom1 0.5421      
local betacnom2 0.5630  
local betaucon1 0.1299 
local betaucon2 0.06747
local beta0 8.534
local sigma 4.322
gen ycon = `beta0' + `betacnom1'*cnom1 + `betacnom2'*cnom2 + `betacbin'*cbin + `betaccon'*ccon + `betaxcon'*xcon + `betaucon1'*ucon1 + ///
           `betaucon2'*ucon2 + `sigma'*rnormal() 

* END OF DO FILE