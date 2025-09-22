/*************************************************************************************************************************************************
                                           SIMULATES A SINGLE DATASET FOR SCENARIO E OF SIMULATION STUDY IV
										   
SCENARIO E: Y CONTINUOUS, X NOMINAL, MULTIPLE U=(CONTINUOUS U1, BINARY U2), C=(CONTINUOUS ccon, BINARY cbin, AND NOMINAL cnom)

ARGUMENTS:
		samplesize       - SAMPLE SIZE OF EACH SIMULATED DATASET
		strength_UC      - MULTIPLIER OF ASSOCIATIONS BETWEEN U AND C (0 INDEPENDENCE BETWEEN U AND C, 1 ASSOCIATIONS OBSERVED IN REAL DATASET,
		                   2 DOUBLE THE OBSERVED ASSOCIATIONS)
		pilevel		     - PREVALENCE OF U2 (1:Pr(U2=1)=59.77%)
						   						   
RETURNS: 
	A SINGLE SIMULATED DATASET CONTAINING Y, X, U, C							   
*************************************************************************************************************************************************/
args samplesize strength_UC pilevel

noisily di "Simulate a dataset of sample size `samplesize' with `strength_UC' x U-C associations"

if `pilevel'==1 {		
	noisily di "Prevalence of U2 is 59.77%"	
	local beta0 8.534
    if `strength_UC'==1 {
		local alpha10 60.67
		local gamma20 -0.5634
	}
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

* CATEGORICAL EXPOSURE xnom
capture drop xnom exp_xb0 exp_xb1
* xnom=0|C
local theta00 4.020
local theta0ccon -2.173
local theta0cbin -0.1836
local theta0cnom1 0.2612
local theta0cnom2 0.5005
* xnom=1|C
local theta10 3.040
local theta1ccon -1.205
local theta1cbin -0.1018
local theta1cnom1 0.1508
local theta1cnom2 0.2860
gen double exp_xb0 = exp(`theta00' + `theta0cnom1'*cnom1 + `theta0cnom2'*cnom2 + `theta0cbin'*cbin + `theta0ccon'*ccon)
gen double exp_xb1 = exp(`theta10' + `theta1cnom1'*cnom1 + `theta1cnom2'*cnom2 + `theta1cbin'*cbin + `theta1ccon'*ccon)
gen double denominator = 1 + exp_xb0 + exp_xb1
gen double pr_xnom0 = exp_xb0/denominator
gen double pr_xnom1 = exp_xb1/denominator
gen double pr_xnom01 = pr_xnom0 + pr_xnom1
replace random  = uniform()
gen xnom = cond(random < pr_xnom0, 0, cond(random < pr_xnom01, 1, 2))
xi i.xnom, noomit   // CREATES DUMMY VARIABLES _Icnom_1 AND _Icnom_2
rename _Ixnom_0 xnom0
rename _Ixnom_1 xnom1

* CONTINUOUS UNMEASURED CONFOUNDER ucon1
capture drop ucon1 
local alpha1xnom0 -7.492
local alpha1xnom1 -4.175
local alpha1cnom1 = `strength_UC'*-1.039
local alpha1cnom2 = `strength_UC'*-.3788
local alpha1ccon = `strength_UC'*3.055
local alpha1cbin = `strength_UC'*0.5880
local eta1 9.554
gen ucon1 = `alpha10' + `alpha1cnom1'*cnom1 + `alpha1cnom2'*cnom2 + `alpha1cbin'*cbin + `alpha1ccon'*ccon + `alpha1xnom0'*xnom0 + ///
            `alpha1xnom1'*xnom1 + `eta1'*rnormal()

* BINARY UNMEASURED CONFOUNDER ubin2
capture drop ubin2 pr_ubin2
local gamma2xnom0 -1.084
local gamma2xnom1 -0.5996
local gamma2cnom1 = `strength_UC'*-0.08549
local gamma2cnom2 = `strength_UC'*-0.1699
local gamma2ccon = `strength_UC'*0.1163
local gamma2cbin = `strength_UC'*0.01151
local gamma2ucon1 = `strength_UC'*0.01258

gen pr_ubin2 = invlogit(`gamma20' + `gamma2cnom1'*cnom1 + `gamma2cnom2'*cnom2 + `gamma2cbin'*cbin + `gamma2ccon'*ccon + `gamma2xnom0'*xnom0 + ///
                        `gamma2xnom1'*xnom1 + `gamma2ucon1'*ucon1)
gen ubin2 = rbinomial(1,pr_ubin2)

* CONTINUOUS OUTCOME ycon
capture drop ycon
local betaxnom0 0
local betaxnom1 0
local betaccon -0.1783
local betacbin -0.2531		
local betacnom1 0.5421      
local betacnom2 0.5630      
local betaucon1 0.1299
local betaubin2 0.5
local sigma 4.322

gen ycon = `beta0' + `betacnom1'*cnom1 + `betacnom2'*cnom2 + `betacbin'*cbin + `betaccon'*ccon + `betaxnom0'*xnom0 + `betaxnom1'*xnom1 + ///
           `betaucon1'*ucon1 + `betaubin2'*ubin2 + `sigma'*rnormal() 

* END OF DO FILE