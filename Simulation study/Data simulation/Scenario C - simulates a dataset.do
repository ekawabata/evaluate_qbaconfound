/***************************************************************************************************************************************************
                                           SIMULATES A SINGLE DATASET FOR SCENARIO C OF SIMULATION STUDY III
										   
SCENARIO C: Y BINARY, X BINARY, SINGLE U=(U1) BINARY, C=(CONTINUOUS ccon, BINARY cbin, AND NOMINAL cnom)

ARGUMENTS:
		samplesize       - SAMPLE SIZE OF EACH SIMULATED DATASET
		strength_UC      - MULTIPLIER OF ASSOCIATIONS BETWEEN U AND C (0 INDEPENDENCE BETWEEN U AND C, 1 ASSOCIATIONS OBSERVED IN REAL DATASET,
		                   2 DOUBLE THE OBSERVED ASSOCIATIONS)
		pilevel		     - PREVALENCE OF U (0:~6.5%, 1:~20.5%, 2:~40%)

RETURNS: 
	A SINGLE SIMULATED DATASET CONTAINING Y, X, U, C							   
*************************************************************************************************************************************************/
args samplesize strength_UC pilevel

noisily di "Simulate a dataset of sample size `samplesize' with `strength_UC' x U-C associations"

clear
set obs `samplesize'

if `pilevel'==1 {		// PR(U=1)=~20.5%
	noisily di "Prevalence of U1 is close to 20.5%"	
	
	local beta0 -2.370		
	
    if `strength_UC'==0 {
		local alpha0 -1.65
	}
	else if `strength_UC'==1 {
		local alpha0 -3.633	
	}
	else if `strength_UC'==2 {
		local alpha0 -5.65
	}
	else if `strength_UC'==3 {
		local alpha0 -7.7
	}
}
else if `pilevel'==0 {  // PR(U=1)=~6.5%
	noisily di "Prevalence of U1 is close to 6.5%"	
	
	local beta0 -2.075		
	
    if `strength_UC'==0 {
		local alpha0 -3  // 6.55
	}
	else if `strength_UC'==1 {
		local alpha0 -5  
	}
	else if `strength_UC'==2 {
		local alpha0 -7.1
	}
	else if `strength_UC'==3 {
		local alpha0 -9.2
	}
}
else if `pilevel'==2 {  // PR(U=1)=~40%
	noisily di "Prevalence of U1 is close to 40%"
	
	local beta0 -2.715
	
    if `strength_UC'==0 {
		local alpha0 -0.65
	}
	else if `strength_UC'==1 {
		local alpha0 -2.6   
	}
	else if `strength_UC'==2 {
		local alpha0 -4.6
	}
	else if `strength_UC'==3 {
		local alpha0 -6.55
	}
}

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
xi i.cnom   
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

* BINARY EXPOSURE xbin
capture drop pr_xbin xbin
local theta0 -4.819
local thetacnom1 -0.3487
local thetacnom2 -0.7207
local thetacbin 0.4609
local thetaccon 1.024
gen pr_xbin = invlogit(`theta0' + `thetacnom1'*cnom1 + `thetacnom2'*cnom2 + `thetacbin'*cbin + `thetaccon'*ccon)
gen xbin = rbinomial(1,pr_xbin)

* BINARY UNMEASURED CONFOUNDER ubin
capture drop pr_ubin ubin
local alphacnom1 = `strength_UC'*-0.1639
local alphacnom2 = `strength_UC'*-0.0004130
local alphacbin = `strength_UC'*0.04386
local alphaccon =`strength_UC'*0.6078
local alphaxbin 1.178
gen pr_ubin = invlogit(`alpha0' + `alphacnom1'*cnom1 + `alphacnom2'*cnom2 + `alphacbin'*cbin + `alphaccon'*ccon + `alphaxbin'*xbin)
gen ubin = rbinomial(1,pr_ubin)

* BINARY OUTCOME ybin
capture drop pr_ybin ybin
local betaxbin 0
local betaccon -0.009095
local betacbin -0.3559
local betacnom1 0.1753
local betacnom2 0.04492
local betaubin 1.5
gen pr_ybin = invlogit(`beta0' + `betacnom1'*cnom1 + `betacnom2'*cnom2 + `betacbin'*cbin + `betaccon'*ccon + `betaxbin'*xbin + `betaubin'*ubin)
gen ybin = rbinomial(1,pr_ybin)

* END OF DO FILE