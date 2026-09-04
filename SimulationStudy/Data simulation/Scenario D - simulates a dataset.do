/*************************************************************************************************************************************************
                                           SIMULATES A SINGLE DATASET FOR SCENARIO D OF SIMULATION STUDY III
										   
SCENARIO D: Y CONTINUOUS, X CONTINUOUS, MULTIPLE U=(U1,U2) BINARY, C=(CONTINUOUS ccon, BINARY cbin, AND NOMINAL cnom)

ARGUMENTS:
		samplesize       - SAMPLE SIZE OF EACH SIMULATED DATASET
		strength_UC      - MULTIPLIER OF ASSOCIATIONS BETWEEN U AND C (0 INDEPENDENCE BETWEEN U AND C, 1 ASSOCIATIONS OBSERVED IN REAL DATASET,
		                   2 DOUBLE THE OBSERVED ASSOCIATIONS)
		pilevel		     - PREVALENCES OF U1 AND U2 (1:Pr(U1=1)~14.1% AND Pr(U2=1)=18%)

RETURNS: 
	A SINGLE SIMULATED DATASET CONTAINING Y, X, U, C							   
*************************************************************************************************************************************************/
args samplesize strength_UC pilevel	

noisily di "Simulate a dataset of sample size `samplesize' with `strength_UC' x U-C associations"

if `pilevel'==1 {		
	noisily di "Prevalence of U1 is close to 14.1 and prevalence of U2 is close to 18%"	
	
	local beta0 8.534
	
    if `strength_UC'==0 {
		local alpha0 -6.695
		local gamma0 -6.395
	}
	else if `strength_UC'==1 {
		local alpha0 -8.118
		local gamma0 -6.406
	}
	else if `strength_UC'==2 {
		local alpha0 -9.563
		local gamma0 -6.42
	}
}
else if `pilevel'==0 {  // 
}
else if `pilevel'==2 {  // 
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

* BINARY UNMEASURED CONFOUNDER ubin1
capture drop pr_ubin1 ubin1
local alphacnom1 = `strength_UC'*-0.1704 
local alphacnom2 = `strength_UC'*-0.02826  
local alphacbin = `strength_UC'*0.1006
local alphaccon = `strength_UC'*0.4337
local alphaxcon 0.2500
gen pr_ubin1 = invlogit(`alpha0' + `alphaxcon'*xcon + `alphacnom1'*cnom1 + `alphacnom2'*cnom2 + `alphacbin'*cbin + `alphaccon'*ccon)
gen ubin1 = rbinomial(1,pr_ubin1)

* BINARY UNMEASURED CONFOUNDER ubin2
capture drop pr_ubin2 ubin2
local gammacnom1 = `strength_UC'*-0.08843  
local gammacnom2 = `strength_UC'*-0.1626  
local gammacbin = `strength_UC'*0.007311
local gammaccon = `strength_UC'*0.02073
local gammaxcon 0.2500
local gammaubin1 = `strength_UC'*0.1495
gen pr_ubin2 = invlogit(`gamma0' + `gammaxcon'*xcon + `gammacnom1'*cnom1 + `gammacnom2'*cnom2 + `gammacbin'*cbin + `gammaccon'*ccon + ///
                        `gammaubin1'*ubin1)
gen ubin2 = rbinomial(1,pr_ubin2)

* CONTINUOUS OUTCOME ycon
capture drop ycon
local betacnom1 0.5421      
local betacnom2 0.5630      
local betacbin -0.2531		
local betaccon -0.1783
local betaxcon 0
local betaubin1 2.5
local betaubin2 1.5
local epsilony 4.322
gen ycon = `beta0' + `betacnom1'*cnom1 + `betacnom2'*cnom2 + `betacbin'*cbin + `betaccon'*ccon + `betaxcon'*xcon + `betaubin1'*ubin1 + /// 
           `betaubin2'*ubin2 + `epsilony'*rnormal() 

drop cnom1 cnom2		   
		   
* END OF DO FILE