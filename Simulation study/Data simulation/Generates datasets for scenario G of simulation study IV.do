****************************************************************************************************************
* Code to simulate survival time-to-event data, scenario G of simulation study IV
* There is uninformative censoring but no delayed entry
* One continuous exposure, Xcon
* One continuous unmeasured confounder Ucon
* Three measured confounders; one continuous (Ccon) one nominal(Cnom) and one binary (Cbin)
* Simulate (1) Ccon, Cnom and Cbin as independent and exogenous
* Simulate (2) X ~ Ccon + Cnom + Cbin
* Simulate (3) Ucon ~ X + Ccon + Cnom + Cbin
* Simulate (4) Y ~ X + Ccon + Cnom + Cbin + Ucon
* Generates 500 simulated datasets, each of N=1000, called "Dataset_1.csv", etc
* Column headings of cnom cbin ccon xcon ucon1 y_t y_d

* Before running to do-file:
	* Install user-written command survsim. 
	* Provide file path (working directory) for location of saved datasets - update line 26 
		* Datasets are saved in a folder called Data. This must not already exist in the working directory

* David Carslake, July 2025
****************************************************************************************************************

*******************
* A. Preliminaries:
*******************
clear all
*local working_directory = [enter filepath of working directory]
mkdir "`working_directory'\ystxcon3"
cd "`working_directory'\ystxcon3"

********************************************
* B. Get some parameters from the real data:
********************************************
* Simulation parameters were derived from University of Bristol teaching dataset pbc1bas.dta
* This in turn is based on data used in Christensen et al 1985 Gastroenterology 89, 1084-1091
*  Outcome = survival time in years from diagnosis with primary biliary cholangitis
*  Xcon = logb0 (Bilirubin at baseline)
*  Ccon = alb0 (Albumin at baseline), Cnom=logigm (IgM at baseline), Cbin = cenc0 (central cholestasis at baseline)
*  U = age (Age at baseline)
*  Time axis = time in study; assume uniform recruitment over 12 years with administrative censoring at the end


********************************************
* C. Define a command to simulate a dataset:
********************************************
global c_args = "b_ccon_cons(real) sd_ccon(real) t1_cbin(real) t1_cnom(real) t2_cnom(real) t3_cnom(real)"
global x_args = "b_xcon_cons(real) b_xcon_ccon(real) b_xcon_cbin(real) b_xcon_cnom2(real) b_xcon_cnom3(real) b_xcon_cnom4(real) sd_xcon(real)"
global u_args = "b_ucon_cons(real) b_ucon_xcon(real) b_ucon_ccon(real) b_ucon_cbin(real) b_ucon_cnom2(real) b_ucon_cnom3(real) b_ucon_cnom4(real) sd_ucon(real)"
global y_args = "b_yst_cons(real) shape_yst(real) b_yst_xcon(real) b_yst_ccon(real) b_yst_cbin(real) b_yst_cnom2(real) b_yst_cnom3(real) b_yst_cnom4(real) b_yst_ucon(real)"
capture program drop Simdata_Xcon
program Simdata_Xcon, rclass
	syntax , fileno(integer) n(integer) maxfu(real) $c_args $x_args $u_args $y_args
	display "Simulating dataset number `fileno' with Xcon"
	quietly{
		* 1. Simulate exogenous Ccon, Cnom and Cbin:
		capture drop *
		set obs `n'
		generate Ccon = rnormal(`b_ccon_cons',`sd_ccon')
		generate Cbin = cond(runiform()<`t1_cbin',0,1)
		generate rands = runiform()
		generate Cnom = .
		replace Cnom = 1 if rands>=0         & rands<`t1_cnom'
		replace Cnom = 2 if rands>=`t1_cnom' & rands<`t2_cnom'
		replace Cnom = 3 if rands>=`t2_cnom' & rands<`t3_cnom'
		replace Cnom = 4 if rands>=`t3_cnom' & rands<.
		drop rands
		tabulate Cnom, generate(Cnom)
		* 2. Simulate X ~ Ccon + Cnom + Cbin
		generate Xcon = `b_xcon_cons' + `b_xcon_ccon'*Ccon + `b_xcon_cbin'*Cbin + `b_xcon_cnom2'*Cnom2 + `b_xcon_cnom3'*Cnom3 + `b_xcon_cnom4'*Cnom4 + rnormal(0,`sd_xcon')
		* 3. Simulate U ~ X + Ccon + Cnom + Cbin
		generate Ucon = `b_ucon_cons' + `b_ucon_xcon'*Xcon + `b_ucon_ccon'*Ccon + `b_ucon_cbin'*Cbin + `b_ucon_cnom2'*Cnom2 + `b_ucon_cnom3'*Cnom3 + `b_ucon_cnom4'*Cnom4 + rnormal(0,`sd_ucon')
		* 4. Simulate Y ~ X + Ccon + Cnom + Cbin + U
		* (The shape parameter gamma in survsim is p in streg. The scale parameter lambda in survsim is _cons, in its exponentiated form, in streg)
		local covariates = `"Xcon `b_yst_xcon' Ccon `b_yst_ccon' Cbin `b_yst_cbin' Cnom2 `b_yst_cnom2' Cnom3 `b_yst_cnom3' Cnom4 `b_yst_cnom4' Ucon `b_yst_ucon'"'
		survsim T_event, lambdas(`b_yst_cons') gammas(`shape_yst') distribution(weibull) covariates(`covariates')
		* 5. Administratively censor some people, assuming uniform recruitment over maxfu years from 2000 and censoring after that:
		generate Date_recruited = 2000+runiform(0,`maxfu')
		generate Date_event = Date_recruited+T_event
		generate Date_censor = 2000+`maxfu'
		generate T_censor = Date_censor-Date_recruited
		generate y_d = cond(Date_event<Date_censor,1,0)
		generate y_t = cond(Date_event<Date_censor,T_event,T_censor)
		drop Date_* T_*
		* 6. Write the data file:
		drop Cnom1 Cnom2 Cnom3 Cnom4
		order Cnom Cbin Ccon Xcon Ucon y_t y_d
		rename Cnom cnom
		rename Cbin cbin
		rename Ccon ccon
		rename Xcon xcon
		rename Ucon ucon
		export delimited "Dataset_`fileno'.csv", delimiter(",") replace
	}
end
* Example code with values taken from real data:
*local c_opts = "b_ccon_cons(-0.058) sd_ccon(5.9) t1_cbin(0.8) t1_cnom(0.25) t2_cnom(0.5) t3_cnom(0.75)"
*local x_opts = "b_xcon_cons(-0.13) b_xcon_ccon(-0.014) b_xcon_cbin(0.50) b_xcon_cnom2(0.063) b_xcon_cnom3(0.084) b_xcon_cnom4(0.058) sd_xcon(0.40)"
*local u_opts = "b_ucon_cons(55) b_ucon_xcon(-2.6) b_ucon_ccon(-0.54) b_ucon_cbin(0.85) b_ucon_cnom2(-1.4) b_ucon_cnom3(-0.29) b_ucon_cnom4(-0.22) sd_ucon(10)"
*local y_opts = "b_yst_cons(0.0062) shape_yst(1.4) b_yst_xcon(2.2) b_yst_ccon(-0.043) b_yst_cbin(0.70) b_yst_cnom2(-0.27) b_yst_cnom3(-0.13) b_yst_cnom4(-0.019) b_yst_ucon(0.04)"
*Simdata_Xcon, fileno(1) n(1000) maxfu(12) `c_opts' `x_opts' `u_opts' `y_opts'


*************************************************
* D. Create the simulated data folders and files:
*************************************************
* (The effect of Xcon on Y was set to zero and the effect of Xcon on Ucon was tripled to increase bias)
local N_datasets = 500
set seed 123456
local c_opts = "b_ccon_cons(-0.058) sd_ccon(5.9) t1_cbin(0.8) t1_cnom(0.25) t2_cnom(0.5) t3_cnom(0.75)"
local x_opts = "b_xcon_cons(-0.13) b_xcon_ccon(-0.014) b_xcon_cbin(0.50) b_xcon_cnom2(0.063) b_xcon_cnom3(0.084) b_xcon_cnom4(0.058) sd_xcon(0.40)"
local u_opts = "b_ucon_cons(55) b_ucon_xcon(-7.8) b_ucon_ccon(-0.54) b_ucon_cbin(0.85) b_ucon_cnom2(-1.4) b_ucon_cnom3(-0.29) b_ucon_cnom4(-0.22) sd_ucon(10)"
local y_opts = "b_yst_cons(0.0062) shape_yst(1.4) b_yst_xcon(0) b_yst_ccon(-0.043) b_yst_cbin(0.70) b_yst_cnom2(-0.27) b_yst_cnom3(-0.13) b_yst_cnom4(-0.019) b_yst_ucon(0.04)"
forvalues ds=1/`N_datasets'{
	Simdata_Xcon, fileno(`ds') n(1000) maxfu(12) `c_opts' `x_opts' `u_opts' `y_opts'	
}

* END OF DO-FILE