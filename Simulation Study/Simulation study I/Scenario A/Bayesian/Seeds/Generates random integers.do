/*******************************************************************************
GENERATES 8,000 RANDOM INTEGERS BETWEEN 1 AND 2147483647
*******************************************************************************/

clear

set obs 8000

* SET SEED
set seed 20241221

* GENERATE RANDOM INTEGERS
gen seed = int(runiform(1, 2147483647))
sum seed

* CHECK FOR DUPLICATES
duplicates report

* SAVE TO A CSV FILE
export delimited "RandomIntegers.csv", replace