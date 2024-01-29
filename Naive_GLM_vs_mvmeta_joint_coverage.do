

**************************************************************
***
**   Do-file to calculate the joint coverage of the different approaches
**   Good examples and descriptions of equations for confidence ellipses here
**   https://geostatisticslessons.com/lessons/errorellipses
**   https://www.stat.cmu.edu/~larry/=stat401/lecture-18.pdf
**   Schubert, P. and Kirchner, M., 2014. Ellipse area calculations and their applicability in posturography. Gait & Posture, 39(1), pp.518-522.
***
**************************************************************



* arguments to be passed to the do-file

args datapath version true_int true_slope adjust_term


* datapath = folder that the simulation results are stored
* version = version of the simulations (called scenario rather than version in the manuscript)
* true_int = true parameter for the effect that the exposure has on the intercept (per sd change in the exposure)
* true_slope = true parameter for the effect that the exposure has on the slope (per sd change in the exposure)
* adjust_term = expected sd of the exposure. This adjustment term is important to put the estimates on the standardised scale of the exposure allowing comparibility between the R^2=10% and R^2=2% simulations 



**********************************************************************************
*** naive approach 
**********************************************************************************



cd "`datapath'"


use naive_intercept_`version', clear

merge 1:1 _n using naive_slope_`version'

drop _merge


*** individual confidence intervals for intercept and slope

gen lower_int = beta_intercept - (abs(invnormal(0.025))*se_intercept)
gen upper_int = beta_intercept + (abs(invnormal(0.025))*se_intercept)



gen lower_slope = beta_slope - (abs(invnormal(0.025))*se_slope)
gen upper_slope = beta_slope + (abs(invnormal(0.025))*se_slope)


count 
local total = r(N)

*** naive approach confidence rectangle

log using joint_coverage_rectangle_naive_`version'.log, replace 

*** note that alternative is to rescale the point estimates and standard errors

count if (lower_slope < (`true_slope'/`adjust_term') & upper_slope > (`true_slope'/`adjust_term')) & (lower_int < (`true_int'/`adjust_term') & upper_int > (`true_int'/`adjust_term'))
di (r(N)/1000)*100


log close

*** adapting old code so renaming variables to be consistent with old file!

rename beta_slope beta_sl
rename beta_intercept beta_int
rename se_slope se_sl
rename se_intercept se_int

*** no covariance estimated in naive approach  - set to zero assuming independence

gen cov_int_sl = 0

*** matrix approach to joint coverage (ellipse)

*** recode all variables to be on sd change scale
*** note alternative approach below starts here

replace se_sl = se_sl*`adjust_term'
replace se_int = se_int*`adjust_term'
replace beta_sl = beta_sl*`adjust_term'
replace beta_int = beta_int*`adjust_term'
replace cov_int_sl = cov_int_sl*(`adjust_term'^2)



forvalues i=1/1000{
	matrix y`i' = (beta_int[`i'] \ beta_sl[`i'])
	mata y`i' = st_matrix("y`i'") 
	matrix true`i' = ( `true_int' \ `true_slope')
	mata true`i' = st_matrix("true`i'") 
	matrix V`i' = (se_int[`i']^2 , cov_int_sl[`i'] \ cov_int_sl[`i'] , se_sl[`i']^2)
	mata V`i' = st_matrix("V`i'")
	mata invV`i' = invsym(V`i')
	mata test`i' = (y`i' - true`i')' * invV`i' * (y`i' - true`i')
	if `i' == 1 {
		mata teststack = test1
	}
	else {
		mata teststack = (teststack \ test`i')
	}
}

getmata teststack, force

*** note alternative approach ends here

/* note alternative approach is to only adjust the true estimates by sd change instead of the variables in the dataset

local true_int = `true_int'/`adjust_term'
local true_slope = `true_slope'/`adjust_term'

forvalues i=1/1000{
	matrix y`i' = (beta_int[`i'] \ beta_sl[`i'])
	mata y`i' = st_matrix("y`i'") 
	matrix true`i' = ( $true_int \ $true_slope)
	mata true`i' = st_matrix("true`i'") 
	matrix V`i' = (se_int[`i']^2 , cov_int_sl[`i'] \ cov_int_sl[`i'] , se_sl[`i']^2)
	mata V`i' = st_matrix("V`i'")
	mata invV`i' = invsym(V`i')
	mata test`i' = (y`i' - true`i')' * invV`i' * (y`i' - true`i')
	if `i' == 1 {
		mata teststack = test1
	}
	else {
		mata teststack = (teststack \ test`i')
	}
}

getmata teststack, force


*/

log using joint_coverage_ellipse_naive_`version'.log, replace 

count if teststack <= invchi2tail(2,0.05)
di (r(N)/1000)*100

log close








*****************************************************
*** bivariate approach
*****************************************************



cd "`datapath'"



use bivariate_`version', clear




rename beta_slope beta_sl
rename beta_intercept beta_int
rename se_slope se_sl
rename se_intercept se_int
rename covariance cov_int_sl


*** recode all variables to be on sd change scale
*** note alternative approach below starts here

replace se_sl = se_sl*`adjust_term'
replace se_int = se_int*`adjust_term'
replace beta_sl = beta_sl*`adjust_term'
replace beta_int = beta_int*`adjust_term'
replace cov_int_sl = cov_int_sl*(`adjust_term'^2)



forvalues i=1/1000{
	matrix y`i' = (beta_int[`i'] \ beta_sl[`i'])
	mata y`i' = st_matrix("y`i'") 
	matrix true`i' = ( `true_int' \ `true_slope')
	mata true`i' = st_matrix("true`i'") 
	matrix V`i' = (se_int[`i']^2 , cov_int_sl[`i'] \ cov_int_sl[`i'] , se_sl[`i']^2)
	mata V`i' = st_matrix("V`i'")
	mata invV`i' = invsym(V`i')
	mata test`i' = (y`i' - true`i')' * invV`i' * (y`i' - true`i')
	if `i' == 1 {
		mata teststack = test1
	}
	else {
		mata teststack = (teststack \ test`i')
	}
}

getmata teststack, force

*** note alternative approach ends here

/* note alternative approach is to only adjust the true estimates by sd change instead of the variables in the dataset

local true_int = `true_int'/`adjust_term'
local true_slope = `true_slope'/`adjust_term'

forvalues i=1/1000{
	matrix y`i' = (beta_int[`i'] \ beta_sl[`i'])
	mata y`i' = st_matrix("y`i'") 
	matrix true`i' = ( $true_int \ $true_slope)
	mata true`i' = st_matrix("true`i'") 
	matrix V`i' = (se_int[`i']^2 , cov_int_sl[`i'] \ cov_int_sl[`i'] , se_sl[`i']^2)
	mata V`i' = st_matrix("V`i'")
	mata invV`i' = invsym(V`i')
	mata test`i' = (y`i' - true`i')' * invV`i' * (y`i' - true`i')
	if `i' == 1 {
		mata teststack = test1
	}
	else {
		mata teststack = (teststack \ test`i')
	}
}

getmata teststack, force


*/

log using joint_coverage_ellipse_bivariate_`version'.log, replace 

count if teststack <= invchi2tail(2,0.05)
di (r(N)/1000)*100

log close




