

**************************************************************
***
**   Do-file to calculate the area of the confidence ellipses (supplementary table 1) 
**   Good examples and descriptions of equations for confidence ellipses here
**   https://geostatisticslessons.com/lessons/errorellipses
**   https://www.stat.cmu.edu/~larry/=stat401/lecture-18.pdf
**   Schubert, P. and Kirchner, M., 2014. Ellipse area calculations and their applicability in posturography. Gait & Posture, 39(1), pp.518-522.
***
**************************************************************


* arguments to be passed to the do-file

args datapath version adjust_term


* datapath = folder that the simulation results are stored
* version = version of the simulations (called scenario rather than version in the manuscript)
* adjust_term = expected sd of the exposure. This adjustment term is important to put the estimates on the standardised scale of the exposure allowing comparibility between the R^2=10% and R^2=2% simulations 



**********************************************************************************
*** naive approach 
**********************************************************************************



cd "`datapath'"


use naive_intercept_`version', clear

merge 1:1 _n using naive_slope_`version'

drop _merge




*** adapting old code so renaming variables to be consistent with old file!

rename beta_slope beta_sl
rename beta_intercept beta_int
rename se_slope se_sl
rename se_intercept se_int

*** no covariance estimated in naive approach - set to zero assuming independence

gen cov_int_sl = 0





replace se_sl = se_sl*`adjust_term'
replace se_int = se_int*`adjust_term'
replace beta_sl = beta_sl*`adjust_term'
replace beta_int = beta_int*`adjust_term'
replace cov_int_sl = cov_int_sl*(`adjust_term'^2)



forvalues i=1/1000{

	matrix V`i' = (se_int[`i']^2 , cov_int_sl[`i'] \ cov_int_sl[`i'] , se_sl[`i']^2)
	mata V`i' = st_matrix("V`i'")
	matrix eigenvalues re`i' im`i' = V`i'
	mata re`i' = st_matrix("re`i'") 
	mata ystar`i' = re`i'[1,1]
	mata xstar`i' = re`i'[1,2]
	if `i' == 1 {
		mata ystarstack = ystar1
		mata xstarstack = xstar1
	}
	else {
		mata ystarstack = (ystarstack \ ystar`i')
		mata xstarstack = (xstarstack \ xstar`i')
	}
}


getmata ystarstack xstarstack, force


log using area_confidence_ellipse_naive_`version'.log, replace 

gen area = _pi*invchi2tail(2,0.05)*(sqrt(ystarstack*xstarstack))
summ area

log close







*****************************************************
*** bivariate approach
*****************************************************




cd "`datapath'"




use bivariate_`version', clear



*** adapting old code so renaming variables to be consistent with old file!

rename beta_slope beta_sl
rename beta_intercept beta_int
rename se_slope se_sl
rename se_intercept se_int
rename covariance cov_int_sl





replace se_sl = se_sl*`adjust_term'
replace se_int = se_int*`adjust_term'
replace beta_sl = beta_sl*`adjust_term'
replace beta_int = beta_int*`adjust_term'
replace cov_int_sl = cov_int_sl*(`adjust_term'^2)



forvalues i=1/1000{

	matrix V`i' = (se_int[`i']^2 , cov_int_sl[`i'] \ cov_int_sl[`i'] , se_sl[`i']^2)
	mata V`i' = st_matrix("V`i'")
	matrix eigenvalues re`i' im`i' = V`i'
	mata re`i' = st_matrix("re`i'") 
	mata ystar`i' = re`i'[1,1]
	mata xstar`i' = re`i'[1,2]
	if `i' == 1 {
		mata ystarstack = ystar1
		mata xstarstack = xstar1
	}
	else {
		mata ystarstack = (ystarstack \ ystar`i')
		mata xstarstack = (xstarstack \ xstar`i')
	}
}


getmata ystarstack xstarstack, force


log using area_confidence_ellipse_bivariate_`version'.log, replace 

gen area = _pi*invchi2tail(2,0.05)*(sqrt(ystarstack*xstarstack))
summ area

log close



