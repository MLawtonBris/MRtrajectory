

**************************************************************
***
** Do-file to create all of the statistics that are reported in tables 2 and 3 of the main manuscript 
***
**************************************************************


* arguments to be passed to the do-file

args datapath version true_int true_slope adjust_term


* datapath = folder that the simulation results are stored
* version = version of the simulations (called scenario rather than version in the manuscript)
* true_int = true parameter for the effect that the exposure has on the intercept (per sd change in the exposure)
* true_slope = true parameter for the effect that the exposure has on the slope (per sd change in the exposure)
* adjust_term = expected sd of the exposure. This adjustment term is important to put the estimates on the standardised scale of the exposure allowing comparibility between the R^2=10% and R^2=2% simulations 



*************************************
*** NAIVE INTERCEPT
*************************************

cd "`datapath'"


use naive_intercept_`version', clear

*** create a log file with summary of results

log using naive_intercept_`version'_summary_results.log, replace 

foreach bob of varlist beta_intercept se_intercept {
replace `bob' = `bob' * `adjust_term'
}

*** summary of estimate for the intercept

summ beta_intercept 

*** mean relative bias (on percentage scale)

di ((r(mean)-`true_int')/(`true_int'))*100

*** generating and summarising the bias

gen intercept_bias = beta_intercept - `true_int'

summ intercept_bias

*** calculate monte carlo standard error for the bias

di r(sd)/sqrt(r(N))

*** summary of the model based standard error

summ se_intercept

*** generate lower and upper part of 95% confidence interval

gen lower_int = beta_intercept - (abs(invnormal(0.025))*se_intercept)
gen upper_int = beta_intercept + (abs(invnormal(0.025))*se_intercept)

*** count number of instances where truth is within the 95% CI and calculate coverage percentage

count if lower_int < `true_int' & upper_int > `true_int'
di (r(N)/1000)*100


log close

*************************************
*** NAIVE SLOPE
*************************************

use naive_slope_`version', clear

*** create a log file with summary of results

log using naive_slope_`version'_summary_results.log, replace 

foreach bob of varlist beta_slope se_slope {
replace `bob' = `bob' * `adjust_term'
}

*** summary of estimate for the slope 

summ beta_slope

*** mean relative bias (on percentage scale)

di ((r(mean)-`true_slope')/(`true_slope'))*100

*** generating and summarising the bias

gen slope_bias = beta_slope - `true_slope'

summ slope_bias

*** calculate monte carlo standard error for the bias

di r(sd)/sqrt(r(N))

*** summary of the model based standard error

summ se_slope

*** generate lower and upper part of 95% confidence interval

gen lower_slope = beta_slope - (abs(invnormal(0.025))*se_slope)
gen upper_slope = beta_slope + (abs(invnormal(0.025))*se_slope)

*** count number of instances where truth is within the 95% CI and calculate coverage percentage

count if lower_slope < `true_slope' & upper_slope > `true_slope'
di (r(N)/1000)*100

log close


*************************************
*** BIVARIATE INTERCEPT AND SLOPE
*************************************

use bivariate_`version', clear

*** create a log file with summary of results

log using bivariate_`version'_summary_results.log, replace 

foreach bob of varlist beta_intercept se_intercept beta_slope se_slope{
replace `bob' = `bob' * `adjust_term'
}

*** INTERCEPT section
*** summary of estimate for the intercept

summ beta_intercept 

*** mean relative bias (on percentage scale)

di ((r(mean)-`true_int')/(`true_int'))*100

*** generating and summarising the bias

gen intercept_bias = beta_intercept - `true_int'

summ intercept_bias

*** calculate monte carlo standard error for the bias

di r(sd)/sqrt(r(N))

*** summary of the model based standard error

summ se_intercept

*** generate lower and upper part of 95% confidence interval

gen lower_int = beta_intercept - (abs(invnormal(0.025))*se_intercept)
gen upper_int = beta_intercept + (abs(invnormal(0.025))*se_intercept)

*** count number of instances where truth is within the 95% CI and calculate coverage percentage

count if lower_int < `true_int' & upper_int > `true_int'
di (r(N)/1000)*100


*** SLOPE section
*** summary of estimate for the slope

summ beta_slope

*** mean relative bias (on percentage scale)

di ((r(mean)-`true_slope')/(`true_slope'))*100

*** generating and summarising the bias

gen slope_bias = beta_slope - `true_slope'

summ slope_bias

*** calculate monte carlo standard error for the bias

di r(sd)/sqrt(r(N))

*** summary of the model based standard error

summ se_slope

*** generate lower and upper part of 95% confidence interval

gen lower_slope = beta_slope - (abs(invnormal(0.025))*se_slope)
gen upper_slope = beta_slope + (abs(invnormal(0.025))*se_slope)

*** count number of instances where truth is within the 95% CI and calculate coverage percentage

count if lower_slope < `true_slope' & upper_slope > `true_slope'
di (r(N)/1000)*100

log close





