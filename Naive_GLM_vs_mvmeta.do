

**************************************************************
***
**  Naive approach to MR procedure - calculating the intercept effect on the 97 SNPs
**  This could use GLM approach or meta-regression approach
**  In GLM approach we want to scale the variance to be 1 first of all for consistency with meta-regeression approach
**  This section of code will create 1,000 new files (one for each simulated dataset) 
**  These new files will have two variables 
**  beta_intercept = "point estimate for effect that exposure has on the intercept"
**  se_intercept = "standard error for effect that exposure has on the intercept"
***
**************************************************************


cd "$datapath"


tempname memhold

postfile `memhold' beta_intercept se_intercept using naive_intercept_${version}, replace


*** loop around all 1000 simulations
foreach bob of numlist 1/1000{
di `bob' 
use simulation_study_snp_${version}_results_dataset_`bob', clear    
glm beta_int beta_exp [iw=1/(se_int^2)], scale(1) nocons    
post `memhold' (_b[beta_exp]) (_se[beta_exp])	
}


postclose `memhold'



**************************************************************
***
**  Naive approach to MR procedure - calculating the slope effect
**  This could use GLM approach or meta-regression approach
**  In GLM approach we want to scale the variance to be 1 first of all for consistency with meta-regeression approach
**  This section of code will create 1,000 new files (one for each simulated dataset)
**  These new files will have two variables 
**  beta_slope = "point estimate for effect that exposure has on the intercept"
**  se_slope = "standard error for effect that exposure has on the intercept"
***
**************************************************************


cd "$datapath"


tempname memhold

postfile `memhold' beta_slope se_slope using naive_slope_${version}, replace


*** loop around all 1000 simulations
foreach bob of numlist 1/1000{
di `bob' 
use simulation_study_snp_${version}_results_dataset_`bob', clear    
glm beta_sl beta_exp [iw=1/(se_sl^2)], scale(1) nocons    
post `memhold' (_b[beta_exp]) (_se[beta_exp])	
}


postclose `memhold'



**************************************************************
***
**  Bivariate approach to MR procedure - mvmeta
**  This section of code will create 1,000 new files (one for each simulated dataset)
**  These new files will have five variables
**  beta_intercept = "point estimate for effect that exposure has on the intercept"
**  se_intercept = "standard error for effect that exposure has on the intercept" 
**  beta_slope = "point estimate for effect that exposure has on the intercept"
**  se_slope = "standard error for effect that exposure has on the intercept"
**  covariance = "covariance between the estimate of the effect of exposure on the intercept and the estimate of the effect the exposure has on the slope"
***
**************************************************************


cd "$datapath"


tempname memhold

postfile `memhold' beta_intercept se_intercept beta_slope se_slope covariance using bivariate_${version}, replace



*** loop around all 1000 simulations
foreach bob of numlist 1/1000{
    di `bob' 
use simulation_study_snp_${version}_results_dataset_`bob', clear 

*** renaming the intercept and slope effects (and covariance matrix) to be in mvmeta format (need numbers not words)

rename beta_int beta1
rename beta_sl beta2
rename beta_exp exposure

gen V11 = se_int^2
gen V22 = se_sl^2
gen V12 = cov_int_sl   

mvmeta beta V exposure, noconstant fixed



post `memhold' ([beta1]_b[exposure]) ([beta1]_se[exposure])	 ([beta2]_b[exposure]) ([beta2]_se[exposure]) (e(V)[1,2])	
}


postclose `memhold'