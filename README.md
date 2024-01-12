# MRtrajectory
This repository contains STATA do-files in support of the paper "Two sample Mendelian Randomisation using an outcome from a multilevel model of disease progression" (M. Lawton et al. 2023).

We are interested in carrying out two sample Mendelian Randomisation (2SMR) to test whether an exposure is causally related to disease progression that is measured as a linear trajectory over time. Disease progression can be modelled using multilevel models where we are interested in the intercept (disease severity at time zero) and the slope (rate of change).  Two sample MR is a technique that uses genetic data as instrumental variables to avoid bias due to confounding. 

In this repository we carry out a simulation study for a multivariate method to carry out two sample MR using multivariate meta-analysis. This type of meta-analysis is often used in meta-analysis of diagnostic studies when researchers are interested in both the sensitivity and specificity of a test. Our aim is to examine bias and coverage of both separate and joint confidence intervals using this approach along with a naive approach that carrys out the 2SMR separately on the intercept and slope.

The files Simulation_ver*_* contain the STATA do-files to run the simulations.  These simulations were carried out using STATA 17.

The "Naive_GLM_vs_mvmeta.do" do-file carries out the naive approach to 2SMR and the multivariate approach to 2SMR on the simulated data.  

The "Naive_GLM_vs_mvmeta_summarise_results_tables_2_and_3.do" summarises the results from the two approaches to 2SMR.  This includes reporting the point estimates, mean model based SE, coverage of the 95% confidence intervals, mean relative bias along with the bias (and its Monte Carlo Standard error) across the 1000 simulations. 

The "Naive_GLM_vs_mvmeta_joint_coverage.do" calculates the joint coverage of the different approaches.  That is constructing a confidence rectangle for the naive approach, a confidence ellipse for the naive approach (assuming independence) and a confidence ellipse for the bivariate approach.  This is reported in table 4 of the main manuscript.

The "Naive_GLM_vs_mvmeta_joint_coverage_calculate_area_of_ellipse.do" calculates the area of the confidence ellipses using the naive and bivariate approach.  This is reported in supplementary table 1.

