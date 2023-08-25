# MRtrajectory
This repository contains STATA do-files in support of the paper "Two sample Mendelian Randomisation using an outcome from a multilevel model of disease progression" (M. Lawton et al. 2023).

We are interested in carrying out two sample Mendelian Randomisation (MR) to test whether an exposure is causally related to disease progression that is measured as a linear trajectory over time. Disease progression can be modelled using multilevel models where we are interested in the intercept (disease severity at time zero) and the slope (rate of change).  Two sample MR is a technique that uses genetic data as instrumental variables to avoid bias due to confounding. 

In this repository we carry out a simulation study for a multivariate method to carry out two sample MR using multivariate meta-analysis which is often used in meta-analysis of diagnostic studies when researchers are interested in both the sensitivity and specificity of a test. Our aim is to examine bias and coverage of both separate and joint confidence intervals using this approach along with a naive approach that carrys out the MR separately on the intercept and slope.

The files XXX contain the STATA do-files to run the simulations.  These simulations were carried out using STATA 17.


