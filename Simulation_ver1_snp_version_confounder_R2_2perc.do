



******************************************************************************
**
**  Uses same parameters as simulation_study_iv_test_v1
**
******************************************************************************


args bob seed

set rng mt64s

set rngstream `bob'

set seed `seed'

use matrix_file, clear

mkmat snp_*, matrix(beta_grs_trans)


clear

set obs 10000


*** text for SNP exposure effects

scalar  scpf_snp_1 = 0.08180
scalar  scpf_snp_2 = 0.05560
scalar  scpf_snp_3 = 0.06010
scalar  scpf_snp_4 = 0.04020
scalar  scpf_snp_5 = 0.04820
scalar  scpf_snp_6 = 0.04470
scalar  scpf_snp_7 = 0.04140
scalar  scpf_snp_8 = 0.03340
scalar  scpf_snp_9 = 0.03150
scalar  scpf_snp_10 = 0.03070
scalar  scpf_snp_11 = 0.03090
scalar  scpf_snp_12 = 0.04510
scalar  scpf_snp_13 = 0.04030
scalar  scpf_snp_14 = 0.03600
scalar  scpf_snp_15 = 0.03110
scalar  scpf_snp_16 = 0.02620
scalar  scpf_snp_17 = 0.02610
scalar  scpf_snp_18 = 0.02420
scalar  scpf_snp_19 = 0.02830
scalar  scpf_snp_20 = 0.02350
scalar  scpf_snp_21 = 0.02970
scalar  scpf_snp_22 = 0.02490
scalar  scpf_snp_23 = 0.06580
scalar  scpf_snp_24 = 0.02270
scalar  scpf_snp_25 = 0.03340
scalar  scpf_snp_26 = 0.02170
scalar  scpf_snp_27 = 0.04770
scalar  scpf_snp_28 = 0.02180
scalar  scpf_snp_29 = 0.02340
scalar  scpf_snp_30 = 0.02300
scalar  scpf_snp_31 = 0.02240
scalar  scpf_snp_32 = 0.01880
scalar  scpf_snp_33 = 0.02290
scalar  scpf_snp_34 = 0.02090
scalar  scpf_snp_35 = 0.02490
scalar  scpf_snp_36 = 0.02110
scalar  scpf_snp_37 = 0.02200
scalar  scpf_snp_38 = 0.02000
scalar  scpf_snp_39 = 0.02210
scalar  scpf_snp_40 = 0.01950
scalar  scpf_snp_41 = 0.02070
scalar  scpf_snp_42 = 0.01910
scalar  scpf_snp_43 = 0.02980
scalar  scpf_snp_44 = 0.04830
scalar  scpf_snp_45 = 0.02020
scalar  scpf_snp_46 = 0.02250
scalar  scpf_snp_47 = 0.01880
scalar  scpf_snp_48 = 0.01920
scalar  scpf_snp_49 = 0.02070
scalar  scpf_snp_50 = 0.02070
scalar  scpf_snp_51 = 0.01850
scalar  scpf_snp_52 = 0.01820
scalar  scpf_snp_53 = 0.01800
scalar  scpf_snp_54 = 0.04920
scalar  scpf_snp_55 = 0.01790
scalar  scpf_snp_56 = 0.01580
scalar  scpf_snp_57 = 0.02770
scalar  scpf_snp_58 = 0.01760
scalar  scpf_snp_59 = 0.03060
scalar  scpf_snp_60 = 0.03080
scalar  scpf_snp_61 = 0.01680
scalar  scpf_snp_62 = 0.01630
scalar  scpf_snp_63 = 0.01410
scalar  scpf_snp_64 = 0.01780
scalar  scpf_snp_65 = 0.02580
scalar  scpf_snp_66 = 0.01920
scalar  scpf_snp_67 = 0.01900
scalar  scpf_snp_68 = 0.01770
scalar  scpf_snp_69 = 0.01640
scalar  scpf_snp_70 = 0.01640
scalar  scpf_snp_71 = 0.01880
scalar  scpf_snp_72 = 0.01740
scalar  scpf_snp_73 = 0.01720
scalar  scpf_snp_74 = 0.03070
scalar  scpf_snp_75 = 0.02010
scalar  scpf_snp_76 = 0.01670
scalar  scpf_snp_77 = 0.02450
scalar  scpf_snp_78 = 0.01820
scalar  scpf_snp_79 = 0.03580
scalar  scpf_snp_80 = 0.01880
scalar  scpf_snp_81 = 0.01870
scalar  scpf_snp_82 = 0.01740
scalar  scpf_snp_83 = 0.01590
scalar  scpf_snp_84 = 0.01750
scalar  scpf_snp_85 = 0.03950
scalar  scpf_snp_86 = 0.01980
scalar  scpf_snp_87 = 0.03170
scalar  scpf_snp_88 = 0.02100
scalar  scpf_snp_89 = 0.01940
scalar  scpf_snp_90 = 0.02170
scalar  scpf_snp_91 = 0.03500
scalar  scpf_snp_92 = 0.01670
scalar  scpf_snp_93 = 0.02330
scalar  scpf_snp_94 = 0.01920
scalar  scpf_snp_95 = 0.01720
scalar  scpf_snp_96 = 0.01970
scalar  scpf_snp_97 = 0.01660


*** text to simulate SNP data (number of effect alleles)

gen snp_1 = rbinomial(2, .415)
gen snp_2 = rbinomial(2, .2357)
gen snp_3 = rbinomial(2, .828)
gen snp_4 = rbinomial(2, .4341)
gen snp_5 = rbinomial(2, .1927)
gen snp_6 = rbinomial(2, .1773)
gen snp_7 = rbinomial(2, .7923)
gen snp_8 = rbinomial(2, .6125)
gen snp_9 = rbinomial(2, .3843)
gen snp_10 = rbinomial(2, .4617)
gen snp_11 = rbinomial(2, .4025)
gen snp_12 = rbinomial(2, .8719)
gen snp_13 = rbinomial(2, .865)
gen snp_14 = rbinomial(2, .8038)
gen snp_15 = rbinomial(2, .7841)
gen snp_16 = rbinomial(2, .4067)
gen snp_17 = rbinomial(2, .6289)
gen snp_18 = rbinomial(2, .4459)
gen snp_19 = rbinomial(2, .6664)
gen snp_20 = rbinomial(2, .5267)
gen snp_21 = rbinomial(2, .1955)
gen snp_22 = rbinomial(2, .32)
gen snp_23 = rbinomial(2, .0397)
gen snp_24 = rbinomial(2, .394)
gen snp_25 = rbinomial(2, .1334)
gen snp_26 = rbinomial(2, .5229)
gen snp_27 = rbinomial(2, .0722)
gen snp_28 = rbinomial(2, .583)
gen snp_29 = rbinomial(2, .7129)
gen snp_30 = rbinomial(2, .6816)
gen snp_31 = rbinomial(2, .7)
gen snp_32 = rbinomial(2, .7227)
gen snp_33 = rbinomial(2, .2867)
gen snp_34 = rbinomial(2, .6462)
gen snp_35 = rbinomial(2, .2106)
gen snp_36 = rbinomial(2, .352)
gen snp_37 = rbinomial(2, .7238)
gen snp_38 = rbinomial(2, .5816)
gen snp_39 = rbinomial(2, .273)
gen snp_40 = rbinomial(2, .555)
gen snp_41 = rbinomial(2, .6416)
gen snp_42 = rbinomial(2, .4285)
gen snp_43 = rbinomial(2, .2029)
gen snp_44 = rbinomial(2, .0655)
gen snp_45 = rbinomial(2, .5526)
gen snp_46 = rbinomial(2, .2651)
gen snp_47 = rbinomial(2, .5479)
gen snp_48 = rbinomial(2, .6198)
gen snp_49 = rbinomial(2, .2833)
gen snp_50 = rbinomial(2, .3033)
gen snp_51 = rbinomial(2, .5749)
gen snp_52 = rbinomial(2, .5748)
gen snp_53 = rbinomial(2, .6864)
gen snp_54 = rbinomial(2, .042)
gen snp_55 = rbinomial(2, .5416)
gen snp_56 = rbinomial(2, .4234)
gen snp_57 = rbinomial(2, .8791)
gen snp_58 = rbinomial(2, .454)
gen snp_59 = rbinomial(2, .1532)
gen snp_60 = rbinomial(2, .0891)
gen snp_61 = rbinomial(2, .4051)
gen snp_62 = rbinomial(2, .4211)
gen snp_63 = rbinomial(2, .3654)
gen snp_64 = rbinomial(2, .6312)
gen snp_65 = rbinomial(2, .8484)
gen snp_66 = rbinomial(2, .3199)
gen snp_67 = rbinomial(2, .293)
gen snp_68 = rbinomial(2, .3958)
gen snp_69 = rbinomial(2, .3934)
gen snp_70 = rbinomial(2, .6116)
gen snp_71 = rbinomial(2, .6879)
gen snp_72 = rbinomial(2, .4782)
gen snp_73 = rbinomial(2, .5246)
gen snp_74 = rbinomial(2, .9013)
gen snp_75 = rbinomial(2, .211)
gen snp_76 = rbinomial(2, .3913)
gen snp_77 = rbinomial(2, .152)
gen snp_78 = rbinomial(2, .6685)
gen snp_79 = rbinomial(2, .9101)
gen snp_80 = rbinomial(2, .3593)
gen snp_81 = rbinomial(2, .748)
gen snp_82 = rbinomial(2, .3653)
gen snp_83 = rbinomial(2, .5085)
gen snp_84 = rbinomial(2, .6092)
gen snp_85 = rbinomial(2, .0517)
gen snp_86 = rbinomial(2, .251)
gen snp_87 = rbinomial(2, .1798)
gen snp_88 = rbinomial(2, .197)
gen snp_89 = rbinomial(2, .7456)
gen snp_90 = rbinomial(2, .8117)
gen snp_91 = rbinomial(2, .9156)
gen snp_92 = rbinomial(2, .5339)
gen snp_93 = rbinomial(2, .1421)
gen snp_94 = rbinomial(2, .7469)
gen snp_95 = rbinomial(2, .4559)
gen snp_96 = rbinomial(2, .173)
gen snp_97 = rbinomial(2, .304)



matrix score total_snp_effects = beta_grs_trans

*** add in a confounder

gen confounder = rnormal(0,0.81147503)

*** see do-file calculating_variance.do-file
*** for details on calculation of residual variance, for proportion of variance explained 10%

gen exp_res = rnormal(0,0.81147503)

*** the minus part is to ensure that the exposure is centred

gen exp = total_snp_effects + confounder + exp_res - 2.3279252


*** create subject ID

gen person = _n

*** create individual level random effects

matrix U = (84.88533, -5.99176\ -5.99176, 4.608265)
drawnorm u0 u1, cov(U)


*** expand and create disease duration

expand 7
bysort person: generate disease_duration = _n - 1

*** level 1 residuals (rnormal is sd not variance)

di sqrt(51.97112)

gen lev1_res = rnormal(0,7.2090998)


*** generate data
*** 10% effect on intercept (per sd) is 2.0.  Sd is 0.51843238 so effect is 2.0/0.518
*** 10% effect on slope (per sd) 0.45

gen UPDRS = 24 + u0 + (2.25 + u1)*disease_duration + (1/0.81147503)*confounder + (0.225/0.81147503)*confounder*disease_duration + (2.0/1.159)*exp + (0.45/1.159)*exp*disease_duration + lev1_res

foreach sam of numlist 1/97{
gen snp_`sam'Xdd =	snp_`sam'*disease_duration
}

*** need postfile loop here to export the results for each snp and snp against disease duration
*** might be more efficient to also export the snp-exposure effects here too by creating scalars above somewhere

tempname memhold

postfile `memhold' beta_exp beta_int se_int beta_sl se_sl cov_int_sl using simulation_study_snp_ver1_confounders_R2_2perc_results_dataset_`bob', replace

foreach sam of numlist 1/97{
mixed UPDRS  disease_duration snp_`sam' snp_`sam'Xdd ||  person: disease_duration, cov(unstruct)
matrix define mat_cov = e(V)
scalar n_int = colnumb(mat_cov,"snp_`sam'")
scalar n_sl = colnumb(mat_cov,"snp_`sam'Xdd")
scalar cov_int_sl = (mat_cov[n_int,n_sl])
post `memhold' (scpf_snp_`sam') (_b[snp_`sam']) (_se[snp_`sam']) (_b[snp_`sam'Xdd]) (_se[snp_`sam'Xdd]) (cov_int_sl) 	
}

postclose `memhold'


