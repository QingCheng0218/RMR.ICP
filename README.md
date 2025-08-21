RMR.ICP
=======
  
  **RMR.ICP** is a package for Robust Mendelian Randomization method to account for Idiosyncratic and Correlated Pleiotropy.

Installation
============
  Install the development version of **RMR.ICP** by use of the 'devtools' package. Note that **RMR.ICP** depends on the 'Rcpp' and 'RcppArmadillo' package, which also requires appropriate setting of Rtools and Xcode for Windows and Mac OS/X, respectively.
```
library(devtools)
install_github("QingCheng0218/RMR.ICP@main")
```

If you have errors installing this package on Linux, try the following commands in R:
  ```
Sys.setlocale(category = "LC_ALL", locale = "English_United States.1252") 
Sys.setenv(LANG="en_US.UTF-8")
Sys.setenv(LC_ALL="en_US.UTF-8")

library(devtools)
install_github("QingCheng0218/RMR.ICP@main")
```

Examples
=========
```
library(RMR.ICP);
# Run RMR.ICP using independent SNPs.
result4indep = RMRICPindep(gamma1_Indep, Gamma2_Indep, se1_Indep, se2_Indep);
beta1 = result4indep$beta.hat # the causal effect estimate, 
beta1_se = result4indep$beta.se # the standard deviation of causal effect.
beta1_pvalue = result4indep$beta.p.value # p-value calculated by RMR-ICP is two-sided based on the posterior distribution of causal effect.
ProEta = result4indep$Etares # $$Pr(\eta=1|text{Data})$$, the posterior probabilities of IVs with idiosyncratic correlated horizontal pleiotropy effects.

# Run RMR.ICP using correlated SNPs.
RMRICPsim = RMRICPSim(gamma1_LD, Gamma2_LD, se1_LD, se2_LD, R, block_inf1);
beta1 = RMRICPsim$beta.hat
beta1_se = RMRICPsim$beta.se
beta1_pvalue = RMRICPsim$beta.p.value
ProEta = RMRICPsim$Etares

```
References
==========

Development
===========
  
  This package is developed and maintained by Qing Cheng (qingcheng0218@gmail.com). 
  
# RMR.ICP
