MR.CUE
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

Usage
=========

References
==========

Development
===========
  
  This package is developed and maintained by Qing Cheng (qingcheng0218@gmail.com). 
  
