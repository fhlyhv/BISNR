# BISNR (R package for BISN)

This R package implements the BISN algorithm proposed in [1]. 

## Dependence
Please make sure to install the following package dependencies before using R package `BISNR`. 
```r
install.packages(c("Rcpp", "RcppArmadillo", "BH", "RcppProgress", "KernSmooth"))
```

## Installation
The R package `BISNR` can be installed from source files in the GitHub repository (R package `devtools` is needed):
```r
library(devtools)
install_github(repo="fhlyhv/BISNR")
```

If it does not work, please download this repository directly as a zip or tar.gz file and install the package locally.

## Example
Below we provide an example of how to call the funtion BADGE to learn time-varying graphical models from the synthetic data associated with the package. To test the algorithm on your own data, please replace data with your own data when calling BADGE. Note that data is a N x P matrix with N time points and P variables.
