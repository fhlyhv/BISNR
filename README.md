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
Below we provide an example of how to call the funtion BISN to learn time-varying graphical models from the synthetic data associated with the package. To test the algorithm on your own data, please replace data with your own data when calling BISN. Note that data is a n x p matrix with n and p variables.
```r
library(BISNR)

# load data and true presicion matrix in the package
data("data")
data("Ktrue")

# learn graphical model structure using BISN
results = BISN(data)

# check performance
precision = sum(results$Ksparse[lower.tri(results$Ksparse)] != 0 & Ktrue[lower.tri(Ktrue)] != 0)/sum(results$Ksparse[lower.tri(results$Ksparse)] != 0)
recall = sum(results$Ksparse[lower.tri(results$Ksparse)] != 0 & Ktrue[lower.tri(Ktrue)] != 0)/sum(Ktrue[lower.tri(Ktrue)] != 0)
f1_score = 2*precision*recall/(precision + recall)
cat("precision = ",precision,", recall = ", recall, ", f1_score = ", f1_score, "run_time = ", results$run_time)
```
