# BISNR (R package for BISN)

This R package implements the BISN algorithm proposed in [1]. The Matlab C++ Mex code that implements the same algorithm can be found at https://github.com/fhlyhv/BISN.

## Dependence
Please make sure to install the following package dependencies before using R package `BISNR`. 
```r
install.packages(c("Rcpp", "RcppArmadillo", "BH", "KernSmooth", "tictoc", "devtools"))
```

## Installation
The R package `BISNR` can be installed from source files in the GitHub repository (R package `devtools` is needed):
```r
library(devtools)
install_github(repo="fhlyhv/BISNR")
```

If it does not work, please download this repository directly as a zip or tar.gz file and install the package locally.

## Example
Below we provide an example of how to call the funtion BISN to learn the graphical model structure from the synthetic data associated with the package. 
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

In the above the example, we simply call the funtion BISN as `results = BISN(data)`, where data is a n x p matrix with n observations for each of the p variables. Missing data can be represented by NA. The resulting results$Ksparse matrix will be a sparse matrix.

You need to reduce the step size upper bound `eta` (see below) if the algorithm divergences (e.g., some very large values suddenly appears). By default, we set `eta = 300`.
```r
results = BISN(data, eta = 100)
```
On the other hand, you may consider increasing `eta` if the algorithm doesn't diverge and you want to speed up the convergence.

Instead of simply estmating L and D from the data, the BISN function here further thresholding <&lambda;<sub>jk</sub>> / (1 +<&lambda;<sub>jk</sub>>) using the method in 
Section V in [3] to yield a sparse graph in an automated manner. However, due to the mean-filed approximation used in BISN, <&lambda;<sub>jk</sub>> of 
elements (j, k) in the bottom-right corner are typically not well estimated when the sample size is small. More specifically, since K<sub>jk</sub> = [LDL<sup>T</sup>]<sub>jk</sub>, the elements in the bottom-right corner of K are the sum of a larger set of elements in L and D than the elements in the top-left corner. Due to the mean field 
approximation, the estimates of <K<sub>jk</sub><sup>2</sup>> is typically corrupted by the variances of elements in L and D. These variances can be large when the sample size is small. As <&lambda;<sub>jk</sub>> is a function of <K<sub>jk</sub><sup>2</sup>>, it cannot be well estimated 
either. To settle this problem, we run the original BISN algorithm again by reversely ordering the data (i.e., setting options.backward_pass = 1) and then average the 
resulting <&lambda;<sub>jk</sub>> with that from the forward pass. Note that `backward_pass = TRUE` by default and it can be set to `FALSE` when the sample 
size is large.

In addition, the kernel density of <&lambda;<sub>jk</sub>> / (1 +<&lambda;<sub>jk</sub>>) might be bumpy when the sample size is small. In this case, it is better to 
plot out the density and choose the threshold of <&lambda;<sub>jk</sub>> manually.

On the other hand, after estimating the sparse precision matrix, the BISN function can further reestimate the non-zero elements in the precision via maximum likelihood by setting `prm_learning = TRUE`. BISN can reliably estimate the non-zero elements when the sample size is relatively large, but it is recommended to reestimate the non-zero elements when the sample size is small. 

In addition to the sparse K matrix, there are other output parameters. Please refer to the help file of the function BISN for more details.

[1] H. Yu, S. Wu, L. Xin, and J. Dauwels. Fast Bayesian Inference of Sparse Networks with Automatic Sparsity Determination. Journal of Machine Learning Research, 2020.

[2] C. Sanderson, R. Curtin. Armadillo: a template-based C++ library for linear algebra. Journal of Open Source Software, 2016.

[3] H. Yu, L. Xin, and J. Dauwels. Variational wishart approximation for graphical model selection: Monoscale and multiscale models. IEEE Transactions on Signal Processing, 2019.
