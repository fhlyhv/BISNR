\name{BISNR-package}
\alias{BISNR-package}
\alias{BISNR}
\docType{package}
\title{
\packageTitle{BISNR}
}
\description{
\packageDescription{BISNR}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{BISNR}
\packageIndices{BISNR}
This R package implements the BISN algorithm proposed in (Yu et al, 2020). For more details of the algorithm, please check the paper.
}
\author{
\packageAuthor{BISNR}

Maintainer: \packageMaintainer{BISNR}
}
\references{
H. Yu, S. Wu, L. Xin, and J. Dauwels. Fast Bayesian Inference of Sparse Networks with Automatic Sparsity Determination. Journal of Machine Learning Research, 2020.

H. Yu, L. Xin, and J. Dauwels. Variational wishart approximation for graphical model selection: Monoscale and multiscale models. IEEE Transactions on Signal Processing, 2019.
}

\keyword{graphical models, structure learning, variational inference, graphical lasso}
\seealso{
~~ Optional links to other man pages, e.g. ~~
~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
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
}
