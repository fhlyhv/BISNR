#' @title Bayesian Inference of Sparse Networks
#' @param data n x p matrix of observed data. n is the sample size and p is the dimension. Missing data can be denoted by NA.
#' @param eta upper bound of the step size (eta >0, 300 by default).
#' @param max_iter maximum number of iterations (1e4 by default).
#' @param tol tolerance to check convergence of the algorithm (1e-2 by default).
#' @param r the decaying factor (0 < r < 1, 0.5 by default)
#' @param s minibatch size (s = p/(1e-3*(p-1)+1)) by default)
#' @param backward_pass boolean value to decide whether to enable the backward pass or not (1 by default). 
#' backward_pass = 1 would run the BISN algorithm again by reversely ordering the data (i.e., the backward pass), 
#' and then average the results from both forward and backward pass. This will improve the estimation accuracy, 
#' especially when the sample size is small.
#' @param prm_learning boolean value to decide whether to reestimate the non-zero elements via maximum 
#' likelihood or not. (0 by default)
#' @return Ksparse p x p matrix with the same zero pattern as the estimated adjacency matrix Adj. 
#' The nonzero elements in Ksparse is reestimated by maximum likelihood if prm_learning = TRUE.
#' @return Adj estimated adjacency matrix using the method in Section V in (Yu et al, Variational wishart approximation for 
#' graphical model selection: Monoscale and multiscale models, 2019).
#' @return Kest p x p full matrix. Kest = ML * MD * ML', where MD and ML denotes the mean of the D and L matrix.
#' @return Lambda p x p estimated Lambda matrix.
#' @return run_time total running time.
#' @export
BISN <- function(data, eta = 300, max_iter = 1e4, tol = 1e-2, r = 0.5, s = ncol(data) / (1e-3 * (ncol(data) - 1) + 1),
                 backward_pass = TRUE, prm_learning = FALSE) {

  if (!is.matrix(data) && !is.data.frame(data)) {
    stop("The input data should be either a matrix or a data frame")
  } else if (is.data.frame(data)) {
    data = as.matrix(data)
  }
  set.seed(0)
  
  tictoc::tic()
  p = ncol(data)
  Lambda = matrix(0, p, p)
  if (backward_pass) Lambda1 = matrix(0, p, p)
  idl = which(lower.tri(Lambda, diag = FALSE))

  id_rc = which(is.na(data), arr.ind = TRUE)

  if (length(id_rc) > 0) {
    if (backward_pass) id_rc1 = which(is.na(data[, p : 1]), arr.ind = TRUE)
    data[is.na(data)] = 0
    cat("forward pass ... \n")
    row_missing = unique(id_rc[, 1])
    results = BISN_missing(data, row_missing, id_rc, eta, max_iter, tol, r, s)
    Kest = (results$ML * rep(results$mD[, 1], rep(p, p))) %*% t(results$ML)
    Lambda[idl] = results$lambda
    rm(results)
    Lambda = Lambda + t(Lambda)

    if (backward_pass) {
      cat("backward pass ...\n")
      data = data[, p : 1]
      row_missing = unique(id_rc1[, 1])
      results = BISN_missing(data, row_missing, id_rc1, eta, max_iter, tol, r, s)
      Kest1 = (results$ML * rep(results$mD[, 1], rep(p, p))) %*% t(results$ML)
      Lambda1[idl] = results$lambda
      rm(results)
      Lambda1 = Lambda1 + t(Lambda1)
      Kest1 = Kest1[p : 1, p : 1]
      Lambda1 = Lambda1[p : 1, p : 1]

      Kest = (Kest + Kest1) / 2;
      Lambda = (Lambda + Lambda1) / 2;

    }

  } else {
    cat("forward pass ... \n")
    results = BISN_obsv(data, eta, max_iter, tol, r, s)
    Kest = (results$ML * rep(results$mD[, 1], rep(p, p))) %*% t(results$ML)
    Lambda[idl] = results$lambda
    rm(results)
    Lambda = Lambda + t(Lambda)

    if (backward_pass) {
      cat("backward pass ...\n")
      data = data[, p : 1]
      results = BISN_obsv(data, eta, max_iter, tol, r, s)
      Kest1 = (results$ML * rep(results$mD[, 1], rep(p, p))) %*% t(results$ML)
      Lambda1[idl] = results$lambda
      rm(results)
      Lambda1 = Lambda1 + t(Lambda1)
      Kest1 = Kest1[p : 1, p : 1]
      Lambda1 = Lambda1[p : 1, p : 1]

      Kest = (Kest + Kest1) / 2;
      Lambda = (Lambda + Lambda1) / 2;

    }
  }

  cat("forward-backward pass is done, ")
  t = tictoc::toc()

  cat("estimate adjacency matrix by thresholding lambda / (1 + lambda)...\n")
  lambda = Lambda[idl]
  ll = lambda / (1 + lambda)
  xy = KernSmooth::bkde(ll,  bandwidth =  KernSmooth::dpik(ll))
  #density(ll, bw = "SJ-ste", kernel = "gaussian", n = 128)
  x = xy$x
  y = xy$y
  thr = x[y == min(y[x<0.6 & x > 1e-2])][1]
  thr = thr / (1 - thr)
  Adj = Lambda < thr
  Ksparse = Kest
  Ksparse[!Adj] = 0
  cat("adjacency marix has been estimated, ")
  t = tictoc::toc()

  if (prm_learning) {
    cat("start reestimating the non-zero elements...\n")
    if (backward_pass) data = data[, p : 1]
    if (length(id_rc) > 0) {
      n = nrow(data)
      id_missing = id_rc[, 1] + (id_rc[, 2] - 1) * n
      obsv_mat = matrix(1, n, p)
      obsv_mat[id_missing] = 0
      data[id_missing] = 0
      S = t(data) %*% data / (t(obsv_mat) %*% obsv_mat - 1)
    } else S = cov(data)
    id_nonzero = which(Adj, arr.ind = TRUE)
    id_nonzero = id_nonzero[id_nonzero[, 1] > id_nonzero[, 2], ]
    Ksparse = QUICParameterLearning(Ksparse, S, id_nonzero[, 1], id_nonzero[, 2], 100, 100);
    cat("reestimating the non-zero elements is done, ")
    t = tictoc::toc()
  }

  list(Ksparse = Ksparse, Adj = Adj, Lambda = Lambda, Kest = Kest, run_time = t$toc - t$tic)

}
