#' MCMC Small Area Estimation
#'
#' Fit multi-level models with possibly correlated random effects using MCMC.
#'
#' Functions to fit multi-level models with Gaussian and binomial or negative
#' binomial likelihoods using MCMC. Models with a linear predictor consisting
#' of various possibly correlated random effects are supported, allowing
#' flexible modeling of temporal, spatial or other kinds of dependence structures.
#' For Gaussian models the variance can be modeled too. By modeling unit-level
#' variances it is possible to account for outliers.
#' The package has been developed with applications to small area estimation
#' in official statistics in mind. The posterior samples for the model
#' parameters can be passed to a prediction function to generate samples from
#' the posterior predictive distribution for user-defined quantities such as
#' finite population domain means. For model assessment, posterior predictive
#' checks and DIC/WAIC criteria are supported.
#'
#' @name mcmcsae-package
#' @aliases mcmcsae
#' @docType package
NULL

#' @importFrom Rcpp evalCpp
#' @useDynLib mcmcsae, .registration=TRUE
NULL

# other namespace imports
#' @importFrom Matrix .updateCHMfactor bandSparse bdiag coerce crossprod diag Diagonal
#'   drop0 invPerm isDiagonal KhatriRao Matrix nnzero rsparsematrix sparseMatrix
#' @importClassesFrom Matrix CHMfactor dCHMsimpl ddiMatrix CsparseMatrix dgCMatrix dsCMatrix
#'   generalMatrix sparseMatrix
#' @importMethodsFrom Matrix %*% as.matrix as.vector Cholesky colSums crossprod
#'   determinant diag rowSums solve t tcrossprod unname
## do not import which() S4 generic from Matrix package as it slows down normal use of which
## @rawNamespace import(Matrix, except = which)
#' @import GIGrvg
#' @importFrom matrixStats colLogSumExps colQuantiles colSds colVars rowVars
#' @importFrom graphics abline axis legend lines matplot pairs par plot
#'   plot.new points segments
#' @importFrom methods as new setMethod show
#' @importFrom stats acf as.formula density fitted make.link mvfft optim
#'   pnorm predict rbeta rbinom rchisq residuals rexp rgamma rnbinom rnorm
#'   runif rWishart sd setNames terms update.formula var weights
#' @importFrom utils modifyList object.size setTxtProgressBar str tail
#'   txtProgressBar
NULL
