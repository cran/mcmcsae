% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/models.R
\name{compute_GMRF_matrices}
\alias{compute_GMRF_matrices}
\title{Compute (I)GMRF incidence, precision and restriction matrices corresponding to a generic model component}
\usage{
compute_GMRF_matrices(
  factor,
  data,
  D = TRUE,
  Q = TRUE,
  R = TRUE,
  cols2remove = NULL,
  remove.redundant.R.cols = TRUE,
  enclos = .GlobalEnv,
  n.parent = 1L,
  ...
)
}
\arguments{
\item{factor}{factor formula of a generic model component,
see \code{\link{gen}}.}

\item{data}{data frame to be used in deriving the matrices.}

\item{D}{if \code{TRUE} compute the incidence matrix.}

\item{Q}{if \code{TRUE} compute the precision matrix.}

\item{R}{if \code{TRUE} compute the restriction matrix.}

\item{cols2remove}{if an integer vector is passed, the dimensions (columns of D,
rows and columns of Q and rows of R) that are removed. This can be useful in the
case of empty domains.}

\item{remove.redundant.R.cols}{whether to test for and remove redundant restrictions from restriction matrix R}

\item{enclos}{enclosure to look for objects not found in \code{data}.}

\item{n.parent}{for internal use; in case of custom factor, the number of frames up
the calling stack in which to evaluate any custom matrices}

\item{...}{further arguments passed to \code{economizeMatrix}.}
}
\value{
A list containing some or all of the components \code{D} (incidence matrix),
 \code{Q} (precision matrix) and \code{R} (restriction matrix).
}
\description{
This function computes incidence, precision and restriction matrices, or
a subset thereof, for a Gaussian Markov Random Field (GMRF).
A GMRF is specified by a formula passed to the \code{factor} argument,
in the same way as for the \code{factor} argument of \code{\link{gen}}.
}
\examples{
n <- 1000
dat <- data.frame(
  x = rnorm(n),
  f1 = factor(sample(1:50, n, replace=TRUE)),
  f2 = factor(sample(1:10, n, replace=TRUE))
)
mats <- compute_GMRF_matrices(~ f1 * RW1(f2), dat)
str(mats)

}
