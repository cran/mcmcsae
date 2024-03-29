% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MatrixUtils.R
\name{aggrMatrix}
\alias{aggrMatrix}
\title{Utility function to construct a sparse aggregation matrix from a factor}
\usage{
aggrMatrix(fac, w = 1, mean = FALSE, facnames = FALSE)
}
\arguments{
\item{fac}{factor variable.}

\item{w}{vector of weights associated with the levels of \code{fac}.}

\item{mean}{if \code{TRUE}, aggregation will produce (weighted) means instead of sums.}

\item{facnames}{whether the factor levels should be used as column names for the aggregation matrix.}
}
\value{
A sparse aggregation matrix of class \code{tabMatrix}.
}
\description{
Utility function to construct a sparse aggregation matrix from a factor
}
\examples{
n <- 1000
f <- sample(1:100, n, replace=TRUE)
x <- runif(n)
M <- aggrMatrix(f)
all.equal(crossprod_mv(M, x), as.vector(tapply(x, f, sum)))

}
