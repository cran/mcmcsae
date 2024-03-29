% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MCMCsim.R
\name{print.mcdraws_summary}
\alias{print.mcdraws_summary}
\title{Print a summary of MCMC simulation results}
\usage{
\method{print}{mcdraws_summary}(x, digits = 3L, max.lines = 10L, tail = FALSE, sort = NULL, ...)
}
\arguments{
\item{x}{an object of class \code{mcdraws_summary} as output by \code{\link{summary.mcdraws}}.}

\item{digits}{number of digits to use, defaults to 3.}

\item{max.lines}{maximum number of elements per vector parameter to display.
If \code{NULL}, all elements are displayed.}

\item{tail}{if \code{TRUE} the last instead of first \code{max.lines} of each component
are displayed.}

\item{sort}{column name on which to sort the output.}

\item{...}{passed on to \code{print.default}.}
}
\description{
Display a summary of an \code{mcdraws} object, as output by \code{\link{MCMCsim}}.
}
\examples{
\donttest{
ex <- mcmcsae_example()
sampler <- create_sampler(ex$model, data=ex$dat)
sim <- MCMCsim(sampler, store.all=TRUE)
print(summary(sim), sort="n_eff")
}

}
