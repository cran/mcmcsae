% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/f_gamma.R, R/f_poisson.R, R/family.R
\name{f_gamma}
\alias{f_gamma}
\alias{f_gaussian_gamma}
\alias{f_poisson}
\alias{mcmcsae-family}
\alias{f_gaussian}
\alias{f_binomial}
\alias{f_negbinomial}
\alias{f_multinomial}
\title{Functions for specifying a sampling distribution and link function}
\usage{
f_gamma(
  link = "log",
  shape.vec = ~1,
  shape.prior = pr_gamma(0.1, 0.1),
  control = set_MH(type = "RWLN", scale = 0.2, adaptive = TRUE)
)

f_gaussian_gamma(link = "identity", var.data, ...)

f_poisson(link = "log", size = 100)

f_gaussian(link = "identity")

f_binomial(link = c("logit", "probit"))

f_negbinomial(link = "logit")

f_multinomial(link = "logit", K = NULL)
}
\arguments{
\item{link}{the name of a link function. Currently the only allowed link functions are:
\code{"identity"} for (log-)Gaussian sampling distributions, \code{"logit"} (default) and \code{"probit"}
for binomial distributions and \code{"log"} for negative binomial sampling distributions.}

\item{shape.vec}{optional formula specification of unequal shape parameter for gamma family}

\item{shape.prior}{prior for gamma shape parameter. Supported prior distributions:
\code{\link{pr_fixed}} with a default value of 1, \code{\link{pr_exp}} and
\code{\link{pr_gamma}}. The current default is \code{pr_gamma(shape=0.1, rate=0.1)}.}

\item{control}{options for the Metropolis-Hastings algorithm employed
in case the shape parameter is to be inferred. Function \code{\link{set_MH}}
can be used to change the default options. The two choices of proposal
distribution type supported are "RWLN" for a random walk proposal on the
log-shape scale, and "gamma" for an approximating gamma proposal, found using
an iterative algorithm. In the latter case, a Metropolis-Hastings accept-reject
step is currently omitted, so the sampling algorithm is an approximate one,
though often quite accurate and efficient.}

\item{var.data}{the (variance) data for the gamma part of family \code{gaussian_gamma}.}

\item{...}{further arguments passed to \code{f_gamma}.}

\item{size}{size or dispersion parameter of the negative binomial distribution used
internally to approximate the Poisson distribution. This should be set to a relatively
large value (default is 100), corresponding to negligible overdispersion, to obtain a
good approximation. However, too large values may cause slow MCMC exploration of the
posterior distribution.}

\item{K}{number of categories for multinomial model; this must be specified for prior predictive sampling.}
}
\value{
A family object.
}
\description{
These functions are intended for use in the \code{family} argument of \code{\link{create_sampler}}.
In future versions these functions may gain additional arguments, but currently the corresponding
functions \code{gaussian} and \code{binomial} can be used as well.
}
\references{
J.W. Miller (2019).
   Fast and Accurate Approximation of the Full Conditional for Gamma Shape Parameters.
   Journal of Computational and Graphical Statistics 28(2), 476-480.
}
