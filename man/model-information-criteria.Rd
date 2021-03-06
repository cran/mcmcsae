% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_eval.R
\encoding{UTF-8}
\name{model-information-criteria}
\alias{model-information-criteria}
\alias{compute_DIC}
\alias{compute_WAIC}
\alias{waic.draws}
\alias{loo.draws}
\title{Compute DIC, WAIC and leave-one-out cross-validation model measures}
\usage{
compute_DIC(x, use.pV = FALSE)

compute_WAIC(x, diagnostic = FALSE, batch.size = NULL, show.progress = TRUE)

\method{waic}{draws}(x, by.unit = FALSE, ...)

\method{loo}{draws}(x, r_eff = FALSE, n.cores = 1L, ...)
}
\arguments{
\item{x}{an object of class \code{draws}.}

\item{use.pV}{whether half the posterior variance of the deviance should be used
as an alternative estimate of the effective number of model parameters for DIC.}

\item{diagnostic}{whether vectors of log-pointwise-predictive-densities and pointwise
contributions to the WAIC effective number of model parameters should be returned.}

\item{batch.size}{number of data units to process per batch.}

\item{show.progress}{whether to show a progress bar.}

\item{by.unit}{if \code{TRUE} the computation is carried out unit-by-unit, which is
slower but uses much less memory.}

\item{...}{Other arguments, passed to \code{\link[loo]{loo}}. Not currently
used by \code{waic.draws}.}

\item{r_eff}{whether to compute relative effective sample size estimates
for the likelihood of each observation. This takes more time, but should
result in a better PSIS approximation. See \code{\link[loo]{loo}}.}

\item{n.cores}{how many cores to use.}
}
\value{
For \code{compute_DIC} a vector with the deviance information criterion and
 effective number of model parameters. For \code{compute_WAIC} a vector with the
 WAIC model selection criterion and WAIC effective number of model parameters.
 Method \code{waic} returns an object of class \code{waic, loo}, see the
 documentation for \code{\link[loo]{waic}} in package \pkg{loo}. 
 Method \code{loo} returns an object of class \code{psis_loo}, see
 \code{\link[loo]{loo}}.
}
\description{
Compute the Deviance Information Criterion (DIC) or
Watanabe-Akaike Information Criterion (WAIC) from an
object of class \code{draws} output by \code{\link{MCMCsim}}.
Method \code{waic.draws} computes WAIC using package \pkg{loo}.
Method \code{loo.draws} also depends on package \pkg{loo} to compute
a Pareto-smoothed importance sampling (PSIS) approximation
to leave-one-out cross-validation.
}
\examples{
\donttest{
ex <- mcmcsae_example(n=100)
sampler <- create_sampler(ex$model, data=ex$dat)
sim <- MCMCsim(sampler, burnin=100, n.iter=300, n.chain=4, store.all=TRUE)
compute_DIC(sim)
compute_WAIC(sim)
if (require(loo)) {
  waic(sim)
  loo(sim, r_eff=TRUE)
}
}

}
\references{
D. Spiegelhalter, N. Best, B. Carlin and A. van der Linde (2002).
   Bayesian Measures of Model Complexity and Fit.
   Journal of the Royal Statistical Society B 64 (4), 583-639.

 S. Watanabe (2010).
   Asymptotic equivalence of Bayes cross validation and widely applicable
   information criterion in singular learning theory.
   Journal of Machine Learning 11, 3571-3594.

 A. Gelman, J. Hwang and A. Vehtari (2014).
   Understanding predictive information criteria for Bayesian models.
   Statistics and Computing 24, 997-1016.

 A. Vehtari, A. Gelman and J. Gabry (2015).
   Pareto smoothed importance sampling.
   arXiv preprint arXiv:1507.02646.

 A. Vehtari, A. Gelman and J. Gabry (2017).
   Practical Bayesian model evaluation using leave-one-out cross-validation and WAIC.
   Statistics and Computing 27, 1413-1432.

 P.-C. Buerkner, J. Gabry and A. Vehtari (2019).
   Bayesian leave-one-out cross-validation for non-factorizable normal models.
   arXiv:1810.10559v3.
}
