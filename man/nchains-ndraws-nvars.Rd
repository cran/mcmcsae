% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MCMCsim.R
\name{nchains-ndraws-nvars}
\alias{nchains-ndraws-nvars}
\alias{nchains}
\alias{ndraws}
\alias{nvars}
\title{Get the number of chains, samples per chain or the number of variables in a simulation object}
\usage{
nchains(obj)

ndraws(obj)

nvars(dc)
}
\arguments{
\item{obj}{an mcdraws object or a draws component (dc) object.}

\item{dc}{a draws component object.}
}
\value{
The number of chains or retained samples per chain or
 the number of variables.
}
\description{
Get the number of chains, samples per chain or the number of variables in a simulation object
}
\examples{
\donttest{
ex <- mcmcsae_example(n=50)
sampler <- create_sampler(ex$model, data=ex$dat)
sim <- MCMCsim(sampler, burnin=100, n.iter=300, thin=2, n.chain=5, store.all=TRUE)
# resolve possible conflict with posterior package:
nchains <- mcmcsae::nchains; ndraws <- mcmcsae::ndraws
nchains(sim); nchains(sim$beta)
ndraws(sim); ndraws(sim$beta)
nvars(sim$beta); nvars(sim$sigma_); nvars(sim$llh_); nvars(sim$v)
plot(sim, "beta")
nchains(subset(sim$beta, chains=1:2))
ndraws(subset(sim$beta, draws=sample(1:ndraws(sim), 100)))
nvars(subset(sim$u, vars=1:2))
}

}
