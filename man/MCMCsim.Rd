% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MCMCsim.R
\name{MCMCsim}
\alias{MCMCsim}
\title{Run a Markov Chain Monte Carlo simulation}
\usage{
MCMCsim(
  sampler,
  from.prior = FALSE,
  n.iter = 1000L,
  n.chain = 3L,
  thin = 1L,
  burnin = if (from.prior) 0L else 250L,
  start = NULL,
  store,
  store.all = FALSE,
  pred = NULL,
  store.mean,
  store.sds = FALSE,
  to.file = NULL,
  filename = "MCdraws_",
  write.single.prec = FALSE,
  verbose = TRUE,
  n.progress = n.iter\%/\%10L,
  trace.convergence = NULL,
  stop.on.convergence = FALSE,
  convergence.bound = 1.05,
  plot.trace = NULL,
  add.to.plot = TRUE,
  plot.type = "l",
  n.cores = 1L,
  cl = NULL,
  seed = NULL,
  export = NULL
)
}
\arguments{
\item{sampler}{sampler object created by \code{\link{create_sampler}}.}

\item{from.prior}{whether to sample from the prior. By default \code{from.prior=FALSE}
and samples are taken from the posterior.}

\item{n.iter}{number of draws after burnin.}

\item{n.chain}{number of independent chains.}

\item{thin}{only every \code{thin}'th draw is kept.}

\item{burnin}{number of draws to discard at the beginning of each chain.}

\item{start}{an optional function to generate starting values or a list containing for each chain
a named list of starting values. It may be used to provide starting values for some or all parameters.
The sampler object's own start function, if it exists, is called to generate any starting values not
provided by the user.}

\item{store}{vector of names of parameters to store MCMC draws for. By default, simulations are
stored for all parameters returned by \code{sampler$store_default}.}

\item{store.all}{if \code{TRUE} simulation vectors of all parameters returned by the sampling
function of \code{sampler} will be stored. The default is \code{FALSE}, and in that case
only simulations for the parameters named in \code{store} are stored.}

\item{pred}{list of character strings defining derived quantities to be computed (and stored) for each draw.}

\item{store.mean}{vector of names of parameters for which only the mean (per chain) is to be stored.
This may be useful for large vector parameters (e.g. regression residuals) for which storing complete
MCMC output would use too much memory. The function \code{sampler$store_mean_default}
exists it provides the default.}

\item{store.sds}{if \code{TRUE} store for all parameters in \code{store.mean}, besides the mean, also
the standard deviation. Default is \code{FALSE}.}

\item{to.file}{vector of names of parameters to write to file.}

\item{filename}{name of file to write parameter draws to.
Each named parameter is written to a separate file, named \code{filename_parametername}.}

\item{write.single.prec}{Whether to write to file in single precision. Default is \code{FALSE}.}

\item{verbose}{if \code{FALSE} no output is sent to the screen during the simulation. \code{TRUE} by default.}

\item{n.progress}{update diagnostics and plots after so many iterations.}

\item{trace.convergence}{vector of names of parameters for which Gelman-Rubin R-hat diagnostics are printed to the screen every \code{n.progress} iterations.}

\item{stop.on.convergence}{if \code{TRUE} stop the simulation if the R-hat diagnostics for all parameters in \code{trace.convergence} are less than \code{convergence.bound}.}

\item{convergence.bound}{threshold used with \code{stop.on.convergence}.}

\item{plot.trace}{character vector of parameter names for which to plot draws
during the simulation. For one or two parameters trace plots will be shown,
and if more parameters are specified the results will be displayed in a pairs
plot. For vector parameters a specific component can be selected using brackets,
e.g. \code{"beta[2]"}.}

\item{add.to.plot}{if \code{TRUE} the plot is updated every \code{n.progress} iterations,
otherwise a new plot (with new scales) is created after every \code{n.progress} iterations.}

\item{plot.type}{default is "l" (lines).}

\item{n.cores}{the number of cpu cores to use. Default is 1, i.e. no parallel computation.
If an existing cluster \code{cl} is provided, \code{n.cores} will be set to the number
of workers in that cluster.}

\item{cl}{an existing cluster can be passed for parallel computation. If \code{NULL} and
\code{n.cores > 1}, a new cluster is created.}

\item{seed}{a random seed (integer). For parallel computation it is used to independently
seed RNG streams for all workers.}

\item{export}{a character vector with names of objects to export to the workers. This may
be needed for parallel execution if expressions in \code{pred} depend on global variables.}
}
\value{
An object of class \code{mcdraws} containing posterior draws as well as some meta information.
}
\description{
Given a sampler object this function runs a MCMC simulation and stores the
posterior draws. A sampler object for a wide class of multilevel models
can be created using \code{\link{create_sampler}}, but users can also define
their own sampler functions, see below.
\code{MCMCsim} allows to choose the parameters for which simulation results
must be stored. It is possible to define derived quantities that will also
be stored. To save memory, it is also possible to only store Monte Carlo
means/standard errors for some large vector parameters, say. Another
way to use less memory is to save the simulation results of large vector
parameters to file.
For parameters specified in \code{plot.trace} trace plots or pair plots of
multiple parameters are displayed during the simulation.
}
\details{
A sampler object is an environment containing data and functions to use
for sampling. The following elements of the sampler object are used by
\code{MCMCsim}:
\describe{
  \item{start}{function to generate starting values.}
  \item{draw}{function to draw samples, typically from a full conditional
    posterior distribution.}
  \item{rprior}{function to draw from a prior distribution.}
  \item{coef.names}{list of vectors of parameter coefficient names, for
    vector parameters.}
  \item{MHpars}{vector of names of parameters that are sampled using a
    Metropolis-Hastings (MH) sampler; acceptance rates are kept for these
    parameters.}
  \item{adapt}{function of acceptance rates of \code{MHpars} to adapt
    MH-kernel, called every 100 iterations during the burn-in period.}
}
}
\examples{
# 1. create a sampler function
sampler <- new.env()
sampler$draw <- function(p) list(x=rnorm(1L), y=runif(1L))
# 2. do the simulation
sim <- MCMCsim(sampler, store=c("x", "y"))
str(sim)
summary(sim)

# example that requires start values or a start function
sampler$draw <- function(p) list(x=rnorm(1L), y=p$x * runif(1L))
sampler$start <- function(p) list(x=rnorm(1L), y=runif(1L))
sim <- MCMCsim(sampler, store=c("x", "y"))
summary(sim)
plot(sim, c("x", "y"))

# example using create_sampler; first generate some data
n <- 100
dat <- data.frame(x=runif(n), f=as.factor(sample(1:4, n, replace=TRUE)))
gd <- generate_data(~ reg(~ x + f, prior=pr_normal(precision=1), name="beta"), data=dat)
dat$y <- gd$y
sampler <- create_sampler(y ~ x + f, data=dat)
sim <- MCMCsim(sampler, burnin=100, n.iter=400, n.chain=2)
(summary(sim))
gd$pars

}
