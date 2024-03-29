% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/samplers.R
\name{create_sampler}
\alias{create_sampler}
\title{Create a sampler object}
\usage{
create_sampler(
  formula,
  data = NULL,
  family = "gaussian",
  ny = NULL,
  ry = NULL,
  r.mod,
  sigma.fixed = NULL,
  sigma.mod = NULL,
  Q0 = NULL,
  formula.V = NULL,
  logJacobian = NULL,
  linpred = NULL,
  compute.weights = FALSE,
  block = NULL,
  prior.only = FALSE,
  control = sampler_control()
)
}
\arguments{
\item{formula}{formula to specify the response variable and additive model components. The model components
form the linear predictor part of the model. A model component on the right hand side can be either
a regression term specified by \code{\link{reg}(...)}, a covariates subject to error term specified
by \code{\link{mec}(...)}, or a generic random effect term specified by \code{\link{gen}(...)}.
See for details the help pages for these model component creation functions.
An offset can be specified as \code{offset(...)}.
Other terms in the formula are collectively interpreted as ordinary regression effects,
treated in the same way as a \code{reg(...)} term, but without the option to change the prior.}

\item{data}{data frame with n rows in which the variables specified in model components can be found.}

\item{family}{character string describing the data distribution. The default is 'gaussian'.
Other options are 'binomial', 'multinomial', 'negbinomial' for the negative binomial distribution,
'poisson', and 'gamma'. See \link{mcmcsae-family} for the related functions
that can be used to specify \code{family} and associated parameters and controls.
For the binomial distribution logistic and probit link functions are supported, the latter only
for binary data. For the negative binomial, Poisson and gamma sampling distributions
a log link function is assumed.
Note that currently \code{family = 'poisson'} is implemented using the negative
binomial distribution with its (reciprocal) overdispersion parameter set
to a very large value.
For categorical or multinomial data, \code{family = "multinomial"} can be used. The implementation
is based on a stick-breaking representation of the multinomial distribution, and the logistic link
function relates each category except the last to a linear predictor. The categories can be
referenced in the model specification formula by 'cat_'.}

\item{ny}{in case \code{family="binomial"} the (vector of) numbers of trials.
It can be either a numeric vector or the name of a variable in \code{data}.
Defaults to a vector of 1s.}

\item{ry}{in case \code{family="negbinomial"} the known, i.e. fixed part of the (reciprocal)
dispersion parameter. It can be specified either as a numeric vector or the name of a
numeric variable in \code{data}. The overall dispersion parameter is the product of \code{ry}
with a positive scalar factor modelled as specified by argument \code{r.mod}. By default
\code{ry} is taken to be 1. For \code{family = "poisson"} a single value can be specified,
determining how well the Poisson distribution is approximated by the negative binomial distribution.
The value should be large enough such that the negative binomial's overdispersion
becomes negligible, but not too large as this might result in slow MCMC mixing. The default is
\code{ry=100} in this case.}

\item{r.mod}{prior specification for a scalar (reciprocal) dispersion parameter
of the negative binomial distribution. The prior can be specified by a call to a prior
specification function. Currently \code{\link{pr_invchisq}}, \code{\link{pr_gig}} and
\code{\link{pr_fixed}} are supported. The default is a chi-squared prior with 1 degree
of freedom. To set the overall dispersion parameter to the value(s) specified by \code{ry},
use \code{r.mod = pr_fixed(value=1)}.}

\item{sigma.fixed}{for Gaussian models, if \code{TRUE} the residual standard deviation parameter 'sigma_' is fixed at 1. In that case
argument \code{sigma.mod} is ignored. This is convenient for Fay-Herriot type models with (sampling) variances assumed to be known.
Default is \code{FALSE}.}

\item{sigma.mod}{prior for the variance parameter of a gaussian sampling distribution.
This can be specified by a call to one of the prior specification functions
\code{\link{pr_invchisq}}, \code{\link{pr_exp}}, \code{\link{pr_gig}} or \code{\link{pr_fixed}} for
inverse chi-squared, exponential, generalized inverse gaussian or degenerate prior distribution,
respectively. The default is an improper prior \code{pr_invchisq(df=0, scale=1)}. A half-t prior on the
standard deviation can be specified using \code{\link{pr_invchisq}} with a chi-squared distributed scale
parameter.}

\item{Q0}{n x n data-level precision matrix for a Gaussian model. It defaults to the unit matrix.
If an n-vector is provided it will be expanded to a (sparse) diagonal matrix with Q0 on its diagonal.
If a name is supplied it will be looked up in \code{data} and subsequently expanded to a diagonal matrix.}

\item{formula.V}{a formula specifying the terms of a variance model in the case of a Gaussian likelihood.
Currently two types of terms are supported: a regression term for the log-variance
specified with \code{\link{vreg}(...)}, and a term \code{\link{vfac}(...)} for multiplicative modeled factors
at a certain level specified by a factor variable. By using unit-level inverse-chi-squared factors the marginal
sampling distribution becomes a Student-t distribution, and by using unit-level exponential factors it becomes
a Laplace or double exponential distribution.}

\item{logJacobian}{if the data are transformed the logarithm of the Jacobian can be supplied so that it
is incorporated in all log-likelihood computations. This can be useful for comparing information criteria
for different transformations. It should be supplied as a vector of the same size as the response variable,
and is currently only supported if \code{family="gaussian"}.
For example, when a log-transformation is used on response vector \code{y}, the vector \code{-log(y)}
should be supplied.}

\item{linpred}{a list of matrices defining (possibly out-of-sample) linear predictors to be simulated.
This allows inference on e.g. (sub)population totals or means. The list must be of the form
\code{list(name_1=X_1, ...)} where the names refer to the model component names and predictions are
computed by summing \code{X_i \%*\% p[[name_i]]}. Alternatively, \code{linpred="fitted"} can be used
as a short-cut for simulations of the full in-sample linear predictor.}

\item{compute.weights}{if \code{TRUE} weights are computed for each element of \code{linpred}. Note that for
a large dataset in combination with vector-valued linear predictors the weights can take up a lot of memory.
By default only means are stored in the simulation carried out using \code{\link{MCMCsim}}.}

\item{block}{DEPRECATED, please use argument \code{control} instead, see also \code{\link{sampler_control}}.
Note that this parameter is now by default set to \code{TRUE}.}

\item{prior.only}{whether a sampler is set up only for sampling from the prior or for sampling from both prior
and posterior distributions. Default \code{FALSE}. If \code{TRUE} there is no need to specify a response in
\code{formula}. This is used by \code{\link{generate_data}}, which samples from the prior predictive
distribution.}

\item{control}{a list with further computational options. These options can
be specified using function \code{\link{sampler_control}}.}
}
\value{
A sampler object, which is the main input for the MCMC simulation
 function \code{\link{MCMCsim}}. The sampler object is an environment with
 precomputed quantities and functions. The main functions are \code{rprior},
 which returns a sample from the prior distributions, \code{draw},
 which returns a sample from the full conditional posterior distributions,
 and \code{start}, which returns a list with starting values for the Gibbs
 sampler. If \code{prior.only} is \code{TRUE}, functions \code{draw} and
 \code{start} are not created.
}
\description{
This function sets up a sampler object, based on the specification of a model. The object contains functions to
draw a set of model parameters from their prior and conditional posterior distributions, and
to generate starting values for the MCMC simulation. The functions share a common environment
containing precomputed quantities such as design matrices based on the model and the data.
The sampler object is the main input for the MCMC simulation function \code{\link{MCMCsim}}.
}
\details{
The right hand side of the \code{formula} argument to \code{create_sampler} can be used to specify
additive model components. Currently four model components are supported: \code{\link{reg}(...)}
for regression or 'fixed' effects, \code{\link{gen}(...)} for generic random effects,
\code{\link{mec}(...)} for measurement in covariates effects, and \code{\link{brt}(...)}
for a Bayesian additive regression trees component. Note that an offset can be added
separately, in the usual way using \code{\link{offset}(...)}.

For gaussian models, \code{formula.V} can be used to specify the variance structure of the model.
Currently two specialized variance model components are supported, \code{\link{vreg}(...)} for regression
effects predicting the log-variance and \code{\link{vfac}(...)} for modeled variance factors.
}
\examples{
# first generate some data
n <- 200
x <- rnorm(n)
y <- 0.5 + 2*x + 0.3*rnorm(n)
# create a sampler for a simple linear regression model
sampler <- create_sampler(y ~ x)
sim <- MCMCsim(sampler)
(summary(sim))

y <- rbinom(n, 1, 1 / (1 + exp(-(0.5 + 2*x))))
# create a sampler for a binary logistic regression model
sampler <- create_sampler(y ~ x, family="binomial")
sim <- MCMCsim(sampler)
(summary(sim))

}
\references{
J.H. Albert and S. Chib (1993).
   Bayesian analysis of binary and polychotomous response data.
   Journal of the American statistical Association 88(422), 669-679.

 D. Bates, M. Maechler, B. Bolker and S.C. Walker (2015).
   Fitting Linear Mixed-Effects Models Using lme4.
   Journal of Statistical Software 67(1), 1-48.

 S.W. Linderman, M.J. Johnson and R.P. Adams (2015).
   Dependent multinomial models made easy: Stick-breaking with the Polya-Gamma augmentation.
   Advances in Neural Information Processing Systems, 3456-3464.

 P.A. Parker, S.H. Holan and R. Janicki (2023).
   Conjugate Modeling Approaches for Small Area Estimation with Heteroscedastic Structure.
   Journal of Survey Statistics and Methodology, smad002.

 N. Polson, J.G. Scott and J. Windle (2013).
   Bayesian Inference for Logistic Models Using Polya-Gamma Latent Variables.
   Journal of the American Statistical Association 108(504), 1339-1349.

 H. Rue and L. Held (2005).
   Gaussian Markov Random Fields.
   Chapman & Hall/CRC.

 M. Zhou and L. Carin (2015).
   Negative Binomial Process Count and Mixture Modeling.
   IEEE Transactions on Pattern Analysis and Machine Intelligence 37(2), 307-320.
}
