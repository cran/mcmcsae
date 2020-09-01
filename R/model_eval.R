#' Extract draws of fitted values or residuals from a draws object
#'
#' For a model created with \code{\link{create_sampler}} and estimated using \code{\link{MCMCsim}},
#' these functions return the posterior draws of fitted values or residuals.
#' In the current implementation the fitted values correspond to the linear predictor
#' and the residuals are computed as the data vector minus the fitted values,
#' regardless of the model's distribution family.
#' For large datasets the returned object can become very large. One may therefore
#' select a subset of draws or chains or use \code{mean.only=TRUE} to
#' return a vector of posterior means only.
#'
#' @examples
#' \donttest{
#' ex <- mcmcsae_example(n=50)
#' sampler <- create_sampler(ex$model, data=ex$dat)
#' sim <- MCMCsim(sampler, burnin=100, n.iter=300, thin=2, store.all=TRUE)
#' fitted(sim, mean.only=TRUE)
#' summary(fitted(sim))
#' residuals(sim, mean.only=TRUE)
#' summary(residuals(sim))
#' bayesplot::mcmc_intervals(as.matrix(subset(residuals(sim), vars=1:20)))
#' }
#'
#' @param object a draws object.
#' @param mean.only if \code{TRUE} only the vector of posterior means is returned. In that case
#'  the subsequent arguments are ignored. Default is \code{FALSE}.
#' @param units the data units (by default all) for which fitted values or residuals should be
#'  computed.
#' @param chains optionally, a selection of chains.
#' @param draws optionally, a selection of draws per chain.
#' @param matrix whether a matrix should be returned instead of a dc object.
#' @param type the type of fitted values: "link" for fitted values on the linear predictor scale
#'  (the default), and "response" for fitted values on the response scale. Returned residuals are
#'  always on the response scale.
#' @param ... currently not used.
#' @return Either a draws component object or a matrix with draws of fitted values or residuals.
#'  The residuals are always on the response scale, whereas fitted values can
#'  be on the scale of the linear predictor or the response depending on \code{type}.
#'  If \code{mean.only=TRUE}, a vector of posterior means.
#' @name residuals-fitted-values
NULL

#' @export
## @method fitted draws
#' @rdname residuals-fitted-values
fitted.draws <- function(object, mean.only=FALSE, units=NULL, chains=seq_len(nchains(object)),
                         draws=seq_len(ndraws(object)), matrix=FALSE, type=c("link", "response"), ...) {
  type <- match.arg(type)
  fitted_res(object, mean.only, units, chains, draws, matrix, type, resid=FALSE)
}

#' @export
## @method residuals draws
#' @rdname residuals-fitted-values
residuals.draws <- function(object, mean.only=FALSE, units=NULL, chains=seq_len(nchains(object)),
                            draws=seq_len(ndraws(object)), matrix=FALSE, ...) {
  fitted_res(object, mean.only, units, chains, draws, matrix, type="response", resid=TRUE)
}

# get fitted values or residuals, internal function
fitted_res <- function(object, mean.only=FALSE, units=NULL, chains=seq_len(nchains(object)),
                       draws=seq_len(ndraws(object)), matrix=FALSE, type, resid, ...) {
  smplr <- object[["_model"]]
  if (mean.only && (!is.null(object[["e_"]]) || !is.null(object[["_means"]][["e_"]]))) {
    if ((resid && smplr$e.is.res) || (!resid && !smplr$e.is.res))
      return(get_means(object, "e_")[[1L]])
    else
      return(smplr$y - get_means(object, "e_")[[1L]])
  }
  if (is.null(units)) {
    units <- seq_len(smplr$n)
    all.units <- TRUE
  } else {
    all.units <- FALSE
  }
  if (is.null(object[["e_"]])) {
    mc <- smplr$mod[[1L]]
    if (all.units) Xi <- mc$X else Xi <- mc$X[units, , drop=FALSE]
    if (matrix) {
      out <- tcrossprod(as.matrix.dc(get_from(object[[mc$name]], chains=chains, draws=draws), colnames=FALSE), Xi)
    } else {
      out <- list()
      for (ch in seq_along(chains))
        out[[ch]] <- tcrossprod(object[[mc$name]][[chains[ch]]][draws, , drop=FALSE], Xi)
    }
    for (mc in smplr$mod[-1L]) {
      if (all.units) Xi <- mc$X else Xi <- mc$X[units, , drop=FALSE]
      if (matrix)
        out <- out + tcrossprod(as.matrix.dc(get_from(object[[mc$name]], chains=chains, draws=draws), colnames=FALSE), Xi)
      else
        for (ch in seq_along(chains))
          out[[ch]] <- out[[ch]] + tcrossprod(object[[mc$name]][[chains[ch]]][draws, , drop=FALSE], Xi)
    }
    if (smplr$use.offset) {
      if (matrix)
        out <- out + rep_each(smplr$offset[units], nrow(out))
      else
        out <- lapply(out, function(x) {x + rep_each(smplr$offset[units], nrow(x))})
    }
  } else {
    out <- get_from(object[["e_"]], chains=chains, draws=draws, vars=units)
    if (smplr$e.is.res) {
      # change to linear predictor as code below assumes that
      out <- lapply(out, function(x) rep_each(smplr$y[units], nrow(x)) - x)
    }
  }
  if (type == "response") {
    if (matrix)
      out <- smplr$f_mean(out)
    else
      out <- lapply(out, smplr$f_mean)
  }
  if (resid) {  # compute residuals at the response scale
    if (matrix)
      out <- rep_each(smplr$y[units], nrow(out)) - out
    else
      out <- lapply(out, function(x) rep_each(smplr$y[units], nrow(x)) - x)
  }
  if (matrix) {
    dimnames(out) <- list(NULL, units)
    if (mean.only) out <- colMeans(out)
  } else {
    attr(out, "labels") <- units
    class(out) <- "dc"
    if (mean.only) out <- colMeans(as.matrix(out))
  }
  out
}

#' Extract weights from a draws object
#'
#' @examples
#' \donttest{
#' # first create a population data frame
#' N <- 1000  # population size
#' pop <- data.frame(x=rnorm(N), area=factor(sample(1:10, N, replace=TRUE)))
#' pop$y <- 1 + 2*pop$x + seq(-1, to=1, length.out=10)[pop$area] + 0.5*rnorm(N)
#' pop$sample <- FALSE
#' pop$sample[sample(seq_len(N), 100)] <- TRUE
#' # a simple linear regression model:
#' sampler <- create_sampler(
#'   y ~ reg(~ x, name="beta"),
#'   linpred=list(beta=rowsum(model.matrix(~ x, pop), pop$area)), compute.weights=TRUE,
#'   data=pop[pop$sample, ]
#' )
#' sim <- MCMCsim(sampler)
#' (summary(sim))
#' str(weights(sim))
#' crossprod_mv(weights(sim), pop$y[pop$sample])
#' summary(sim$linpred_)
#' # a multilevel model:
#' sampler <- create_sampler(
#'   y ~ reg(~ x, name="beta") + gen(factor = ~ area, name="v"),
#'   linpred=list(beta=rowsum(model.matrix(~ x, pop), pop$area), v=diag(10)), compute.weights=TRUE,
#'   data=pop[pop$sample, ]
#' )
#' sim <- MCMCsim(sampler)
#' (summary(sim))
#' str(weights(sim))
#' crossprod_mv(weights(sim), pop$y[pop$sample])
#' summary(sim$linpred_)
#' }
#'
#' @export
## @method weights draws
#' @param object an object of class \code{draws}.
#' @param ... currently not used.
#' @return A vector with (simulation means of) weights.
weights.draws <- function(object, ...) get_means(object, "weights_")[[1L]]


#' Compute DIC, WAIC and leave-one-out cross-validation model measures
#'
#' Compute the Deviance Information Criterion (DIC) or
#' Watanabe-Akaike Information Criterion (WAIC) from an
#' object of class \code{draws} output by \code{\link{MCMCsim}}.
#' Method \code{waic.draws} computes WAIC using package \pkg{loo}.
#' Method \code{loo.draws} also depends on package \pkg{loo} to compute
#' a Pareto-smoothed importance sampling (PSIS) approximation
#' to leave-one-out cross-validation.
#'
#' @examples
#' \donttest{
#' ex <- mcmcsae_example(n=100)
#' sampler <- create_sampler(ex$model, data=ex$dat)
#' sim <- MCMCsim(sampler, burnin=100, n.iter=300, n.chain=4, store.all=TRUE)
#' compute_DIC(sim)
#' compute_WAIC(sim)
#' if (require(loo)) {
#'   waic(sim)
#'   loo(sim, r_eff=TRUE)
#' }
#' }
#'
#' @param x an object of class \code{draws}.
#' @param use.pV whether half the posterior variance of the deviance should be used
#'  as an alternative estimate of the effective number of model parameters for DIC.
#' @param diagnostic whether vectors of log-pointwise-predictive-densities and pointwise
#'  contributions to the WAIC effective number of model parameters should be returned.
#' @param batch.size number of data units to process per batch.
#' @param show.progress whether to show a progress bar.
#' @param by.unit if \code{TRUE} the computation is carried out unit-by-unit, which is
#'  slower but uses much less memory.
#' @param r_eff whether to compute relative effective sample size estimates
#'  for the likelihood of each observation. This takes more time, but should
#'  result in a better PSIS approximation. See \code{\link[loo]{loo}}.
#' @param n.cores how many cores to use.
#' @param ... Other arguments, passed to \code{\link[loo]{loo}}. Not currently
#'  used by \code{waic.draws}.
#' @return For \code{compute_DIC} a vector with the deviance information criterion and
#'  effective number of model parameters. For \code{compute_WAIC} a vector with the
#'  WAIC model selection criterion and WAIC effective number of model parameters.
#'  Method \code{waic} returns an object of class \code{waic, loo}, see the
#'  documentation for \code{\link[loo]{waic}} in package \pkg{loo}. 
#'  Method \code{loo} returns an object of class \code{psis_loo}, see
#'  \code{\link[loo]{loo}}.
#' @encoding UTF-8
#' @references
#'  D. Spiegelhalter, N. Best, B. Carlin and A. van der Linde (2002).
#'    Bayesian Measures of Model Complexity and Fit.
#'    Journal of the Royal Statistical Society B 64 (4), 583-639.
#'
#'  S. Watanabe (2010).
#'    Asymptotic equivalence of Bayes cross validation and widely applicable
#'    information criterion in singular learning theory.
#'    Journal of Machine Learning 11, 3571-3594.
#'
#'  A. Gelman, J. Hwang and A. Vehtari (2014).
#'    Understanding predictive information criteria for Bayesian models.
#'    Statistics and Computing 24, 997-1016.
#'
#'  A. Vehtari, A. Gelman and J. Gabry (2015).
#'    Pareto smoothed importance sampling.
#'    arXiv preprint arXiv:1507.02646.
#'
#'  A. Vehtari, A. Gelman and J. Gabry (2017).
#'    Practical Bayesian model evaluation using leave-one-out cross-validation and WAIC.
#'    Statistics and Computing 27, 1413-1432.
#'
#'  P.-C. Buerkner, J. Gabry and A. Vehtari (2019).
#'    Bayesian leave-one-out cross-validation for non-factorizable normal models.
#'    arXiv:1810.10559v3.
#' @name model-information-criteria
NULL

#' @export
#' @rdname model-information-criteria
compute_DIC <- function(x, use.pV=FALSE) {
  if (!inherits(x, "draws")) stop("not an object of class 'draws'")
  if (x[["_info"]]$from.prior) stop("cannot compute DIC for draws from prior")
  post.means <- get_means(x)
  if (is.null(post.means[["e_"]]))
    stop("cannot compute DIC: missing simulation means of ", if (x[["_model"]]$e.is.res) "residuals" else "linear predictor", "'e_'")
  if (x[["_model"]]$modeled.Q && (x[["_model"]]$family$family == "gaussian")) {
    if (is.null(post.means[["Q_"]])) stop("cannot compute DIC: missing simulation means of 'Q_'")
    if (x[["_model"]]$Q0.type == "symm") {
      # reconstruct mean sparse precision matrix
      post.means[["QM_"]] <- block_scale_dsCMatrix(x[["_model"]]$Q0, post.means[["Q_"]])
    }
  }
  if (is.null(post.means[["llh_"]]))
    stop("cannot compute DIC; missing simulation means of 'llh_'")
  deviance_at_mean <- -2 * x[["_model"]]$llh(post.means)
  mean_deviance <- -2 * post.means[["llh_"]]
  if (use.pV) {
    # larger Monte Carlo error(?), but guaranteed to be positive
    pDIC <- 2 * c(var(as.matrix.dc(x[["llh_"]], colnames=FALSE)))
  } else {  # default
    pDIC <- mean_deviance - deviance_at_mean
  }
  c(DIC=deviance_at_mean + 2 * pDIC, p_DIC=pDIC)
}

get_lppd_function <- function(x) {
  llh_i <- x[["_model"]]$llh_i
  if (is.null(llh_i)) stop("pointwise log-likelihood function not implemented; cannot compute WAIC")
  if (!("e_" %in% par_names(x) || all(names(x[["_model"]]$mod) %in% par_names(x))))
    stop("WAIC can only be computed if all coefficients are stored. Please use 'store.all=TRUE' in MCMCsim.")
  if (x[["_model"]]$modeled.Q && x[["_model"]]$family$family == "gaussian" && !all(names(x[["_model"]]$Vmod) %in% par_names(x)))
    stop("WAIC can only be computed if all modeled variance factors are stored. Please use 'store.all=TRUE' in MCMCsim.")
  llh_i
}

#' @export
#' @rdname model-information-criteria
compute_WAIC <- function(x, diagnostic=FALSE, batch.size=NULL, show.progress=TRUE) {
  if (!inherits(x, "draws")) stop("not an object of class 'draws'")
  if (x[["_info"]]$from.prior) stop("cannot compute WAIC for draws from prior")
  n <- x[["_model"]]$n
  if (is.null(batch.size)) {
    # determine batch size to avoid too much memory use; now by default truncated to ~ 80MB chunks
    batch.size <- min(n, 10000000L %/% (nchains(x) * ndraws(x)))
  }
  if ((x[["_model"]]$Q0.type == "symm") && (batch.size != n)) {
    warning("For non-diagonal precision matrix, use 'batch.size' equal to the number of
      observations for correct leave-one-out predictive densities", immediate.=TRUE)
  }
  batch <- seq_len(batch.size)
  n.batch <- n %/% batch.size + (n %% batch.size > 0L)
  llh_i <- get_lppd_function(x)
  lppd_sum <- pWAIC1 <- pWAIC2 <- 0
  if (diagnostic) lppd <- pwaic <- rep.int(NA_real_, n)
  show.progress <- show.progress && n.batch > 2L
  if (show.progress) pb <- txtProgressBar(min=1L, max=n.batch, style=3L)
  for (i in seq_len(n.batch)) {
    if ((i == n.batch) && (n %% batch.size))
      ind <- (n - (n %% batch.size) + 1L):n
    else
      ind <- batch + (i - 1L) * batch.size
    logprob_i <- llh_i(x, ind)
    nr <- dim(logprob_i)[1L]; nc <- dim(logprob_i)[2L]
    lppd_i <- colLogSumExps(logprob_i) - log(nr)
    pwaic_i <- colVars(logprob_i)
    if (diagnostic) {
      lppd[ind] <- lppd_i
      pwaic[ind] <- pwaic_i
    }
    lppd_sum <- lppd_sum + sum(lppd_i)
    pWAIC1 <- pWAIC1 + sum(.colMeans(logprob_i, nr, nc))
    pWAIC2 <- pWAIC2 + sum(pwaic_i)
    if (show.progress) setTxtProgressBar(pb, i)
  }
  if (show.progress) close(pb)
  pWAIC1 <- 2 * (lppd_sum - pWAIC1)
  WAIC1 <- -2 * lppd_sum + 2 * pWAIC1
  WAIC2 <- -2 * lppd_sum + 2 * pWAIC2
  if (diagnostic)
    list(lppd=lppd, pwaic=pwaic, WAIC1=WAIC1, p_WAIC1=pWAIC1, WAIC2=WAIC2, p_WAIC2=pWAIC2)
  else
    c(WAIC1=WAIC1, p_WAIC1=pWAIC1, WAIC2=WAIC2, p_WAIC2=pWAIC2)
}

#' @method waic draws
#' @export
#' @rdname model-information-criteria
#' @importFrom loo waic
# based on waic.stanreg from package rstanarm
waic.draws <- function(x, by.unit=FALSE, ...) {
  if ((x[["_model"]]$Q0.type == "symm") && by.unit) {
    warning("For non-diagonal precision matrix, use 'batch.size' equal to the number of
      observations for correct leave-one-out predictive densities", immediate.=TRUE)
  }
  llh_i <- get_lppd_function(x)
  if (by.unit) {
    data <- data.frame(i=seq_len(x[["_model"]]$n))
    f <- function(data_i, draws) llh_i(draws, data_i$i)
    loo::waic.function(f, data=data, draws=x)
  } else {
    loo::waic.matrix(as.matrix(llh_i(x)))
  }
}

#' @method loo draws
#' @export
#' @rdname model-information-criteria
#' @importFrom loo loo
# TODO
#  - for large data and number of draws, use loo.function
#    and batch processing as in compute_WAIC
loo.draws <- function(x, r_eff=FALSE, n.cores=1L, ...) {
  llh_i <- get_lppd_function(x)
  llh <- llh_i(x)
  if (r_eff)
    r_eff <- loo::relative_eff(exp(llh), chain_id=rep_each(seq_len(nchains(x)), ndraws(x)), cores=n.cores)
  else
    r_eff <- NULL
  loo::loo(llh, r_eff=r_eff, cores=n.cores, ...)
}
