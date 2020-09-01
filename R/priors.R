
# Q0: precision matrix of dimension q
set_prior_precision <- function(Q0=NULL, q, sparse=NULL) {
  if (is.null(Q0)) {
    Q0 <- Cdiag(rep.int(0, q))  # default: noninformative, improper prior
  } else {
    if (is.vector(Q0)) {
      if (length(Q0) == 1L) {
        Q0 <- Cdiag(rep.int(Q0, q))
      } else {
        if (length(Q0) == q)
          Q0 <- Cdiag(Q0)
        else
          stop("wrong input for precision matrix 'Q0'")
      }
    } else {
      if (!identical(dim(Q0), c(q, q))) stop("wrong input for precision matrix 'Q0'")
    }
    Q0 <- economizeMatrix(Q0, symmetric=TRUE, sparse=sparse)
  }
  if (isUnitDiag(Q0)) Q0 <- expandUnitDiag(Q0)  # block sampler expects x-slot
  Q0
}

#' Create an object containing information about exponential prior distributions
#'
#' @export
#' @param scale scale parameter of length 1 or \code{n}.
#' @param n dimension, if known. For internal use only.
#' @param post whether conditional posterior sampling function should be
#'  created. For internal use only.
#' @return An environment with information about the prior and possibly conditional
#'  posterior distribution(s), to be used by other package functions.
# TODO sample precision instead of variance parameters; the default value
#      for p then leads to inverse Gaussian posterior, which may be faster
#      to sample from
pr_exp <- function(scale=1, n=NULL, post=FALSE) {
  if (!all(scale > 0)) stop("scale parameter must be positive")
  if (!is.null(n) && !(length(scale) %in% c(1L, n))) stop("scale parameter has wrong length")
  if (!is.null(n)) {
    n <- as.integer(n)
    rprior <- function() scale * rexp(n)
    if (post) {
      # TODO set a=2/scale here, and default p=-1/2
      # b is typically updated
      draw <- function(p, a, b) Crgig(n, p, a, b)
    }
  }
  type <- "exp"
  environment()
}

#' Create an object containing information about inverse chi-squared priors
#' with possibly modeled degrees of freedom and scale parameters
#'
#' @export
#' @param df degrees of freedom parameter. This can be a numeric scalar or
#'  vector of length \code{n}, the dimension of the parameter vector.
#'  Alternatively, for a scalar degrees of freedom parameter,
#'  \code{df="modeled"} or \code{df="modelled"} assign a default (gamma) prior
#'  to the degrees of freedom parameter. For more control of this gamma prior a
#'  list can be passed with some of the following components:
#'  \describe{
#'    \item{alpha0}{shape parameter of the gamma distribution}
#'    \item{beta0}{rate parameter of the gamma distribution}
# TODO MH parameters below should not be part of the prior...
#'    \item{proposal}{"RW" for random walk Metropolis-Hastings
#'      or "mala" for Metropolis-adjusted Langevin}
#'    \item{tau}{(starting) scale of Metropolis-Hastings update}
#'    \item{adapt}{whether to adapt the scale of the proposal distribution
#'      during burnin to achieve better acceptance rates.}
#'   }
#' @param scale scalar or vector scale parameter. Alternatively,
#'  \code{scale="modeled"} or \code{scale="modelled"} puts a default
#'  chi-squared prior on the scale parameter. For more control on this
#'  chi-squared prior a list can be passed with some of the following components:
#'  \describe{
#'    \item{df}{degrees of freedom (scalar or vector)}
#'    \item{scale}{scale (scalar or vector)}
#'    \item{common}{whether the modeled scale parameter of the inverse chi-squared
#'      distribution is (a scalar parameter) common to all \code{n} parameters.}
#'  }
#' @param n dimension, if known. For internal use only.
#' @param post whether conditional posterior sampling function should be
#'  created. For internal use only.
#' @return An environment with information about the prior and possibly conditional
#'  posterior distribution(s), to be used by other package functions.
# TODO sample precisions instead of variances
pr_invchisq <- function(df=1, scale=1, n=NULL, post=FALSE) {
  if (is.character(df) && df %in% c("modeled", "modelled")) df <- list()
  if (is.character(scale) && scale %in% c("modeled", "modelled")) scale <- list()
  if (!is.null(n)) {
    n <- as.integer(n)
    rprior <- function() {}
    if (post) {
      # function draw to sample from full conditional posterior
      # assumes the invchisq prior is for variance parameters of a gaussian distribution
      draw <- function(df.data, SSR) {}
    }
  }
  if (is.list(df)) {
    if (is.null(df$alpha0)) df$alpha0 <- 2  # default in prior Gamma(alpha0, beta0)
    if (is.null(df$beta0)) df$beta0 <- 0.1  # default in prior Gamma(alpha0, beta0)
    # TODO tau, proposal, adapt: move to another object
    if (is.null(df$tau)) df$tau <- 1  # (starting) scale of MH update
    if (is.null(df$proposal)) df$proposal <- "RW"
    if (is.null(df$adapt)) df$adapt <- TRUE
    if (!is.null(n)) {
      rprior_df <- function() rgamma(1L, df$alpha0, df$beta0)
      formals(rprior) <- c(alist(df=), formals(rprior))  # different signature without check NOTE
      if (post) {
        draw_df <- function(df.current, Q.current) {}  # Q.current is current precision parameter (vector)
        switch(df$proposal,
          RW = draw_df <- add(draw_df, bquote(draw_df_MH_RW(.(as.numeric(n)), df.current, Q.current, df))),
          mala = draw_df <- add(draw_df, bquote(draw_df_MH_mala(.(as.numeric(n)), df.current, Q.current, df)))
        )
        formals(draw) <- c(alist(df=), formals(draw))
      }
    }
  } else {
    if (!is.null(n)) {
      if (!(length(df) %in% c(1L, n))) stop("degrees of freedom parameter has wrong length")
      # do not enforce df > 0 to allow improper prior for data scale variance
    }
  }
  if (is.list(scale)) {
    defaults <- list(df=1, scale=1, common=FALSE)
    if (!all(names(scale) %in% names(defaults))) stop("invalid 'scale' options list")
    scale <- modifyList(defaults, scale)
    rm(defaults)
    if (!is.null(n)) {
      if (!(length(scale$df) %in% c(1L, n))) stop("degrees of freedom parameter has wrong length")
      if (!(length(scale$scale) %in% c(1L, n))) stop("scale parameter has wrong length")
      if (scale$common && n == 1L) scale$common <- FALSE
      psi0 <- scale$df / scale$scale
      if (scale$common) {
        if (length(scale$df) != 1L || length(scale$scale) != 1L) stop("scalar 'df' and 'scale' expected in common scale model")
        rprior <- add(rprior, quote(scale <- rchisq_scaled(1L, scale$df, psi=psi0)))
        rprior <- add(rprior, bquote(1 / rchisq_scaled(.(n), df, scale)))
      } else {
        rprior <- add(rprior, bquote(draw_betaprime(.(n), 0.5*scale$df, 0.5*df, df/psi0)))
      }
      if (post) {
        formals(draw) <- c(formals(draw), alist(Q=))
        # draw 1/kappa from its full conditional posterior
        if (scale$common) {
          if (is.list(df) || length(df) == 1L)
            draw <- add(draw, bquote(kappa_inv <- rchisq_scaled(1L, .(n) * df + scale$df, psi = psi0 + sum(df * Q))))
          else
            draw <- add(draw, quote(kappa_inv <- rchisq_scaled(1L, sum(df) + scale$df, psi = psi0 + sum(df * Q))))
        } else {
          draw <- add(draw, bquote(kappa_inv <- rchisq_scaled(.(n), df + scale$df, psi = psi0 + df * Q)))
        }
        draw <- add(draw, quote(psi <- df * kappa_inv))
      }
    }
  } else {
    if (!is.null(n)) {
      if (!(length(scale) %in% c(1L, n))) stop("scale parameter has wrong length")
      if (is.list(df)) {
        rprior <- add(rprior, bquote(1 / rchisq_scaled(.(n), df, scale)))
      } else {
        psi0 <- df * scale
        rprior <- add(rprior, bquote(1 / rchisq_scaled(.(n), df, psi=psi0)))
      }
      if (post) {
        draw <- add(draw, quote(psi <- df * scale))
      }
    }
  }
  if (!is.null(n) && post) {
    draw <- add(draw, bquote(1 / rchisq_scaled(.(n), df + df.data, psi=psi + SSR)))
  }
  type <- "invchisq"
  environment()
}

#' Create an object containing information about an inverse Wishart prior,
#' possibly with modeled scale matrix
#'
#' @export
#' @param df Degrees of freedom parameter. This should be a scalar numeric value.
#'  The default value is the dimension (\code{n}) plus one.
#' @param scale Either a (known) scale matrix, or
#'  \code{scale="modeled"} or \code{scale="modelled"}, which puts default
#'  chi-squared priors on the diagonal elements of the inverse Wishart scale matrix.
#'  For more control on these chi-squared priors a list can be passed with some of the
#'  following components:
#'  \describe{
#'    \item{df}{degrees of freedom (scalar or vector) of the chi-squared distribution(s)}
#'    \item{scale}{scale parameter(s) of the chi-squared distribution(s)}
#'    \item{common}{whether the modeled scale parameter of the inverse chi-squared
#'      distribution is (a scalar parameter) common to all \code{n} diagonal elements.}
#'  }
#' @param n dimension, if known. For internal use only.
#' @return An environment with information about the prior distribution used,
#'  to be used by other package functions.
#' @references
#'  A. Huang and M.P. Wand (2013).
#'    Simple marginally noninformative prior
#'    distributions for covariance matrices.
#'    Bayesian Analysis 8, 439-452.
pr_invwishart <- function(df=NULL, scale=NULL, n=NULL) {
  if (is.null(n)) {
    if (is_a_matrix(scale)) n <- dim(scale)[1L]
  }
  if (!is.null(n) && n == 1L) stop("pr_invwishart does not support scalar variance")
  if (is.null(df)) {
    if (!is.null(n)) df <- n + 1  # default number of degrees of freedom
  } else {
    if (length(df) != 1L) stop("degrees of freedom parameter for inverse Wishart must be scalar")
    # to allow improper priors, do not enforce df > n - 1
  }
  if (is.null(scale)) {
    if (!is.null(n)) scale <- diag(n)
  } else {
    if (is.character(scale) && scale %in% c("modeled", "modelled")) scale <- list()
    if (is.list(scale)) {  # Huang-Wand prior
      if (is.null(scale$df)) scale$df <- 1
      if (is.null(scale$scale)) scale$scale <- 1
      if (is.null(scale$common)) scale$common <- FALSE  # by default (HW prior)
      if (!is.null(n)) {
        if (!(length(scale$df) %in% c(1L, n))) stop("degrees of freedom parameter has wrong length")
        if (!(length(scale$scale) %in% c(1L, n))) stop("scale parameter has wrong length")
        if (scale$common) {
          if (length(scale$df) != 1L || length(scale$scale) != 1L) stop("scalar 'df' and 'scale' expected in common scale model")
        }
      }
    } else {
      if (!is.null(n) && !identical(dim(scale), c(n, n))) stop("incompatible scale matrix")
    }
  }
  list(type="invwishart", df=df, scale=scale)
}
