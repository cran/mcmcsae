#' Create a model component object for a regression (fixed effects) component
#' in the linear predictor with measurement errors in quantitative covariates
#'
#' This function is intended to be used on the right hand side of the
#' \code{formula} argument to \code{\link{create_sampler}} or
#' \code{\link{generate_data}}. It creates an additive regression term in the
#' model's linear predictor. By default, the prior for the regression
#' coefficients is improper uniform. If \code{b0} or \code{Q0} are specified
#' the prior becomes normal with mean \code{b0} (default 0) and variance
#' (matrix) \code{sigma_^2 Q0^-1} where \code{sigma_^2} is the overall scale
#' parameter of the model, if any. Covariates are assumed to be measured subject
#' to normally distributed errors with zero mean and variance specified using
#' the \code{formula} or \code{V} arguments. Note that this means that \code{formula}
#' should only contain quantitative variables, and no intercept.
#'
#' @examples
#' \donttest{
#' # example of Ybarra and Lohr (2008)
#' m <- 50
#' X <- rnorm(m, mean=5, sd=3)  # true covariate values
#' v <- rnorm(m, sd=2)
#' theta <- 1 + 3*X + v  # true values
#' psi <- rgamma(m, shape=4.5, scale=2)
#' e <- rnorm(m, sd=sqrt(psi))  # sampling error
#' y <- theta + e  # direct estimates
#' C <- c(rep(3, 10), rep(0, 40))  # measurement error for first 10 values
#' W <- X + rnorm(m, sd=sqrt(C))  # covariate subject to measurement error
#'
#' # fit Ybarra-Lohr model
#' sampler <- create_sampler(
#'   y ~ 1 + mec(~ 0 + W, V=C) + gen(factor=~local_),
#'   Q0=1/psi, sigma.fixed=TRUE, linpred="fitted"
#' )
#' sim <- MCMCsim(sampler, n.iter=800, n.chain=2, store.all=TRUE, verbose=FALSE)
#' (summ <- summary(sim))
#' plot(X, W, xlab="true X", ylab="inferred X")
#' points(X, summ$mec2_X[, "Mean"], col="green")
#' abline(0, 1, col="red")
#' legend("topleft", legend=c("prior mean", "posterior mean"), col=c("black", "green"), pch=c(1,1))
#' }
#'
#' @export
#' @param formula a formula specifying the predictors subject to measurement error
#'  and possibly their variances as well. In the latter case the formula syntax
#'  \code{~ (x1 | V.x1) + (x2 | V.x2) + ...} should be used where \code{x1, x2, ...}
#'  are the names of (quantitative) predictors and \code{V.x1, V.x2, ...} are the names
#'  of the variables holding the corresponding measurement error variances.
#'  If only the predictors are specified
#'  the formula has the usual form \code{~ x1 + x2 + ...}. In that case variances
#'  should be specified using argument \code{V}.
#'  All variable names are looked up in the data frame
#'  passed as \code{data} argument to \code{\link{create_sampler}} or
#'  \code{\link{generate_data}}, or in \code{environment(formula)}.
#' @param sparse whether the model matrix associated with \code{formula} should
#'  be sparse. The default is to base this on a simple heuristic.
#' @param X a (possibly sparse) design matrix can be specified directly, as an
#'  alternative to the creation of one based on \code{formula}. If \code{X} is
#'  specified \code{formula} is ignored.
#' @param V measurement error variance; can contain zeros
# TODO in general case this can be a list of q x q precision matrices where q
#      is the number of covariates
#' @param Q0 prior precision matrix for the regression effects. The default is a
#'  zero matrix corresponding to a noninformative improper prior.
#'  It can be specified as a scalar value, as a numeric vector of appropriate
#'  length, or as a matrix object.
#' @param b0 prior mean for the regression effect. Defaults to a zero vector.
#'  It can be specified as a scalar value or as a numeric vector of
#'  appropriate length.
#' @param R optional constraint matrix for equality restrictions R'x = r where
#'  \code{x} is the vector of regression effects.
#' @param r right hand side for the equality constraints.
#' @param S optional constraint matrix for inequality constraints S'x >= s where
#'  x is the vector of regression effects.
#' @param s right hand side for the inequality constraints.
#' @param lower as an alternative to \code{s}, \code{lower} and \code{upper} may be specified
#'  for two-sided constraints lower <= S'x <= upper.
#' @param upper as an alternative to \code{s}, \code{lower} and \code{upper} may be specified
#'  for two-sided constraints lower <= S'x <= upper.
#' @param name the name of the model component. This name is used in the output of the
#'  MCMC simulation function \code{\link{MCMCsim}}. By default the name will be 'reg'
#'  with the number of the model term attached.
#' @param perm whether permutation should be used in the Cholesky decomposition used for
#'  updating the model component's coefficient. Default is based on a simple heuristic.
#' @param debug if \code{TRUE} a breakpoint is set at the beginning of the posterior
#'  draw function associated with this model component. Mainly intended for developers.
#' @param e for internal use only.
#' @return an object with precomputed quantities and functions for sampling from
#'  prior or conditional posterior distributions for this model component. Only intended
#'  for internal use by other package functions.
#' @references
#'  L.M. Ybarra and S.L. Lohr (2008).
#'    Small area estimation when auxiliary information is measured with error.
#'    Biometrika 95(4), 919-931.
#'
#'  S. Arima, G.S. Datta and B. Liseo (2015).
#'    Bayesian estimators for small area models when auxiliary information is measured with error.
#'    Scandinavian Journal of Statistics 42(2), 518-529.
mec <- function(formula = ~ 1, sparse=NULL, X=NULL, V=NULL, Q0=NULL, b0=NULL,
                R=NULL, r=NULL, S=NULL, s=NULL, lower=NULL, upper=NULL,
                name="", perm=NULL,
                debug=FALSE, e=parent.frame()) {

  type <- "mec"
  if (name == "") stop("missing model component name")
  remove.redundant <- FALSE

  if (is.null(X)) {
    # TODO: warn if formula contains explicit intercept or categorical terms
    vs <- as.list(attr(terms(formula), "variables")[-1L])
    if (all(sapply(vs, length) == 3L)) {
      if (!all(sapply(vs, function(x) identical(x[[1L]], as.name("|"))))) stop("invalid 'formula' argument")
      V.in.formula <- TRUE
      formula.X <- as.formula(paste0("~ 0 + ", paste(sapply(vs, function(x) deparse(x[[2L]])), collapse=" + ")))
      formula.V <- as.formula(paste0("~ 0 + ", paste(sapply(vs, function(x) deparse(x[[3L]])), collapse=" + ")))
      X <- model_Matrix(formula.X, data=e$data, remove.redundant=FALSE, sparse=sparse)
      V <- model_Matrix(formula.V, data=e$data, remove.redundant=FALSE, sparse=sparse)
    } else {
      if (all(sapply(vs, length) <= 2L)) {
        V.in.formula <- FALSE
        formula.X <- as.formula(paste0("~ 0 + ", paste(sapply(vs, deparse), collapse=" + ")))
        X <- model_Matrix(formula.X, data=e$data, remove.redundant=FALSE, sparse=sparse)
      } else {
        stop("invalid 'formula' argument")
      }
    }
  } else {
    V.in.formula <- FALSE
  }
  if (nrow(X) != e$n) stop("design matrix with incompatible number of rows")
  e$coef.names[[name]] <- colnames(X)
  X <- unname(X)
  q <- ncol(X)
  in_block <- name %in% unlist(e$block)

  if (!V.in.formula) {
    if (is.null(V)) stop("no measurement error variances supplied in mec component")
    if (is.vector(V)) {
      if (q != 1L) stop("'V' should be a matrix")
      V <- matrix(V, ncol=1L)
    }
    if (nrow(V) != e$n || ncol(V) != q) stop("wrong size of variance matrix 'V'")
    # TODO allow n x q matrix (uncorrelated measurement errors), or list of n q x q symmetric matrices
  }
  # determine units with measurement error
  V <- unname(V)
  if (any(V < 0)) stop("negative measurement error variance(s)")
  i.me <- which(apply(V, 1L, function(x) all(x > 0)))  # TODO add tolerance + allow list of matrix
  nme <- length(i.me)
  if (nme == e$n) i.me <- 1:e$n  # use ALTREP (?)
  if (q == 1L) {
    QME <- 1 / V[i.me, ]
    rm(V)
  } else {
    V <- V[i.me, ]
  }

  if (e$Q0.type == "unit") {
    if (q == 1L) v0 <- 1 else v0 <- rep.int(1, nme)
  } else {
    if (e$Q0.type == "symm") stop("not supported: measurement error component and non-diagonal sampling variance")
    v0 <- 1/diag(e$Q0)[i.me]
  }

  Q0 <- set_prior_precision(Q0, q, sparse=if (in_block) TRUE else NULL)  # mc_block uses x-slot of precision matrix
  informative.prior <- !is_zero_matrix(Q0)
  if (!informative.prior || is.null(b0) || all(b0 == 0)) {
    b0 <- 0
    Q0b0 <- 0
    zero.mean <- TRUE
  } else {
    if (in_block) stop("not yet implemented: block sampling with non-zero prior mean")
    if (length(b0) == 1L) {
      Q0b0 <- Q0 %m*v% rep.int(as.numeric(b0), q)
    } else {
      if (length(b0) == q)
        Q0b0 <- Q0 %m*v% as.numeric(b0)
      else
        stop("wrong input for prior mean 'b0' in model component", name)
    }
    zero.mean <- FALSE
  }

  if (!zero.mean && e$compute.weights) stop("weights cannot be computed if some coefficients have non-zero prior means")
  log_prior_0 <- -0.5 * q * log(2*pi)
  if (informative.prior) {
    d <- diag(suppressWarnings(chol(Q0, pivot=TRUE)))
    log_prior_0 <- log_prior_0 + sum(log(d[d > 0]))
    rm(d)
  }

  name_X <- paste0(name, "_X")
  linpred <- function(p) p[[name_X]] %m*v% p[[name]]

  # function that creates an oos prediction function (closure) based on new data
  # this computes the contribution of this component to the linear predictor
  make_predict <- function(newdata) {
    if (!V.in.formula) stop("for out-of-sample prediction with mec components make sure to specify the",
                            "measurement error variances using the formula argument")
    Xnew <- model_Matrix(formula.X, data=newdata, remove.redundant=FALSE, sparse=sparse)
    dimnames(Xnew) <- list(NULL, NULL)
    Vnew <- model_Matrix(formula.V, data=newdata, remove.redundant=FALSE, sparse=sparse)
    dimnames(Vnew) <- list(NULL, NULL)
    if (any(Vnew < 0)) stop("negative measurement error variance(s) in 'newdata'")
    i.me.new <- which(apply(Vnew, 1L, function(x) all(x > 0)))  # TODO add tolerance + allow list of matrix
    if (length(i.me.new) == nrow(newdata)) i.me.new <- seq_len(nrow(newdata))
    rm(newdata)
    if (nrow(Vnew) < 0.5 * nrow(Xnew)) {
      rm(i.me.new)
      pred <- function(p) (Xnew + sqrt(Vnew) * Crnorm(length(Vnew))) %m*v% p[[name]]
    } else {
      Vnew <- Vnew[i.me.new, ]
      pred <- function(p) {
        out <- Xnew
        out[i.me.new, ] <- out[i.me.new, ] + sqrt(Vnew) * Crnorm(length(Vnew))
        out %m*v% p[[name]]
      }
    }
    pred
  }

  rprior <- function(p) {}
  # X, QME are assumed to be part of the prior information
  if (nme == e$n)
    rprior <- add(rprior, bquote(p[[.(name_X)]] <- X + sqrt(.(if (q == 1L) quote(1/QME) else quote(V))) * Crnorm(length(QME))))
  else {
    rprior <- add(rprior, bquote(p[[.(name_X)]] <- X))
    rprior <- add(rprior, bquote(p[[.(name_X)]][i.me, ] <- sqrt(.(if (q == 1L) quote(1/QME) else quote(V))) * Crnorm(length(QME))))
  }
  rprior <- add(rprior, bquote(p[[.(name)]] <- b0 + drawMVN_Q(Q0, sd=.(if (e$sigma.fixed) 1 else quote(p[["sigma_"]])))))
  rprior <- add(rprior, quote(p))
  # TODO TMVN prior


  if (!e$prior.only) {

    if (q > 1L) base_tcrossprod <- base::tcrossprod

    draw <- function(p) {}
    if (debug) draw <- add(draw, quote(browser()))
    if (e$single.block && length(e$mod) == 1L) {
      draw <- add(draw, quote(p$e_ <- e$y_eff()))
    } else {
      if (e$e.is.res)
        draw <- add(draw, bquote(p$e_ <- p[["e_"]] + p[[.(name_X)]] %m*v% p[[.(name)]]))
      else
        draw <- add(draw, bquote(p$e_ <- p[["e_"]] - p[[.(name_X)]] %m*v% p[[.(name)]]))
    }

    # 1. draw covariates p[[name_X]]
    # for now we assume that X is dense
    draw <- add(draw, bquote(scale <- .(if (e$sigma.fixed) quote(v0) else quote(v0 * p[["sigma_"]]^2))))
    if (q == 1L) {
      draw <- add(draw, bquote(Vscaled <- 1 / (QME * scale + p[[.(name)]]^2)))
      if (nme == e$n) {
        draw <- add(draw, bquote(p[[.(name_X)]] <- X + Vscaled * p[[.(name)]] * (p[["e_"]] - p[[.(name)]] * X) + sqrt(scale * Vscaled) * Crnorm(nme)))
      } else {
        draw <- add(draw, bquote(p[[.(name_X)]][i.me] <- X[i.me] + Vscaled * p[[.(name)]] * (p[["e_"]][i.me] - p[[.(name)]] * X[i.me]) + + sqrt(scale * Vscaled) * Crnorm(nme)))
      }
    } else {
      # V-form (independent m.e.)
      draw <- add(draw, bquote(
        for (i in i.me) {
          Vi <- V[i, ]
          Cb <- Vi * p[[.(name)]]
          f <- 1 / (scale[i] + dotprodC(p[[.(name)]], Cb))
          Vx <- add_diagC(- f * base_tcrossprod(Cb), Vi)
          Xi <- X[i, ]
          p[[.(name_X)]][i, ] <- Xi + Cb * f * (p[["e_"]][i] - dotprodC(p[[.(name)]], Xi)) + crossprod_mv(chol.default(Vx), rnorm(q))
        }
      ))
      # TODO use Rcpp, and use Q-form at least for cases where all QME[i, ] > 0, and handle list-of-matrix (correlated m.e.)
    }

    if (!in_block) {
      # 2. draw coefficients p[[name]]
      if (all(Q0b0 == 0))
        draw <- add(draw, bquote(Xy <- crossprod_mv(p[[.(name_X)]], e$Q_e(p))))
      else
        draw <- add(draw, bquote(Xy <- crossprod_mv(p[[.(name_X)]], e$Q_e(p)) + Q0b0))

      # TODO modeled.Q case: use p[["Q_"]]
      #if (e$modeled.Q) {  # in this case both XX and Xy must be updated in each iteration
      if (e$Q0.type == "symm") {
        # TODO this case not (yet) supported!
        draw <- add(draw, bquote(XX <- crossprod_sym(p[[.(name_X)]], p[["QM_"]])))
      } else {
        #draw <- add(draw, quote(XX <- crossprod_sym(p[[name_X]], p[["Q_"]])))
        if (informative.prior) {
          # TODO matsum template? or can we always assume dense matrices?
          draw <- add(draw, bquote(XX <- crossprod_sym(p[[.(name_X)]], e$Q0) + Q0))
        } else
          draw <- add(draw, bquote(XX <- crossprod_sym(p[[.(name_X)]], e$Q0)))
      }
      if (informative.prior) {
        # TODO instead of runif(n, ...) account for the structure in Qmod; lambda may not vary per unit!
        # TODO here we need M1 to be dense without zeros
        mat_sum <- make_mat_sum(M0=Q0, M1=crossprod_sym(X, crossprod_sym(Diagonal(x=runif(e$n, 0.9, 1.1)), e$Q0)))
        MVNsampler <- create_TMVN_sampler(
          Q=mat_sum(crossprod_sym(X, crossprod_sym(Diagonal(x=runif(e$n, 0.9, 1.1)), e$Q0))),
          update.Q=TRUE, name=name, R=R, r=r, S=S, s=s, lower=lower, upper=upper
        )
        draw <- add(draw, bquote(p[[.(name)]] <- MVNsampler$draw(p, .(if (e$sigma.fixed) 1 else quote(p[["sigma_"]])), Q=mat_sum(XX), Xy=Xy)[[.(name)]]))
      } else {
        MVNsampler <- create_TMVN_sampler(
          Q=crossprod_sym(X, crossprod_sym(Diagonal(x=runif(e$n, 0.9, 1.1)), e$Q0)),
          update.Q=TRUE, name=name, R=R, r=r, S=S, s=s, lower=lower, upper=upper
        )
        draw <- add(draw, bquote(p[[.(name)]] <- MVNsampler$draw(p, .(if (e$sigma.fixed) 1 else quote(p[["sigma_"]])), Q=XX, Xy=Xy)[[.(name)]]))
      }
    }  # END if (!in_block)

    # compute full residuals/fitted values
    if (e$e.is.res)
      draw <- add(draw, bquote(p$e_ <- p[["e_"]] - p[[.(name_X)]] %m*v% p[[.(name)]]))
    else
      draw <- add(draw, bquote(p$e_ <- p[["e_"]] + p[[.(name_X)]] %m*v% p[[.(name)]]))
    draw <- add(draw, quote(p))

    start <- function(p) {
      if (!in_block) p <- MVNsampler$start(p)
      p[[name_X]] <- X
      p
    }
  }  # END if (!e$prior.only)

  rm(R, r, S, s, lower, upper)

  #rm(e)
  environment()
}
