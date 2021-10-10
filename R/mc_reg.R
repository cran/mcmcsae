#' Create a model component object for a regression (fixed effects) component
#' in the linear predictor
#'
#' This function is intended to be used on the right hand side of the
#' \code{formula} argument to \code{\link{create_sampler}} or
#' \code{\link{generate_data}}. It creates an additive regression term in the
#' model's linear predictor. By default, the prior for the regression
#' coefficients is improper uniform. If \code{b0} or \code{Q0} are specified
#' the prior becomes normal with mean \code{b0} (default 0) and variance
#' (matrix) \code{sigma_^2 Q0^-1} where \code{sigma_^2} is the overall scale
#' parameter of the model, if any.
#'
#' @examples
#' \donttest{
#' data(iris)
#' # default: flat priors on regression coefficients
#' sampler <- create_sampler(Sepal.Length ~
#'     reg(~ Petal.Length + Species, name="beta"),
#'   data=iris
#' )
#' sim <- MCMCsim(sampler, burnin=100, n.iter=400)
#' summary(sim)
#' # (weakly) informative normal priors on regression coefficients
#' sampler <- create_sampler(Sepal.Length ~
#'     reg(~ Petal.Length + Species, Q0=1e-2, name="beta"),
#'   data=iris
#' )
#' sim <- MCMCsim(sampler, burnin=100, n.iter=400)
#' summary(sim)
#' # binary regression
#' sampler <- create_sampler(Species == "setosa" ~
#'     reg(~ Sepal.Length, Q=0.1, name="beta"),
#'   family="binomial", data=iris)
#' sim <- MCMCsim(sampler, burnin=100, n.iter=400)
#' summary(sim)
#' pred <- predict(sim)
#' str(pred)
#' # example with equality constrained regression effects
#' n <- 500
#' df <- data.frame(x=runif(n))
#' df$y <- rnorm(n, 1 + 2*df$x)
#' R <- matrix(1, 2, 1)
#' r <- 3
#' sampler <- create_sampler(y ~ reg(~ 1 + x, R=R, r=r, name="beta"), data=df)
#' sim <- MCMCsim(sampler)
#' summary(sim)
#' plot(sim, "beta")
#' summary(transform_dc(sim$beta, fun=function(x) crossprod_mv(R, x) - r))
#' }
#'
#' @export
#' @param formula a formula specifying the predictors to be used in the model,
#'  in the same way as the right hand side of the \code{formula} argument of
#'  R's \code{lm} function. Variable names are looked up in the data frame
#'  passed as \code{data} argument to \code{\link{create_sampler}} or
#'  \code{\link{generate_data}}, or in \code{environment(formula)}.
#' @param remove.redundant whether redundant columns should be removed from the
#'  design matrix. Default is \code{FALSE}. But note that treatment contrasts are
#'  automatically applied to all factor variables in \code{formula}.
# TODO allow to pass contrasts.arg to model_matrix (e.g. "contr.none", and later "contr.sparse")
#' @param sparse whether the model matrix associated with \code{formula} should
#'  be sparse. The default is to base this on a simple heuristic.
#' @param X a (possibly sparse) design matrix can be specified directly, as an
#'  alternative to the creation of one based on \code{formula}. If \code{X} is
#'  specified \code{formula} is ignored.
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
reg <- function(formula = ~ 1, remove.redundant=FALSE, sparse=NULL, X=NULL,
                Q0=NULL, b0=NULL, R=NULL, r=NULL, S=NULL, s=NULL, lower=NULL, upper=NULL,
                name="", perm=NULL,
                debug=FALSE, e=parent.frame()) {

  type <- "reg"
  if (name == "") stop("missing model component name")

  if (is.null(X)) {
    if (e$family$family == "multinomial") {
      edat <- new.env(parent = .GlobalEnv)
      X <- NULL
      for (k in seq_len(e$Km1)) {
        edat$cat_ <- factor(rep.int(e$cats[k], e$n0), levels=e$cats[-length(e$cats)])
        X <- rbind(X, model_matrix(formula, data=e$data, sparse=sparse, enclos=edat))
      }
      rm(edat)
    } else {
      X <- model_matrix(formula, e$data, sparse=sparse)
    }
    if (remove.redundant) X <- remove_redundancy(X)
  }
  X <- economizeMatrix(X, sparse=sparse, strip.names=FALSE)

  if (nrow(X) != e$n) stop("design matrix with incompatible number of rows")
  e$coef.names[[name]] <- colnames(X)
  X <- unname(X)
  q <- ncol(X)
  in_block <- name %in% unlist(e$block)

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

  linpred <- function(p) X %m*v% p[[name]]

  make_predict <- function(newdata) {
    nnew <- nrow(newdata)
    if (e$family$family == "multinomial") {
      edat <- new.env(parent = .GlobalEnv)
      Xnew <- NULL
      for (k in seq_len(e$Km1)) {
        edat$cat_ <- factor(rep.int(e$cats[k], nnew), levels=e$cats[-length(e$cats)])
        Xnew <- rbind(Xnew, model_matrix(formula, data=newdata, sparse=sparse, enclos=edat))
      }
      rm(edat)
    } else {
      Xnew <- model_matrix(formula, newdata, sparse=sparse)
    }
    # for prediction do not (automatically) remove redundant columns!
    if (remove.redundant)
      Xnew <- Xnew[, e$coef.names[[name]], drop=FALSE]
    # check that X has the right number of columns
    if (ncol(Xnew) != q) stop("'newdata' yields ", ncol(Xnew), " predictor column(s) for model term '", name, "' versus ", q, " originally")
    # TODO check the category names as well; allow oos categories; insert missing categories in newdata
    economizeMatrix(Xnew, sparse=sparse, strip.names=TRUE)
  }

  rprior <- function(p) {}
  rprior <- add(rprior, bquote(b0 + drawMVN_Q(Q0, sd=.(if (e$sigma.fixed) 1 else quote(p[["sigma_"]])))))
  # TODO TMVN prior

  if (!in_block && !e$prior.only) {

    draw <- function(p) {}
    if (debug) draw <- add(draw, quote(browser()))
    if (e$single.block) {
      draw <- add(draw, quote(p$e_ <- e$y_eff()))
    } else {
      if (e$e.is.res)
        draw <- add(draw, bquote(p$e_ <- p[["e_"]] + X %m*v% p[[.(name)]]))
      else
        draw <- add(draw, bquote(p$e_ <- p[["e_"]] - X %m*v% p[[.(name)]]))
    }
    if (e$single.block && !e$modeled.Q && e$family$link != "probit") {
      # single regression component, no variance modeling, Xy fixed
      Xy <- crossprod_mv(X, e$Q0 %m*v% e$y_eff()) + Q0b0
    } else {  # Xy updated in each iteration
      if (all(Q0b0 == 0))
        draw <- add(draw, quote(Xy <- crossprod_mv(X, e$Q_e(p))))
      else
        draw <- add(draw, quote(Xy <- crossprod_mv(X, e$Q_e(p)) + Q0b0))
    }

    if (e$modeled.Q) {  # in this case both XX and Xy must be updated in each iteration
      if (e$Q0.type == "symm")
        draw <- add(draw, quote(XX <- crossprod_sym(X, p[["QM_"]])))
      else
        draw <- add(draw, quote(XX <- crossprod_sym(X, p[["Q_"]])))
      if (informative.prior) {
        # TODO instead of runif(n, ...) account for the structure in Qmod; lambda may not vary per unit!
        mat_sum <- make_mat_sum(M0=Q0, M1=crossprod_sym(X, crossprod_sym(Diagonal(x=runif(e$n, 0.9, 1.1)), e$Q0)))
        MVNsampler <- create_TMVN_sampler(
          Q=mat_sum(crossprod_sym(X, crossprod_sym(Diagonal(x=runif(e$n, 0.9, 1.1)), e$Q0))),
          update.Q=TRUE, name=name, R=R, r=r, S=S, s=s, lower=lower, upper=upper
        )
        draw <- add(draw, bquote(p <- MVNsampler$draw(p, .(if (e$sigma.fixed) 1 else quote(p[["sigma_"]])), Q=mat_sum(XX), Xy=Xy)))
      } else {
        MVNsampler <- create_TMVN_sampler(
          Q=crossprod_sym(X, crossprod_sym(Diagonal(x=runif(e$n, 0.9, 1.1)), e$Q0)),
          update.Q=TRUE, name=name, R=R, r=r, S=S, s=s, lower=lower, upper=upper
        )
        draw <- add(draw, bquote(p <- MVNsampler$draw(p, .(if (e$sigma.fixed) 1 else quote(p[["sigma_"]])), Q=XX, Xy=Xy)))
      }
    } else {  # precision matrix XX + Q0 not updated
      if (e$single.block && e$family$link != "probit") {
        MVNsampler <- create_TMVN_sampler(
          Q=crossprod_sym(X, e$Q0) + Q0, Xy=Xy,
          R=R, r=r, S=S, s=s, lower=lower, upper=upper, perm=perm,
          name=name
        )
        draw <- add(draw, bquote(p[[.(name)]] <- MVNsampler$draw(p, .(if (e$sigma.fixed) 1 else quote(p[["sigma_"]])))[[.(name)]]))
        rm(Xy)
      } else {
        MVNsampler <- create_TMVN_sampler(
          Q=crossprod_sym(X, e$Q0) + Q0,
          update.mu=TRUE,
          R=R, r=r, S=S, s=s, lower=lower, upper=upper, perm=perm,
          name=name
        )
        draw <- add(draw, bquote(p[[.(name)]] <- MVNsampler$draw(p, .(if (e$sigma.fixed) 1 else quote(p[["sigma_"]])), Xy=Xy)[[.(name)]]))
      }
    }
    if (e$e.is.res)
      draw <- add(draw, bquote(p$e_ <- p[["e_"]] - X %m*v% p[[.(name)]]))
    else
      draw <- add(draw, bquote(p$e_ <- p[["e_"]] + X %m*v% p[[.(name)]]))
    draw <- add(draw, quote(p))

    start <- function(p) MVNsampler$start(p)
  }  # END if (!in_block)

  rm(R, r, S, s, lower, upper)

  #rm(e)
  environment()
}

# experimental: display the reg component's prior
show_reg <- function(mc) {
  mod_str <- paste(mc$name, "~")
  mod_str <- paste(mod_str,
    if (mc$informative.prior) {
      if (mc$e$sigma.fixed) {
        "N(b0, Q0^-1)"
      } else {
        "N(b0, sigma_^2 Q0^-1)"
      }
    } else {
      "1"
    }
  )
  if (!is.null(mc$MVNsampler)) {
    # TODO MVNsampler currently only used for posterior sampling
    if (mc$MVNsampler$eq) mod_str <- paste0(mod_str, ",  R'", mc$name, " = r")
    if (mc$MVNsampler$ineq) mod_str <- paste0(mod_str, ",  S'", mc$name, " >= s")
  }
  cat(mod_str, "\n")
}
