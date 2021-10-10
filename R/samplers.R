#' Create a sampler object
#'
#' This function sets up a sampler object, based on the specification of a model. The object contains functions to
#' draw a set of model parameters from their prior and conditional posterior distributions, and
#' to generate starting values for the MCMC simulation. The functions share a common environment
#' containing precomputed quantities such as design matrices based on the model and the data.
#' The sampler object is the main input for the MCMC simulation function \code{\link{MCMCsim}}.
#'
#' The right hand side of the \code{formula} argument to \code{create_sampler} can be used to specify
#' additive model components. Currently two specialized model components are supported, \code{\link{reg}(...)}
#' and \code{\link{gen}(...)} for regression and generic random effects components, respectively.
#'
#' For gaussian models, \code{formula.V} can be used to specify the variance structure of the model.
#' Currently two specialized variance model components are supported, \code{\link{vreg}(...)} for regression
#' effects predicting the log-variance and \code{\link{vfac}(...)} for modeled variance factors.
#'
#' Further computational options can be set using the \code{control} parameter, which should be
#' passed as a list with possible elements
#' \describe{
#'   \item{add.outer.R}{whether to add the outer product of the constraint matrix for a better conditioned solve system
#'     for blocks. This is done by default when using blocked Gibbs sampling for blocks with constraints.}
#'   \item{recompute.e}{when \code{FALSE} residuals or linear predictors are only computed at the start of the simulation.
#'     This may give a modest speedup but in some cases may be less accurate due to roundoff error accumulation.
#'     Default is \code{TRUE}.}
#'   \item{CG}{for a conjugate gradient iterative algorithm instead of Cholesky updates for sampling
#'     the model's coefficients. This must be a list with possible components \code{max.it},
#'     \code{stop.criterion}, \code{verbose}, \code{preconditioner} and \code{scale},
#'     see \code{\link{setup_CG_sampler}}}. This is currently an experimental feature.
#' }
#'
#' @examples
#' # first generate some data
#' n <- 200
#' x <- rnorm(n)
#' y <- 0.5 + 2*x + 0.3*rnorm(n)
#' # create a sampler for a simple linear regression model
#' sampler <- create_sampler(y ~ x)
#' sim <- MCMCsim(sampler)
#' (summary(sim))
#'
#' y <- rbinom(n, 1, 1 / (1 + exp(-(0.5 + 2*x))))
#' # create a sampler for a binary logistic regression model
#' sampler <- create_sampler(y ~ x, family="binomial")
#' sim <- MCMCsim(sampler)
#' (summary(sim))
#'
#' @export
#' @param formula formula to specify the response variable and additive model components. The model components
#'  form the linear predictor part of the model. A model component on the right hand side can be either
#'  a regression term specified by \code{\link{reg}(...)}, a covariates subject to error term specified
#'  by \code{\link{mec}(...)}, or a generic random effect term specified by \code{\link{gen}(...)}.
#'  See for details the help pages for these model component creation functions.
#'  An offset can be specified as \code{offset(...)}.
#'  Other terms in the formula are collectively interpreted as ordinary regression effects,
#'  treated in the same way as a \code{reg(...)} term, but without the option to change the prior.
#' @param data data frame with n rows in which the variables specified in model components (if any) can be found.
#' @param family character string describing the data distribution. The default is 'gaussian'.
#'  Other options are 'binomial' for the binomial distribution and 'negbinomial' for the negative binomial distribution,
#'  and "poisson" for the Poisson distribution.
#'  For the binomial distribution logistic and probit link functions are supported, the latter only
#'  for binary data. For the negative binomial and Poisson distributions a log link function is assumed.
#'  Note that posterior inference based on the Poisson distribution is implemented only approximately,
#'  as a special case of the negative binomial distribution.
#'  For categorical or multinomial data, \code{family = "multinomial"} can be used. The implementation
#'  is based on a stick-breaking representation of the multinomial distribution, and the logistic link
#'  function relates each category except the last to a linear predictor. The categories can be
#'  referenced in the model specification formula by 'cat_'.
#' @param ny in case \code{family="binomial"} the (vector of) numbers of trials.
#'  It can be either a numeric vector or the name of a variable in \code{data}.
#'  Defaults to a vector of 1s.
#' @param ry in case \code{family="negbinomial"} the known, i.e. fixed part of the (reciprocal)
#'  dispersion parameter. It can be specified either as a numeric vector or the name of a
#'  numeric variable in \code{data}. The overall dispersion parameter is the product of \code{ry}
#'  with a positive scalar factor modelled as specified by argument \code{r.mod}. By default
#'  \code{ry} is taken to be 1. For \code{family = "poisson"} a single value can be specified,
#'  determining how well the Poisson distribution is approximated by the negative binomial distribution.
#'  The value should be large enough such that the negative binomial's overdispersion
#'  becomes negligible, but not too large as this might result in slow MCMC mixing. The default is
#'  \code{ry=100} in this case.
#' @param r.mod prior specification for a scalar (reciprocal) dispersion parameter
#'  of the negative binomial distribution. The prior can be specified by a call to a prior
#'  specification function. Currently \code{\link{pr_invchisq}}, \code{\link{pr_gig}} and
#'  \code{\link{pr_fixed}} are supported. The default is a chi-squared prior with 1 degree
#'  of freedom. To set the overall dispersion parameter to the value(s) specified by \code{ry},
#'  use \code{r.mod = pr_fixed(value=1)}.
#' @param sigma.fixed for Gaussian models, if \code{TRUE} the residual standard deviation parameter 'sigma_' is fixed at 1. In that case
#'  argument \code{sigma.mod} is ignored. This is convenient for Fay-Herriot type models with (sampling) variances assumed to be known.
#'  Default is \code{FALSE}.
#' @param sigma.mod prior for the variance parameter of a gaussian sampling distribution.
#'  This can be specified by a call to one of the prior specification functions
#'  \code{\link{pr_invchisq}}, \code{\link{pr_exp}}, \code{\link{pr_gig}} or \code{\link{pr_fixed}} for
#'  inverse chi-squared, exponential, generalized inverse gaussian or degenerate prior distribution,
#'  respectively. The default is an improper prior \code{pr_invchisq(df=0, scale=1)}. A half-t prior on the
#'  standard deviation can be specified using \code{\link{pr_invchisq}} with a chi-squared distributed scale
#'  parameter.
#' @param Q0 n x n data-level precision matrix for a Gaussian model. It defaults to the unit matrix.
#'  If an n-vector is provided it will be expanded to a (sparse) diagonal matrix with Q0 on its diagonal.
#'  If a name is supplied it will be looked up in \code{data} and subsequently expanded to a diagonal matrix.
#' @param formula.V a formula specifying the terms of a variance model in the case of a Gaussian likelihood.
#'  Currently two types of terms are supported: a regression term for the log-variance
#'  specified with \code{\link{vreg}(...)}, and a term \code{\link{vfac}(...)} for multiplicative modeled factors
#'  at a certain level specified by a factor variable. By using unit-level inverse-chi-squared factors the marginal
#'  sampling distribution becomes a Student-t distribution, and by using unit-level exponential factors it becomes
#'  a Laplace or double exponential distribution.
#' @param logJacobian if the data are transformed the logarithm of the Jacobian can be supplied so that it
#'  is incorporated in all log-likelihood computations. This can be useful for comparing information criteria
#'  for different transformations. It should be supplied as a vector of the same size as the response variable,
#'  and is currently only supported if \code{family="gaussian"}.
#'  For example, when a log-transformation is used on response vector \code{y}, the vector \code{-log(y)}
#'  should be supplied.
#' @param linpred a list of matrices defining (possibly out-of-sample) linear predictors to be simulated.
#'  This allows inference on e.g. (sub)population totals or means. The list must be of the form
#'  \code{list(name_1=X_1, ...)} where the names refer to the model component names and predictions are
#'  computed by summing \code{X_i \%*\% p[[name_i]]}. Alternatively, \code{linpred="fitted"} can be used
#'  for simulations of the full in-sample linear predictor.
#' @param compute.weights if \code{TRUE} weights are computed for each element of \code{linpred}. Note that for
#'  a large dataset in combination with vector-valued linear predictors the weights can take up a lot of memory.
#'  By default only means are stored in the simulation carried out using \code{\link{MCMCsim}}.
#' @param block if \code{TRUE} all coefficients are sampled in a single block. Alternatively, a list of
#'  character vectors indicating which coefficients should be sampled together in blocks.
#' @param prior.only whether a sampler is set up only for sampling from the prior or for sampling from both prior
#'  and posterior distributions. Default \code{FALSE}. If \code{TRUE} there is no need to specify a response in
#'  \code{formula}. This is used by \code{\link{generate_data}}, which samples from the prior predictive
#'  distribution.
#' @param control a list with further computational options, see details section.
#' @return A sampler object, which is the main input for the MCMC simulation
#'  function \code{\link{MCMCsim}}. The sampler object is an environment with
#'  precomputed quantities and functions. The main functions are \code{rprior},
#'  which returns a sample from the prior distributions, \code{draw},
#'  which returns a sample from the full conditional posterior distributions,
#'  and \code{start}, which returns a list with starting values for the Gibbs
#'  sampler. If \code{prior.only} is \code{TRUE}, functions \code{draw} and
#'  \code{start} are not created.
#' @references
#'  J.H. Albert and S. Chib (1993).
#'    Bayesian analysis of binary and polychotomous response data.
#'    Journal of the American statistical Association 88(422), 669-679.
#'
#'  D. Bates, M. Maechler, B. Bolker and S.C. Walker (2015).
#'    Fitting Linear Mixed-Effects Models Using lme4.
#'    Journal of Statistical Software 67(1), 1-48.
#'
#'  S.W. Linderman, M.J. Johnson and R.P. Adams (2015).
#'    Dependent multinomial models made easy: Stick-breaking with the Polya-gamma augmentation.
#'    Advances in Neural Information Processing Systems, 3456–3464.
#'
#'  N. Polson, J.G. Scott and J. Windle (2013).
#'    Bayesian Inference for Logistic Models Using Polya-Gamma Latent Variables.
#'    Journal of the American Statistical Association 108(504), 1339-1349.
#'
#'  H. Rue and L. Held (2005).
#'    Gaussian Markov Random Fields.
#'    Chapman & Hall/CRC.
#'    
#'  M. Zhou and L. Carin (2015).
#'    Negative Binomial Process Count and Mixture Modeling.
#'    IEEE Transactions on Pattern Analysis and Machine Intelligence 37(2), 307-320.
create_sampler <- function(formula, data=NULL, family="gaussian",
                           ny=NULL, ry=NULL, r.mod,
                           sigma.fixed=NULL,
                           sigma.mod=NULL, Q0=NULL, formula.V=NULL, logJacobian=NULL,
                           linpred=NULL,
                           compute.weights=FALSE, block=compute.weights,
                           prior.only=FALSE,
                           control=NULL) {

  if (is.character(family)) {
    family <- match.arg(family, c("gaussian", "binomial", "negbinomial", "poisson", "multinomial"))
    family <- eval(call(paste0("f_", family)))
  }

  if (missing(formula)) stop("a formula specifying response and linear predictor model components must be specified")
  if (!inherits(formula, "formula")) stop("'formula' must be a formula")
  formula <- standardize_formula(formula, data=data)  # pass data to interpret '.' in formula
  if (has_response(formula)) {
    y <- get_response(formula, data)
    if (family$family == "multinomial") {
      if (is.matrix(y)) {
        cats <- colnames(y)
      } else {  # categorical/multinoulli
        y <- as.factor(y)
        cats <- levels(y)
        y <- model_matrix(~ 0 + y, sparse=FALSE)
      }
      if (is.null(cats)) cats <- as.character(seq_len(ncol(y)))
      if (is.null(family$K)) {
        family$K <- length(cats)
      } else {
        family$K <- as.integer(family$K)
        if (!identical(family$K, length(cats))) stop("observed number of categories (", length(cats),
          ") does not match specified value K=", family$K)
      }
      Km1 <- family$K - 1L
      if (!is.null(ny)) warning("argument 'ny' ignored for multinomial family")
      # construct ny variable according to stick-breaking representation
      ny <- rowSums(y)
      ny <- c(ny, ny - as.vector(rowCumsums(y[, seq_len(Km1 - 1L), drop=FALSE])))
      y <- as.vector(y[, -length(cats)])
    } else {
      if (is.logical(y)) y <- as.integer(y)
      if (is.character(y)) y <- factor(y)
      if (is.factor(y)) {
        if (nlevels(y) > 2L) stop("response factor variable with more than two levels; please use family='multinomial'")
        y <- as.integer(y) - 1L
      }
      if (!is.numeric(y)) stop("non-numeric response variable")
    }
    if (anyNA(y)) stop(sum(is.na(y)), " missing(s) in response variable")
    n <- length(y)
    if (family$family == "multinomial") {
      n0 <- n %/% Km1
      if (is.null(data)) data <- n0
    } else {
      if (is.null(data)) data <- n
    }
  } else {
    if (!prior.only) {
      warning("no left hand side in 'formula': setting up prior samplers only", immediate.=TRUE)
      prior.only <- TRUE
    }
    if (is.null(data)) {
      if (!length(all.vars(formula))) stop("no way to determine the number of cases")
      n <- length(get(all.vars(formula)[1L]))
      data <- n
    } else
      n <- nrow(data)
    if (family$family == "multinomial") {
      if (is.null(family$K)) stop("for prior predictive multinomial sampling use family=f_multinomial(K) to specify the number of categories")
      Km1 <- family$K - 1L
      cats <- as.character(seq_len(family$K))
      n0 <- n
      n <- n * Km1
    }
  }

  if (n < 1L) stop("empty data")

  mod <- to_mclist(formula)
  if (!length(mod)) stop("empty 'formula'")

  if (family$family == "gaussian") {
    if (!is.null(ny)) warning("argument 'ny' ignored for family 'gaussian'")
    if (prior.only || n == 1L) {
      scale_y <- 1
    } else {
      # scale of y used to generate default starting values
      # divide by number of model components + 1: each component explains part of the total variation
      scale_y <- sd(y) / (length(mod) + 1L)
      if (scale_y == 0 || !is.finite(scale_y)) scale_y <- 1
    }
    if (!isTRUE(sigma.fixed)) sigma.fixed <- FALSE
    f_mean <- identity  # mean function acting on linear predictor
    
    
    
  } else {  # binomial, multinomial, negative-binomial or Poisson likelihood
    if (!is.null(Q0) || !is.null(formula.V)) stop("'Q0' and 'formula.V' arguments cannot be used with (negative) binomial or multinomial likelihood")
    if (identical(sigma.fixed, FALSE) || !is.null(sigma.mod)) warning("arguments 'sigma.fixed' and 'sigma.mod' ignored for (negative) binomial or multinomial likelihood")
    sigma.mod <- NULL
    sigma.fixed <- TRUE
    modeled.Q <- family$link != "probit"
    modeled.r <- FALSE
    scale_y <- 2.5
    if (family$family == "binomial") {
      if (length(ny) == 1L)
        ny.input <- ny  # for use in predict
      else
        ny.input <- NULL
      ny <- check_ny(ny, data)
      if (family$link == "probit") {
        if (!all(ny == 1L)) warning("only binary data with 'ny = 1' supported by binomial probit model")
        ny <- 1L
      } else {
        ny <- as.numeric(ny)  # both CrPGapprox and BayesLogit::rpg require doubles
      }
      f_mean <- function(eta) ny / (1 + exp(-eta))  # mean function acting on linear predictor
      if (!prior.only && any(y > ny)) stop("'y' cannot be larger than 'ny'")
    } else if (family$family == "multinomial") {
      ny <- check_ny(ny, n)
    } else {  # negative binomial or Poisson likelihood
      if (!is.null(ny)) warning("argument 'ny' ignored for negative binomial or Poisson family")
      if (family$family == "negbinomial") {
        if (length(ry) == 1L)
          ry.input <- ry  # for use in predict
        else
          ry.input <- NULL
        ry <- check_ry(ry, data)
        if (missing(r.mod)) r.mod <- pr_invchisq(df=1, scale=1)
        modeled.r <- TRUE
        switch(r.mod$type,
          fixed = {
            if (r.mod$value <= 0) stop("negative binomial dispersion parameter must be positive")
            if (r.mod$value == 1) {
              r.mod <- NULL
              modeled.r <- FALSE
            } else
              r.mod <- pr_fixed(r.mod$value, n=1L, !prior.only)
          },
          invchisq = {
            if (is.list(r.mod$df)) stop("modeled degrees of freedom parameter not supported in 'r.mod'")
            r.mod <- pr_invchisq(r.mod$df, r.mod$scale, n=1L, post = !prior.only)
          },
          gig = {
            r.mod <- pr_gig(r.mod$a, r.mod$b, r.mod$p, n=1L, post = !prior.only)
          },
          stop("unsupported prior")
        )
        if (modeled.r) {
          if (!prior.only) y <- as.numeric(y)  # for C++ rCRT function; alternatively force to be integer
          f_mean <- function(eta) stop("TBI: mean function for negative binomial with modeled dispersion parameter")
        } else {
          if (!prior.only) ny <- y + ry  # used in Polya-Gamma full conditional
          f_mean <- function(eta) ry * exp(eta)  # mean function acting on linear predictor
        }
      } else {  # Poisson family; posterior inference approximated using negative binomial
        if (!is.null(ry)) {
          if (!is.numeric(ry) || length(ry) != 1L || is.na(ry) || ry <= 0) stop("'ry' must be a positive scalar")
        } else {
          # choose a default value large enough for negligible overdispersion, but not so large
          #   that MCMC mixing becomes too slow
          ry <- 100
        }
        if (!prior.only) ny <- y + ry  # used in Polya-Gamma full conditional
        f_mean <- function(eta) ry * exp(eta)  # mean function acting on linear predictor
      }
    }
    if (!prior.only && !family$link == "probit") {
      if (.opts$PG.approx) {
        mPG <- as.integer(.opts$PG.approx.m)
        if (!length(mPG) %in% c(1L, n)) stop("invalid value for option 'PG.approx.m'")
        rPolyaGamma <- function(b, c) CrPGapprox(n, b, c, mPG)
      } else {
        if (!requireNamespace("BayesLogit", quietly=TRUE)) stop("please install package 'BayesLogit' and try again")
        if (family$family == "binomial") {
          ny[ny == 0] <- .Machine$double.eps  # to prevent generation of NA by BayesLogit::rpg
        }
        rpg <- BayesLogit::rpg
        rPolyaGamma <- function(b, c) rpg(n, b, c)
      }
    }
  }

  # data-level variance prior/model
  if (is.null(Q0)) {
    Q0 <- CdiagU(n)
  } else if (is.character(Q0)) {
    Q0 <- Cdiag(data[[Q0]])
  } else if (is.vector(Q0)) {
    if (length(Q0) == 1L)
      Q0 <- Cdiag(rep.int(Q0, n))
    else
      Q0 <- Cdiag(Q0)
  }
  Q0 <- economizeMatrix(Q0, symmetric=TRUE)
  if (!identical(dim(Q0), c(n, n))) stop("incompatible 'Q0'")
  if (isDiagonal(Q0)) {
    if (isUnitDiag(Q0))
      Q0.type <- "unit"
    else
      Q0.type <- "diag"
  } else {
    Q0.type <- "symm"  # non-diagonal precision matrix
  }
  scale_sigma <- scale_y * mean(sqrt(diag(Q0)))

  # data-level variance model
  if (is.null(formula.V)) {
    Vmod <- NULL
    if (family$family == "gaussian") modeled.Q <- FALSE
  } else {
    if (!inherits(formula.V, "formula")) stop("'formula.V' must be a formula")
    formula.V <- standardize_formula(formula.V, c("vreg", "vfac"), "vreg", data=data)  # pass data to interpret '.'
    # TODO deal with offset: divide Q0 by offset (and then set Q0.type)
    Vmod <- to_Vmclist(formula.V)
    if (!length(Vmod)) stop("empty 'formula.V'")
    if (any(names(mod) %in% names(Vmod))) stop("names used in 'formula' and 'formula.V' should be unique")
    modeled.Q <- TRUE
  }

  e.is.res <- family$link == "identity"
  offset <- get_offset(formula, data)
  if (!is.null(offset) || family$family == "poisson") {
    if (family$family == "poisson") {
      # Poisson approximated as negative binomial with large ry and offset -log(ry)
      if (is.null(offset))
        offset <- -log(ry)
      else
        offset <- offset - log(ry)
    }
    if (e.is.res)
      y_eff <- function() y - offset
    else
      y_eff <- function() offset
    use.offset <- TRUE
  } else {
    if (e.is.res)
      y_eff <- function() y
    else
      y_eff <- function() 0
    use.offset <- FALSE
    rm(offset)
  }

  # check Gibbs blocks
  if (is.logical(block)) {
    if (length(block) != 1L) stop("unexpected input for 'block' argument")
    if (block)
      block <- list(names(mod))  # all in a single block
    else
      block <- list()
  } else {
    if (!is.list(block)) stop("'block' should be either a scalar logical or a list of component name vectors")
    if (!length(block)) stop("'block' must contain at least one character vector")
    for (b in seq_along(block)) {
      if (length(block[[b]]) < 2L) warning("block with only one model component")
      if (!all(block[[b]] %in% names(mod))) {
        stop("invalid name(s) '", paste(setdiff(block[[b]], names(mod)), collapse="', '"), "' in block ", b)
      }
    }
    if (any(duplicated(unlist(block)))) stop("a model component name can only appear once in 'block'")
  }
  single.block <- length(mod) %in% c(1L, length(unlist(block)))

  if (!prior.only) {
    if (is.null(control)) control <- list()
    if (!is.list(control)) stop("'control' must be a list")
    control.defaults <- list(add.outer.R = length(block) > 0L, recompute.e=TRUE, CG=NULL)
    if (!all(names(control) %in% names(control.defaults))) stop("unrecognized 'control' parameters")
    control <- modifyList(control.defaults, control)
    rm(control.defaults)
    if (isTRUE(control$CG)) control$CG <- list()  # TODO allow different control$CG for each block
    if (is.list(control$CG) && !length(block))
      warning("conjugate gradient algorithm currently only used inside blocks; see argument 'block'", immediate.=TRUE)

    if (single.block) {  # computation of working response simplifies
      control$recompute.e <- FALSE  # already computed in the coefficient sampling function
      y <- as.numeric(y)  # some C++ matrix algebra functions expect double
    }
  }

  coef.names <- vector(mode="list", length(mod))
  names(coef.names) <- names(mod)
  types <- get_types(formula)
  for (k in seq_along(mod)) {
    mc <- mod[[k]]
    mc$name <- names(mod)[k]
    mc$e <- environment()
    mc <- as.list(mc)[-1L]
    mod[[mc$name]] <- switch(types[k],
      reg = do.call(reg, mc, envir=parent.frame()),
      mec = do.call(mec, mc, envir=parent.frame()),
      gen = do.call(gen, mc, envir=parent.frame())
    )
  }

  # sigma prior
  if (sigma.fixed) {
    if (!is.null(sigma.mod)) warning("'sigma.mod' ignored because 'sigma.fixed=TRUE'")
    sigma.mod <- NULL
  } else {
    if (is.null(sigma.mod)) sigma.mod <- pr_invchisq(df=0, scale=1)
    switch(sigma.mod$type,
      fixed = {
        if (sigma.mod$value <= 0) stop("standard deviation parameter must be positive")
        sigma.mod <- pr_fixed(sigma.mod$value, n=1L, !prior.only)
      },
      invchisq = {
        if (is.list(sigma.mod$df)) stop("modeled degrees of freedom parameter not supported in 'sigma.mod'")
        sigma.mod <- pr_invchisq(sigma.mod$df, sigma.mod$scale, n=1L, !prior.only)
      },
      exp = {
        sigma.mod <- pr_exp(sigma.mod$scale, n=1L, !prior.only)
      },
      gig = {
        sigma.mod <- pr_gig(sigma.mod$a, sigma.mod$b, sigma.mod$p, n=1L, !prior.only)
      },
      stop("unsupported prior")
    )
  }

  # compose following 3 functions:
  # draw: function to draw parameters from their full conditionals
  # rprior: function to draw parameters from their priors
  # start: function to create starting values (including for residuals/linear predictor e_)
  rprior <- function(p) {p <- list()}
  if (family$family == "negbinomial" && modeled.r)
    rprior <- add(rprior, quote(p[["negbin_r_"]] <- 1 / r.mod$rprior()))

  if (modeled.Q && family$family == "gaussian") {
    types <- get_types(formula.V, c("vreg", "vfac"))
    for (k in seq_along(Vmod)) {
      mc <- Vmod[[k]]
      mc$name <- names(Vmod)[k]
      mc$e <- environment()
      mc <- as.list(mc)[-1L]
      switch(types[k],
        vreg = {
          Vmod[[mc$name]] <- do.call(vreg, mc, envir=parent.frame())
          rprior <- add(rprior, bquote(p[[.(mc$name)]] <- Vmod[[.(k)]]$rprior(p)))
        },
        vfac = {
          Vmod[[mc$name]] <- do.call(vfac, mc, envir=parent.frame())
          rprior <- add(rprior, bquote(p <- Vmod[[.(k)]]$rprior(p)))
        }
      )
    }

    # compute product of inverse variance factors
    compute_Qfactor <- function(p) {
      out <- Vmod[[1L]]$compute_Qfactor(p)
      for (mc in Vmod[-1L]) out <- out * mc$compute_Qfactor(p)
      out
    }
    # compute data-level precision matrix from Q0 and scale factor computed by compute_Qfactor
    compute_Q <- function(p, Qfactor=NULL) {
      if (is.null(Qfactor)) Qfactor <- compute_Qfactor(p)
      switch(Q0.type,
        unit = p$Q_ <- Qfactor,
        diag = p$Q_ <- Q0@x * Qfactor,
        symm = {  # store full precision matrix
          p$QM_ <- block_scale_dsCMatrix(Q0, Qfactor)
          p$Q_ <- Qfactor  # by default in store.mean for use in compute_DIC
        }
      )
      p
    }
  }  # END if (modeled.Q && family$family == "gaussian")

  if (!sigma.fixed) rprior <- add(rprior, quote(p$sigma_ <- sqrt(sigma.mod$rprior())))

  for (k in seq_along(mod)) {
    mc <- mod[[k]]
    switch(mc$type,
      reg = rprior <- add(rprior, bquote(p[[.(mc$name)]] <- mod[[.(k)]]$rprior(p))),
      mec = rprior <- add(rprior, bquote(p <- mod[[.(k)]]$rprior(p))),
      gen = rprior <- add(rprior, bquote(p <- mod[[.(k)]]$rprior(p)))
    )
  }

  # build a list of indices for vector representation of all likelihood parameters (for optimization)
  vec_list <- list()
  n_vec_list <- 0L
  # single block sampler assumes coefficients are at the start of the parameter vector
  for (k in seq_along(mod)) {
    vec_list[[names(mod)[k]]] <- (n_vec_list + 1L):(n_vec_list + mod[[k]][["q"]])
    n_vec_list <- n_vec_list + mod[[k]][["q"]]
  }
  if (!sigma.fixed) {
    vec_list[["sigma_"]] <- n_vec_list + 1L
    n_vec_list <- n_vec_list + 1L
  }
  if (!is.null(Vmod)) {
    for (k in seq_along(Vmod)) {
      vec_list[[names(Vmod)[k]]] <- (n_vec_list + 1L):(n_vec_list + Vmod[[k]][["q"]])
      n_vec_list <- n_vec_list + Vmod[[k]][["q"]]
    }
  }
  if (family$family == "negbinomial" && modeled.r) {
    vec_list[["negbin_r_"]] <- n_vec_list + 1L
    n_vec_list <- n_vec_list + 1L
  }

  vec2list <- function(x) {
    pars <- names(vec_list)
    out <- list()
    for (k in seq_along(vec_list)) out[[pars[k]]] <- x[vec_list[[k]]]
    out
  }

  list2vec <- function(p) {
    pars <- names(vec_list)
    out <- rep.int(NA_real_, n_vec_list)
    for (k in seq_along(vec_list)) out[vec_list[[k]]] <- p[[pars[k]]]
    out
  }

  if (is.null(linpred)) {
    do.linpred <- FALSE
  } else {
    if (identical(linpred, "fitted")) {
      linpred <- NULL
    } else {
      if (!all(names(linpred) %in% names(mod))) stop("not all names of 'linpred' are model component names")
      for (k in names(linpred))
        linpred[[k]] <- economizeMatrix(linpred[[k]], strip.names=FALSE)
    }
    do.linpred <- TRUE
    rprior <- add(rprior, quote(p$linpred_ <- lin_predict(p, linpred)))
  }

  lin_predict <- function(p, linpred) {
    out <- if (use.offset) offset else 0
    if (is.null(linpred))
      for (mc in mod) out <- out + mc$linpred(p)
    else
      for (k in names(linpred))
        out <- out + if (is.function(linpred[[k]])) linpred[[k]](p) else linpred[[k]] %m*v% p[[k]]
    out
  }
  # generate from predictive distribution; lp is output of lin_predict
  switch(family$family,
    gaussian =
      rpredictive <- function(p, lp, cholQ=NULL, var=NULL, V=NULL, ny, ry) {
        sigma <- if (sigma.fixed) 1 else p[["sigma_"]]
        if (is.null(var)) {  # in-sample, cholQ must be supplied
          if (modeled.Q) {
            if (is.null(p[["Q_"]])) p <- compute_Q(p)
            if (Q0.type == "symm")
              cholQ$update(p[["QM_"]])
            else
              cholQ$update(p[["Q_"]])
          }
          lp + drawMVN_cholQ(cholQ, sd=sigma)
        } else {
          if (modeled.Q)
            for (mc in Vmod) var <- var * V[[mc$name]](p)
          lp + sigma * sqrt(var) * Crnorm(length(lp))
        }
      },
    binomial =
      rpredictive <- function(p, lp, cholQ=NULL, var=NULL, V=NULL, ny, ry)
        rbinom(length(lp), size=ny, prob=family$linkinv(lp)),
    multinomial = {
      # long vector format
      rpredictive <- function(p, lp, cholQ=NULL, var=NULL, V=NULL, ny, ry) {
        ptilde <- family$linkinv(lp)
        out <- integer(length(lp))
        n0 <- length(lp) %/% Km1
        ind <- seq_len(n0)
        # assume ny has length n0 or 1
        temp <- rbinom(n0, size=ny, prob=ptilde[ind])
        out[ind] <- temp
        for (k in 2:Km1) {
          ny <- ny - temp
          ind <- ind + n0
          temp <- rbinom(n0, size=ny, prob=ptilde[ind])
          out[ind] <- temp
        }
        out
      }
      # short vector format, only possible for categorica/multinoulli data
      rpredictive_cat <- function(p, lp, cholQ=NULL, var=NULL, V=NULL, ny, ry) {
        pSB <- family$linkinv(lp)
        n0 <- length(lp) %/% Km1
        out <- rep.int(family$K, n0)  # baseline category
        out[ny == 0L] <- NA_integer_
        ind <- seq_len(n0)
        # assume ny is 1 or 0
        temp <- rbinom(n0, size=ny, prob=pSB[ind])
        out[temp == 1L] <- 1L
        for (k in 2:Km1) {
          ny <- ny - temp
          ind <- ind + n0
          temp <- rbinom(n0, size=ny, prob=pSB[ind])
          out[temp == 1L] <- k
        }
        out
      }
    },
    negbinomial =
      # NB definition of rnbinom has p <-> 1-p
      if (modeled.r) {
        rpredictive <- function(p, lp, cholQ=NULL, var=NULL, V=NULL, ny, ry)
          rnbinom(length(lp), size=ry*p[["negbin_r_"]], prob=1/(1 + exp(lp)))
      } else {
        rpredictive <- function(p, lp, cholQ=NULL, var=NULL, V=NULL, ny, ry)
          rnbinom(length(lp), size=ry, prob=1/(1 + exp(lp)))
      },
    poisson =
      rpredictive <- function(p, lp, cholQ=NULL, var=NULL, V=NULL, ny, ry)
        rpois(length(lp), lambda=exp(lp))
  )

  rprior <- add(rprior, quote(p))  # return state p

  # TODO logprior for optimization of logposterior = logprior + llh
  log_prior <- function(p) {
    out <- 0
    sigma <- if (sigma.fixed) 1 else p[["sigma_"]]
    for (mc in mod)
      out <- out + mc$log_prior_0 +
        switch(mc$type,
          reg = {
            delta.beta <- p[[mc$name]] - mc$b0
            - mc$q * log(sigma) - 0.5 * dotprodC(delta.beta, mc$Q0 %m*v% delta.beta) / sigma^2
          },
          mec = {
            stop("TBI: log-prior for measurement error in covariates component")
          },
          gen = {
            stop("TBI: log-prior for generic model component")
          }
        )
    # TODO add log-priors of Vmod components
    out
  }

  # What parameters to store by default? This information is used by MCMCsim.
  store_default <- function(prior.sampler=FALSE) {
    out <- if (prior.sampler) NULL else "llh_"
    if (!sigma.fixed) out <- c(out, "sigma_")
    if (do.linpred) out <- c(out, "linpred_")
    if (family$family == "negbinomial" && modeled.r) out <- c(out, "negbin_r_")
    for (mc in mod)
      switch(mc$type,
        reg = out <- c(out, mc$name),
        mec = out <- c(out, mc$name),
        gen = {
          out <- c(out, mc$name_sigma)
          if (mc$var == "unstructured") out <- c(out, mc$name_rho)
          if (mc$gl) out <- c(out, mc$name_gl)
          if (mc$Leroux_update) out <- c(out, mc$name_Leroux)
          if (!is.null(mc$priorA) && is.list(mc$priorA$df)) out <- c(out, mc$name_df)
        }
      )
    for (mc in Vmod)
      switch(mc$type,
        vreg = out <- c(out, mc$name),
        vfac = if (mc$prior$type == "invchisq" && is.list(mc$prior$df)) out <- c(out, mc$name_df)
      )
    out
  }
  store_mean_default <- function(prior.sampler=FALSE) {
    if (prior.sampler) {
      out <- NULL
    } else {
      out <- "e_"
      if (family$family == "gaussian" && modeled.Q) out <- c(out, "Q_")
      if (compute.weights) out <- c(out, "weights_")
    }
    out
  }

  rm(types, mc, k)

  if (prior.only) {
    return(environment())  # environment including rprior(), but not draw() and start()
  }


  if (compute.weights && !do.linpred) stop("weights can only be computed for a linear predictor specified by argument 'linpred'")
  if (compute.weights && !single.block) stop("'compute.weights=TRUE' cannot be combined with 'block=FALSE'")

  if (length(block)) {
    mbs <- list()
    for (b in seq_along(block)) mbs[[b]] <- create_mc_block(mod[block[[b]]])
    rm(b)
  }

  # compute residuals/linear predictor e_ for state p
  if (e.is.res) {
    # residuals
    compute_e <- function(p) {
      if (use.offset)
        e <- y - offset
      else
        e <- y
      for (mc in mod) e <- e - mc$linpred(p)
      e
    }
  } else {
    # linear predictor
    compute_e <- function(p) {
      e <- mod[[1L]]$linpred(p)
      for (mc in mod[-1L]) e <- e + mc$linpred(p)
      if (use.offset) e + offset else e
    }
  }

  # function Q_e computes matrix-vector product of Q and working response e_
  if (family$family == "gaussian") {
    if (modeled.Q) {
      Q_e <- switch(Q0.type,
        unit=, diag = function(p) p[["Q_"]] * p[["e_"]],
        # in case of non-diagonal Q0, full precision matrix stored as p[["QM_"]]
        symm = function(p) p[["QM_"]] %m*v% p[["e_"]]
      )
    } else {
      Q_e <- switch(Q0.type,
        unit = function(p) p[["e_"]],
        diag = function(p) Q0@x * p[["e_"]],
        symm = function(p) Q0 %m*v% p[["e_"]]
      )
    }
  } else {  # binomial or negative binomial
    # in this case Q0 is ignored
    if (family$link == "probit") {
      if (single.block && !use.offset)
        Q_e <- function(p) p[["z_"]]
      else
        Q_e <- function(p) p[["z_"]] - p[["e_"]]
    } else if (family$family == "negbinomial" && modeled.r) {
      if (single.block && !use.offset)
        Q_e <- function(p) 0.5 * (y - ry*p[["negbin_r_"]])
      else
        Q_e <- function(p) 0.5 * (y - ry*p[["negbin_r_"]]) - p[["Q_"]] * p[["e_"]]
    } else {
      if (family$family %in% c("negbinomial", "poisson"))
        y_shifted <- 0.5 * (y - ry)
      else  # binomial or multinomial
        y_shifted <- y - 0.5 * ny
      if (single.block && !use.offset)
        Q_e <- function(p) y_shifted
      else
        Q_e <- function(p) y_shifted - p[["Q_"]] * p[["e_"]]
    }
  }

  if (control$recompute.e) {
    # compute residuals/linear predictor at each call of sampler (reduces build-up of rounding error?)
    draw <- function(p) {p$e_ <- compute_e(p)}
  } else {
    draw <- function(p) {}
  }
  start <- function(p=list()) {
    if (!is.list(p)) stop("input to 'start' function must be a list")
    if (is.null(p[["e_"]])) p$e_ <- Crnorm(n, sd=scale_y)
    if (length(p[["e_"]]) != n) stop("wrong length for 'e_' start value")
  }
  MHpars <- NULL
  adapt <- function(ar) {}

  if (family$family == "negbinomial" && modeled.r) {
    if (r.mod$type == "fixed") {
      start <- add(start, quote(if (is.null(p[["negbin_r_"]])) p$negbin_r_ <- 1/r.mod$rprior()))
    } else {
      # draw latent L_i (i=1:n) from its CRT f.c.
      mCRT <- .opts$CRT.approx.m
      draw <- add(draw, quote(L <- CrCRT(y, ry*p[["negbin_r_"]], mCRT)))
      start <- add(start, quote(
        if (is.null(p[["negbin_r_"]])) p$negbin_r_ <- runif(1L, 0.1, 10)
        else {
          if (length(p[["negbin_r_"]]) != 1L) stop("wrong length for 'negbin_r_' start value")
          p[["negbin_r_"]] <- as.numeric(p[["negbin_r_"]])
          if (is.na(p[["negbin_r_"]]) || p[["negbin_r_"]] <= 0) stop("'negbin_r_' start value must be a positive number")
        }
      ))
    }
    # draw dispersion parameter r
    switch(r.mod$type,
      fixed = {
        draw <- add(draw, quote(p$negbin_r_ <- 1/r.mod$draw()))
      },
      invchisq = {
        # TODO check ry appearance in modeled scale inv-chi2 and gig cases
        if (is.list(r.mod$scale))
          if (length(ry) == 1L)
            draw <- add(draw, quote(p$negbin_r_ <- 1 / r.mod$draw(2*sum(L), 2*ry*sum(log1pexpC(p[["e_"]])), p[["negbin_r_"]])))
          else
            draw <- add(draw, quote(p$negbin_r_ <- 1 / r.mod$draw(2*sum(L), 2*sum(ry*log1pexpC(p[["e_"]])), p[["negbin_r_"]])))
        else
          if (length(ry) == 1L)
            draw <- add(draw, quote(p$negbin_r_ <- 1 / r.mod$draw(2*sum(L), 2*ry*sum(log1pexpC(p[["e_"]])))))
          else
            draw <- add(draw, quote(p$negbin_r_ <- 1 / r.mod$draw(2*sum(L), 2*sum(ry*log1pexpC(p[["e_"]])))))
      },
      gig = {
        if (length(ry) == 1L)
          draw <- add(draw, quote(p$negbin_r_ <- 1 / r.mod$draw(p = r.mod$p - sum(L), a = r.mod$a, b = r.mod$b + 2*ry*sum(log1pexpC(p[["e_"]])))))
        else
          draw <- add(draw, quote(p$negbin_r_ <- 1 / r.mod$draw(p = r.mod$p - sum(L), a = r.mod$a, b = r.mod$b + 2*sum(ry*log1pexpC(p[["e_"]])))))
      }
    )
    draw <- add(draw, quote(ny <- y + ry*p[["negbin_r_"]]))  # used in Polya-Gamma full conditional for latent precision vector
    start <- add(start, quote(ny <- y + ry*p[["negbin_r_"]]))
  }

  if (family$family %in% c("binomial", "multinomial", "negbinomial", "poisson") && family$link != "probit") {
    draw <- add(draw, quote(p[["Q_"]] <- rPolyaGamma(ny, p[["e_"]])))
    draw <- add(draw, quote(p$llh_ <- llh(p)))
    start <- add(start, quote(if (is.null(p[["Q_"]])) p$Q_ <- rPolyaGamma(ny, p[["e_"]])))
    start <- add(start, quote(if (length(p[["Q_"]]) != n) stop("wrong length for 'Q_' start value")))
  }

  if (family$family == "binomial" && family$link == "probit") {
    draw <- add(draw, quote(p[["z_"]] <- CrTNprobit(p[["e_"]], y)))
    draw <- add(draw, quote(p$llh_ <- llh(p)))
    start <- add(start, quote(if (is.null(p[["z_"]])) p$z_ <- CrTNprobit(p[["e_"]], y)))
    start <- add(start, quote(if (length(p[["z_"]]) != n) stop("wrong length for 'z_' start value")))
  }

  if (family$family == "gaussian") {
    if (modeled.Q) {
      if (length(Vmod) > 1L) {
        # recompute precision Q for numerical stability
        draw <- add(draw, quote(p <- compute_Q(p)))
      }
      for (mc in Vmod) {
        switch(mc$type, 
          vreg = MHpars <- c(MHpars, mc$name),
          vfac = {
            if (mc$prior$type == "invchisq" && is.list(mc$prior$df)) {
              MHpars <- c(MHpars, mc$name_df)
              if (mc$prior$df$adapt)
                adapt <- add(adapt, bquote(Vmod[[.(mc$name)]]$adapt(ar)))
            }
          }
        )
        draw <- add(draw, bquote(p <- Vmod[[.(mc$name)]]$draw(p)))
        start <- add(start, bquote(p <- Vmod[[.(mc$name)]]$start(p)))
      }
      start <- add(start, quote(p <- compute_Q(p)))
    }  # END if (modeled.Q)
    draw <- add(draw, quote(SSR <- dotprodC(p[["e_"]], Q_e(p))))
    draw <- add(draw, quote(p$llh_ <- llh(p, SSR)))
  }  # END if (family$family == "gaussian")

  if (!sigma.fixed) {
    draw_sigma <- function(p, SSR) {}
    # compute df.data + update SSR with contributions from reg and gen components
    df.data <- n
    for (k in seq_along(mod)) {
      mc <- mod[[k]]
      switch(mc$type,
        reg=, mec = {
          if (mc$informative.prior) {
            if (is.null(mc$R))
              df.data <- df.data + mc$q
            else
              df.data <- df.data + mc$q - ncol(mc$R)
            draw_sigma <- add(draw_sigma, bquote(delta.beta <- p[[.(mc$name)]] - mod[[.(k)]]$b0))
            draw_sigma <- add(draw_sigma, bquote(SSR <- SSR + dotprodC(delta.beta, mod[[.(k)]]$Q0 %m*v% delta.beta)))
          }
        },
        gen = {
          if (mc$usePX && mc$PX$data.scale) {
            df.data <- df.data + if (mc$PX$vector) mc$q0 else 1L
            draw_sigma <- add(draw_sigma, bquote(SSR <- SSR + dotprodC(p[[.(mc$name_xi)]], mod[[.(k)]]$PX$Q0 %m*v% p[[.(mc$name_xi)]])))
          }
          if (mc$gl && mc$glp$informative.prior) {
            df.data <- df.data + mc$glp$q
            draw_sigma <- add(draw_sigma, bquote(delta.beta <- p[[.(mc$name_gl)]] - mod[[.(k)]]$glp$b0))
            draw_sigma <- add(draw_sigma, bquote(SSR <- SSR + dotprodC(delta.beta, mod[[.(k)]]$glp$Q0 %m*v% delta.beta)))
          }
        },
        stop("unrecognized model type")
      )
    }
    switch(sigma.mod$type,
      fixed = {
        draw_sigma <- add(draw_sigma, quote(p$sigma_ <- sqrt(sigma.mod$draw())))
      },
      invchisq = {
        if (is.list(sigma.mod$scale))
          draw_sigma <- add(draw_sigma, quote(p$sigma_ <- sqrt(sigma.mod$draw(df.data, SSR, 1 / p[["sigma_"]]^2))))
        else
          draw_sigma <- add(draw_sigma, quote(p$sigma_ <- sqrt(sigma.mod$draw(df.data, SSR))))
      },
      exp = {
        draw_sigma <- add(draw_sigma, quote(p$sigma_ <- sqrt(sigma.mod$draw(1 - 0.5*df.data, 2/sigma.mod$scale, SSR))))
      },
      gig = {
        draw_sigma <- add(draw_sigma, quote(p$sigma_ <- sqrt(sigma.mod$draw(sigma.mod$p - 0.5*df.data, sigma.mod$a, sigma.mod$b + SSR))))
      }
    )
    draw_sigma <- add(draw_sigma, quote(p))

    draw <- add(draw, quote(p <- draw_sigma(p, SSR)))

    if (sigma.mod$type == "fixed") {
      start <- add(start, quote(if (is.null(p[["sigma_"]])) p$sigma_ <- sqrt(sigma.mod$rprior())))
    } else {
      #start <- add(start, quote(if (is.null(p[["sigma_"]])) p$sigma_ <- runif(1L, 0.1 * scale_y, scale_y) * mean(sqrt(diag(Q0)))))
      start <- add(start, quote(if (is.null(p[["sigma_"]])) p$sigma_ <- runif(1L, 0.1 * scale_sigma, scale_sigma)))
      start <- add(start, quote(if (length(p[["sigma_"]]) != 1L) stop("wrong length for 'sigma_' start value")))
    }
  }  # END if (!sigma.fixed)

  for (k in seq_along(mod)) {
    mc <- mod[[k]]
    switch(mc$type,
      reg = {
        if (!mc$in_block) {
          start <- add(start, bquote(p <- mod[[.(k)]]$start(p)))
          draw <- add(draw, bquote(p <- mod[[.(k)]]$draw(p)))
        }
      },
      mec = {
        start <- add(start, bquote(p <- mod[[.(k)]]$start(p)))
        draw <- add(draw, bquote(p <- mod[[.(k)]]$draw(p)))
      },
      gen = {
        if (mc$gl && mc$usePX) MHpars <- c(MHpars, mc$name_xi)
        if (mc$Leroux_update) MHpars <- c(MHpars, mc$name_Leroux)
        start <- add(start, bquote(p <- mod[[.(k)]]$start(p)))
        draw <- add(draw, bquote(p <- mod[[.(k)]]$draw(p)))
      }
    )
  }

  if (length(block)) {
    for (k in seq_along(mbs)) {
      draw <- add(draw, bquote(p <- mbs[[.(k)]]$draw(p)))
      start <- add(start, bquote(p <- mbs[[.(k)]]$start(p)))
    }
  }

  if (do.linpred) {
    draw <- add(draw, quote(p$linpred_ <- lin_predict(p, linpred)))
  }

  draw <- add(draw, quote(p))  # return state p

  if (!control$recompute.e && !single.block) {
    # adding this sometimes gives bad starting values in case of single.block (no need anyway in that case)
    start <- add(start, quote(p$e_ <- compute_e(p)))
  }
  start <- add(start, quote(p))
  if (length(body(adapt)) <= 1L) rm(adapt)

  # data and functions for log-likelihood and related computations (DIC, WAIC)
  # llh can also be used as an objective function to optimize
  if (!is.null(logJacobian)) {
    if (family$family != "gaussian") warning("argument 'logJacobian' ignored for non-gaussian distributions")
    logJacobian <- as.numeric(logJacobian)
    if (length(logJacobian) != n) stop("'logJacobian' should be a vector of the same size as the response vector")
    if (anyNA(logJacobian)) stop("missing values in 'logJacobian'")
  }
  switch(family$family,
    gaussian = {
      # constant term llh_0 of log-likelihood
      llh_0 <- -0.5 * n * log(2*pi)
      if (!is.null(logJacobian)) llh_0 <- llh_0 + sum(logJacobian)
      if (!modeled.Q || Q0.type == "symm") {
        llh_0 <- llh_0 + 0.5 * as.numeric(determinant(1 * Q0, logarithm=TRUE)$modulus)  # 1 * Q0 to ensure no chol object is stored with Q0 and possibly other refs to Q0
      }
      # log-likelihood function; SSR need not be recomputed if it is provided
      llh <- function(p, SSR=NULL) {
        out <- llh_0
        if (modeled.Q) out <- out + 0.5 * sum(log(p[["Q_"]]))
        if (is.null(SSR)) {
          if (is.null(p[["e_"]])) p$e_ <- compute_e(p)  # compute residuals (for DIC)
          SSR <- dotprodC(p[["e_"]], Q_e(p))
        }
        if (sigma.fixed)
          out - 0.5 * SSR
        else
          out - n * log(p[["sigma_"]]) - 0.5 * SSR / p[["sigma_"]]^2
      }
      # for WAIC computation: compute log-likelihood for each observation/batch of observations, vectorized over parameter draws and observations
      llh_0_i <- -0.5 * log(2*pi)
      llh_i <- function(draws, i=seq_len(n)) {
        nc <- nchains(draws)
        ni <- ndraws(draws)
        all.units <- length(i) == n
        if (is.null(draws[["e_"]])) {
          res_i <- matrix(rep_each(y_eff()[i], nc*ni), nc*ni, length(i))
          for (mc in mod) {
            if (all.units) Xi <- mc$X else Xi <- mc$X[i, , drop=FALSE]
            res_i <- res_i - tcrossprod(as.matrix.dc(draws[[mc$name]], colnames=FALSE), Xi)
          }
        } else {
          res_i <- as.matrix.dc(get_from(draws[["e_"]], vars=i))
        }
        if (sigma.fixed) {
          q <- 1
        } else {
          sigma <- as.numeric(as.matrix.dc(draws[["sigma_"]], colnames=FALSE))
          q <- 1/sigma^2
        }
        if (modeled.Q) {
          q <- matrix(q, nc*ni, length(i))
          for (mc in Vmod) {
            if (all.units) Xi <- mc$X else Xi <- mc$X[i, , drop=FALSE]
            q <- q * switch(mc$type,
              vreg = exp(-tcrossprod(as.matrix.dc(draws[[mc$name]], colnames=FALSE), Xi)),
              vfac = tcrossprod(1 / as.matrix.dc(draws[[mc$name]], colnames=FALSE), Xi)
            )
          }
        }
        switch(Q0.type,
          diag = q <- q * rep_each(Q0@x[i], each=nc*ni),
          symm = {
            dQ0 <- rep_each(diag(Q0)[i], each=nc*ni)
            q <- q * dQ0
            if (length(i) == n) {
              # pointwise loo-llh for non-factorizable models: p(y_i|y_{-i},theta)
              # if length(i) < n we currently ignore this correction; should warn in the calling function
              res_i <- (1 / dQ0) * (res_i %m*m% Q0)
            }
          }
        )
        if (is.null(logJacobian))
          llh_0_i + 0.5 * ( log(q) - q * res_i^2 )
        else
          llh_0_i + matrix(rep_each(logJacobian[i], nc*ni), nc*ni, length(i)) + 0.5 * ( log(q) - q * res_i^2 )
      }
    },
    binomial=, multinomial=, negbinomial=, poisson = {
      # NB for Poisson we use the actual negative binomial likelihood used to approximate Poisson
      if (family$family %in% c("negbinomial", "poisson")) {
        if (modeled.r)
          llh_0 <- - sum(lgamma(y + 1))
        else
          llh_0 <- sum(negbinomial_coef(ry, y))
      } else {
        llh_0 <- sum(binomial_coef(ny, y))  # zero in case of binary data
      }
      # argument SSR only for gaussian llh, but needed here as well to prevent check NOTE
      llh <- function(p, SSR=NULL) {
        if (is.null(p[["e_"]])) p$e_ <- compute_e(p)  # for use in DIC
        neg_fitted <- -p[["e_"]]
        if (modeled.r) {  # negative binomial with unknown dispersion parameter
          r <- ry * p[["negbin_r_"]]
          ny <- y + r
          llh_0 + sum(lgamma(ny) - lgamma(r)) + sum(r * neg_fitted - ny * log1pexpC(neg_fitted))
        } else {
          if (family$link == "probit") {
            sum(pnorm((1 - 2*y) * neg_fitted, log.p=TRUE))
          } else {  # logistic binomial/multinomial or negative binomial model
            # faster than dbinom since normalization constant computed only once
            llh_0 + sum((ny - y) * neg_fitted - ny * log1pexpC(neg_fitted))
            #p_eta <- 1 / (1 + exp(neg_fitted))
            #sum(dbinom(y, ny, p_eta, log=TRUE))
          }
        }
      }
      llh_i <- function(draws, i=seq_len(n)) {
        nc <- nchains(draws)
        ni <- ndraws(draws)
        nr <- nc*ni
        all.units <- length(i) == n
        if (is.null(draws[["e_"]])) {
          mc <- mod[[1L]]
          if (all.units) Xi <- mc$X else Xi <- mc$X[i, , drop=FALSE]
          neg_fitted_i <- -tcrossprod(as.matrix.dc(draws[[mc$name]], colnames=FALSE), Xi)
          for (mc in mod[-1L]) {
            if (all.units) Xi <- mc$X else Xi <- mc$X[i, , drop=FALSE]
            neg_fitted_i <- neg_fitted_i - tcrossprod(as.matrix.dc(draws[[mc$name]], colnames=FALSE), Xi)
          }
          if (use.offset)
            if (length(offset) == 1L)
              neg_fitted_i <- neg_fitted_i - offset
            else
              neg_fitted_i <- neg_fitted_i - rep_each(offset[i], nr)
        } else {
          neg_fitted_i <- -as.matrix.dc(get_from(draws[["e_"]], vars=i))
        }

        if (family$family %in% c("negbinomial", "poisson")) {
          if (modeled.r) {
            r <- as.numeric(as.matrix.dc(draws[["negbin_r_"]], colnames=FALSE))
            if (length(ry) == 1L) {
              if (ry != 1) r <- ry * r
            } else {
              r <- r * rep_each(ry[i], nr)
            }
            yi <- rep_each(y[i], nr)
            nyi <- yi + r
            negbinomial_coef(r, yi) + r * neg_fitted_i - nyi * log1pexpC(neg_fitted_i)
          } else {
            if (length(ry) == 1L)
              rep_each(negbinomial_coef(ry, y[i]), nr) + ry * neg_fitted_i - rep_each(ny[i], nr) * log1pexpC(neg_fitted_i)
            else
              rep_each(negbinomial_coef(ry[i], y[i]), nr) + rep_each(ry[i], nr) * neg_fitted_i - rep_each(ny[i], nr) * log1pexpC(neg_fitted_i)
          }
        } else {  # binomial or multinomial likelihood
          if (family$link == "logit") {
            if (length(ny) == 1L)  # typically binary regression, ny=1
              rep_each(binomial_coef(ny, y[i]), nr) + rep_each(ny - y[i], nr) * neg_fitted_i - ny * log1pexpC(neg_fitted_i)
            else  # general (negative) binomial regression, ny has length n
              rep_each(binomial_coef(ny[i], y[i]), nr) + rep_each(ny[i] - y[i], nr) * neg_fitted_i - rep_each(ny[i], nr) * log1pexpC(neg_fitted_i)
          } else {  # probit
            pnorm(rep_each(1 - 2*y[i], nr) * neg_fitted_i, log.p=TRUE)
          }
        }
      }
    }
  )  # END switch(family$family,

  # log-likelihood function for optimization: function of a vector instead of list
  llh_opt <- function(x) {
    p <- vec2list(x)
    llh(p)
  }

  # remove quantities no longer needed from the closure
  rm(k, mc, data)

  # return the function environment, including draw, rprior, start functions
  environment()
}
