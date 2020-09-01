#' Create a model component object for a generic random effects component in the linear predictor
#'
#' This function is intended to be used on the right hand side of the \code{formula} argument to
#' \code{\link{create_sampler}} or \code{\link{generate_data}}.
#'
#' @export
#' @param formula a model formula specifying the effects that vary over the levels of
#'  the factor variable(s) specified by argument \code{factor}. Defaults to \code{~1},
#'  corresponding to random intercepts. If \code{X} is specified \code{formula} is ignored.
#'  Variable names are looked up in the data frame passed as \code{data} argument to
#'  \code{\link{create_sampler}} or \code{\link{generate_data}}, or in \code{environment(formula)}.
#' @param factor a formula with factors by which the effects specified in the \code{formula}
#'  argument vary. Often only one such factor is needed but multiple factors are allowed so that
#'  interaction terms can be modeled conveniently. The formula must take the form
#'  \code{ ~ f1(fac1, ...) * f2(fac2, ...) ...}, where
#'  \code{fac1, fac2} are factor variables and \code{f1, f2} determine the
#'  correlation structure assumed between levels of each factor, and the \code{...} indicate
#'  that for some correlation types further arguments can be passed. Correlation structures
#'  currently supported include \code{iid} for independent identically distributed effects,
#'  \code{RW1} and \code{RW2} for random walks of first or second order over the factor levels,
#'  \code{AR1} for first-order autoregressive effects, \code{season} for seasonal effects,
#'  \code{spatial} for spatial (CAR) effects and \code{custom} for supplying a custom
#'  precision matrix corresponding to the levels of the factor. For further details about
#'  the correlation structures, and further arguments that can be passed, see
#'  \code{\link{correlation}}. Argument \code{factor} is ignored if \code{X} is specified.
#'  The factor variables are looked up in the data frame passed as \code{data} argument to
#'  \code{\link{create_sampler}} or \code{\link{generate_data}}, or in \code{environment(formula)}.
#' @param remove.redundant whether redundant columns should be removed from the model matrix
#'  associated with \code{formula}. Default is \code{FALSE}.
#' @param drop.empty.levels whether to remove factor levels without observations.
#' @param X A (possibly sparse) design matrix. This can be used instead of \code{formula} and \code{factor}.
#' @param var the (co)variance structure among the varying effects defined by \code{formula}
#'  over the levels of the factors defined by \code{factor}.
#'  The default is \code{"unstructured"}, meaning that a full covariance matrix
#'  parameterization is used. For uncorrelated effects with different variances use
#'  \code{var="diagonal"}. For uncorrelated and equal variances use \code{var="scalar"}.
#'  In the case of a single varying effect there is no difference between these choices.
#' @param prior the prior specification for the variance parameters of the random effects.
#'  These can currently be specified by a call to \code{\link{pr_invwishart}} in case
#'  \code{var="unstructured"} or by a call to \code{\link{pr_invchisq}} otherwise.
#'  See the documentation of those prior specification functions for more details.
#' @param Q0 precision matrix associated with \code{formula}. This can only be used in
#'  combination with \code{var="scalar"}.
#' @param PX whether parameter expansion should be used. Default is \code{TRUE}, which
#'  applies parameter expansion with default options. Alternative options can be specified
#'  by supplying a list with one or more of the following components:
#'  \describe{
#'    \item{vector}{whether a redundant multiplicative expansion parameter is used for each varying effect
#'      specified by \code{formula}. The default is \code{TRUE} except when \code{var="scalar"}.
#'      If \code{FALSE} a single redundant multiplicative parameter is used.}
#'    \item{data.scale}{whether the data level scale is used as a variance factor for the expansion
#'      parameters. Default is \code{TRUE}.}
#'    \item{mu0}{location (vector) parameter for parameter expansion. By default \code{0}.}
#'    \item{Q0}{precision (matrix) parameter for parameter expansion. Default is the identity matrix.}
#     \item{sparse}{UNDOCUMENTED}
#'  }
#' @param GMRFmats list of incidence/precision/constraint matrices. This can be specified
#'  as an alternative to \code{factor}. It should be a list such as is usually returned
#'  by \code{\link{compute_GMRF_matrices}}. Can be used together with argument \code{X}
#'  as a flexible alternative to \code{formula} and \code{factor}.
#' @param priorA prior distribution for scale factors at the variance scale associated with \code{QA}.
#'  In case of IGMRF models the scale factors correspond to the innovations.
#'  The default \code{NULL} means not to use any local scale factors.
#'  A prior can currently be specified using \code{\link{pr_invchisq}} or \code{\link{pr_exp}}.
#' @param Leroux this option alters the precision matrix determined by \code{factor} by taking a
#'  weighted average of it with the identity matrix. If \code{TRUE} the model gains an additional parameter,
#'  the 'Leroux' parameter, being the weight of the original, structured, precision matrix in the weighted
#'  average. By default a uniform prior for the weight and a uniform Metropolis-Hastings proposal density
#'  are employed. This default can be changed by supplying a list with elements a, b, and a.star, b.star,
#'  implying a beta(a, b) prior and a beta(a.star, b.star) independence proposal density. A third option is
#'  to supply a single number between 0 and 1, which is then used as a fixed value for the Leroux parameter.
#' @param R0 an optional equality restriction matrix acting on the coefficients defined by \code{formula}, for each
#'  level defined by \code{factor}. If c is the number of restrictions, \code{R0} is a
#'  q0 x c matrix where q0 is the number of columns of the design matrix derived
#'  from \code{formula}. Together with \code{RA} it defines the set of equality constraints
#'  to be imposed on the vector of coefficients. Only allowed in combination with \code{var="scalar"}.
#' @param RA an optional equality restriction matrix acting on the coefficients defined by \code{factor},
#'  for each effect defined by \code{formula}. If c is the number of restrictions, \code{RA} is a
#'  l x c matrix where l is the number of levels defined by \code{factor}.
#'  Together with \code{R0} this defines the set of equality constraints to be imposed on the vector
#'  of coefficients.
#'  If \code{constr=TRUE}, additional constraints are imposed, corresponding to the
#'  null-vectors of the singular precision matrix in case of an intrinsic Gaussian Markov Random Field.
#' @param constr whether constraints corresponding to the null-vectors of the precision matrix
#'  are to be imposed on the vector of coefficients. By default this is \code{TRUE} for
#'  improper or intrinsic GMRF model components, i.e. components with a singular precision matrix
#'  such as random walks or CAR spatial components.
#' @param S0 an optional inequality restriction matrix acting on the coefficients defined by \code{formula}, for each
#'  level defined by \code{factor}. If c is the number of restrictions, \code{S0} is a
#'  q0 x c matrix where q0 is the number of columns of the design matrix derived
#'  from \code{formula}. Together with \code{SA} it defines the set of inequality constraints
#'  to be imposed on the vector of coefficients.
#   TODO does this work, also for var != "scalar"?
#' @param SA an optional inequality restriction matrix acting on the coefficients defined by \code{factor},
#'  for each effect defined by \code{formula}. If c is the number of restrictions, \code{SA} is a
#'  l x c matrix where l is the number of levels defined by \code{factor}.
#'  Together with \code{S0} this defines the set of constraints to be imposed on the vector
#'  of coefficients.
# TODO R,S instead of R0,RA and S0,SA, and rhs r,s
#' @param formula.gl a formula of the form \code{~ glreg(...)} for group-level predictors
#'  around which the random effect component is hierarchically centered.
#'  See \code{\link{glreg}} for details.
#' @param name the name of the model component. This name is used in the output of the MCMC simulation
#'  function \code{\link{MCMCsim}}. By default the name will be 'gen' with the number of the model term attached.
#' @param sparse whether the model matrix associated with \code{formula} should be sparse.
#'  The default is based on a simple heuristic based on storage size.
#' @param perm whether permutation should be used in the Cholesky decomposition used for updating
#'  the model component's coefficient vector. Default is based on a simple heuristic.
#' @param debug if \code{TRUE} a breakpoint is set at the beginning of the posterior
#'  draw function associated with this model component. Mainly intended for developers.
#' @param e for internal use only.
#' @return An object with precomputed quantities and functions for sampling from
#'  prior or conditional posterior distributions for this model component. Only intended
#'  for internal use by other package functions.
#' @references
#'  J. Besag and C. Kooperberg (1995).
#'    On Conditional and Intrinsic Autoregression.
#'    Biometrika 82(4), 733-746.
#'
#'  C.M. Carvalho, N.G. Polson and J.G. Scott (2010).
#'    The horseshoe estimator for sparse signals.
#'    Biometrika 97(2), 465-480.
#'
#'  L. Fahrmeir, T. Kneib and S. Lang (2004).
#'    Penalized Structured Additive Regression for Space-Time Data:
#'    a Bayesian Perspective.
#'    Statistica Sinica 14, 731-761.
#'
#'  A. Gelman (2006).
#'    Prior distributions for variance parameters in hierarchical models.
#'    Bayesian Analysis 1(3), 515-533.
#'
#'  A. Gelman, D.A. Van Dyk, Z. Huang and W.J. Boscardin (2008).
#'    Using Redundant Parameterizations to Fit Hierarchical Models.
#'    Journal of Computational and Graphical Statistics 17(1), 95-122.
#'
#'  B. Leroux, X. Lei and N. Breslow (1999).
#'    Estimation of Disease Rates in Small Areas: A New Mixed Model for Spatial Dependence.
#'    In M. Halloran and D. Berry (Eds.), Statistical Models in Epidemiology,
#'    the Environment and Clinical Trials, 135-178.
#'
#'  T. Park and G. Casella (2008).
#'    The Bayesian Lasso.
#'    Journal of the American Statistical Association 103(482), 681-686.
#'
#'  H. Rue and L. Held (2005).
#'    Gaussian Markov Random Fields.
#'    Chapman & Hall/CRC.
gen <- function(formula = ~ 1, factor=NULL,
                remove.redundant=FALSE, drop.empty.levels=FALSE, X=NULL,
                var=NULL, prior=NULL, Q0=NULL, PX=TRUE,
                GMRFmats=NULL, priorA=NULL, Leroux=FALSE,
                R0=NULL, RA=NULL, constr=NULL,
                S0=NULL, SA=NULL,
                formula.gl=NULL,
                name="", sparse=NULL, perm=NULL,
                debug=FALSE, e=parent.frame()) {

  type <- "gen"  # for generic (random effects)
  if (name == "") stop("missing model component name")

  # check supplied var component (and if NULL leave them to be filled in later)
  if (!is.null(var)) var <- match.arg(var, c("unstructured", "diagonal", "scalar", "fixed"))
  if (!is.null(factor) && !inherits(factor, "formula")) stop("element 'factor' of a model component must be a formula")
  if (is.list(Leroux)) {
    if (!all(names(Leroux) %in% c("a", "b", "a.star", "b.star"))) stop("invalid Leroux model component parameter list")
    if (is.null(Leroux$a)) Leroux$a <- 1
    if (is.null(Leroux$b)) Leroux$b <- 1
    if (is.null(Leroux$a.star)) Leroux$a.star <- 1
    if (is.null(Leroux$b.star)) Leroux$b.star <- 1
    if (any(unlist(Leroux) < 0)) stop("parameters of beta distribution cannot be negative")
    Leroux_type <- "general"
  } else {
    if (is.numeric(Leroux)) {
      if (length(Leroux) != 1L) stop("only a single number can be specified for a Leroux parameter")
      if (Leroux < 0 || Leroux > 1) stop("Leroux parameter must be between 0 and 1")
      Leroux_type <- "fixed"
    } else {
      if (!is.logical(Leroux) || (length(Leroux) != 1L)) stop("wrong input for 'Leroux'")
      Leroux_type <- if (Leroux) "default" else "none"
    }
  }
  Leroux_update <- Leroux_type %in% c("general", "default")

  if (is.null(X)) {
    X <- compute_X(formula, factor, remove.redundant=remove.redundant,
      drop.empty.levels=drop.empty.levels, sparse=sparse, data=e$data)
    factor.cols.removed <- attr(X, "factor.cols.removed")
    if (!is.null(factor.cols.removed)) attr(X, "factor.cols.removed") <- NULL
  } else {
    factor.cols.removed <- NULL
  }
  if (nrow(X) != e$n) stop("design matrix with incompatible number of rows")
  e$coef.names[[name]] <- colnames(X)
  X <- economizeMatrix(X, strip.names=TRUE)
  q <- ncol(X)
  in_block <- name %in% unlist(e$block)

  self <- environment()

  # fast GMRF prior currently only for IGMRF without user-defined constraints
  fastGMRFprior <- is.null(priorA) && allow_fastGMRFprior(self) && is.null(RA) && is.null(R0) && is.null(SA) && is.null(S0)
  # by default, IGMRF constraints are imposed, unless the GMRF is proper
  # NB constr currently only refers to QA, not Q0
  if (is.null(constr)) constr <- !is_proper_GMRF(self)
  if (is.null(GMRFmats)) {
    GMRFmats <- compute_GMRF_matrices(factor, e$data,
      D=fastGMRFprior || !is.null(priorA),
      Q=!(fastGMRFprior || !is.null(priorA)),
      R=constr, sparse=if (in_block) TRUE else NULL,
      cols2remove=factor.cols.removed, drop.zeros=TRUE
    )
  }
  if (fastGMRFprior || !is.null(priorA)) {  # compute lD x l matrix DA
    DA <- GMRFmats$D
    l <- ncol(DA)
    if (is.null(priorA))
      QA <- economizeMatrix(crossprod(DA), sparse=if (in_block) TRUE else NULL, symmetric=TRUE)
  } else {  # compute precision matrix QA only
    QA <- economizeMatrix(GMRFmats$Q, sparse=if (in_block) TRUE else NULL, symmetric=TRUE)
    l <- nrow(QA)
  }
  if (q %% l != 0L) stop("incompatible dimensions of design and precision matrices")
  q0 <- q %/% l  # --> q = q0 * l

  if (!is.null(priorA)) {  # scaled precision matrix QA = DA' QD DA with QD modeled
    lD <- nrow(DA)
    name_omega <- paste0(name, "_omega")
    if (Leroux_type != "none") stop("not implemented: scaled Leroux type correlation")
    switch(priorA$type,
      invchisq = {
        priorA <- pr_invchisq(priorA$df, priorA$scale, lD, !e$prior.only)
        df.data.omega <- q0  # TODO is this always the correct value?
      },
      exp = priorA <- pr_exp(priorA$scale, lD, !e$prior.only)
    )
  }

  if (q0 == 1L) {
    var <- "scalar"  # single varying effect
  } else {
    if (is.null(var))
      if (l == 1L) 
        var <- "diagonal"  # single group, default ARD prior
      else
        var <- "unstructured"
  }

  PX.defaults <- list(mu0=0, Q0=1, data.scale=TRUE, vector=var %in% c("diagonal", "unstructured"), sparse=NULL)
  if (is.list(PX)) {
    usePX <- TRUE
    if (!all(names(PX) %in% names(PX.defaults))) stop("invalid 'PX' options list")
    PX <- modifyList(PX.defaults, PX)
  } else {
    if (is.logical(PX) && length(PX) == 1L) {
      usePX <- PX
      PX <- PX.defaults
    } else
      stop("wrong input for 'PX'")
  }
  rm(PX.defaults)
  if (usePX) {
    if (is.null(PX$sparse)) PX$sparse <- (q0 > 10L)
    if (q0 == 1L) PX$vector <- FALSE
    if (!PX$vector) PX$sparse <- FALSE
    PX$dim <- if (PX$vector) q0 else 1L

    if (PX$vector) {
      if (length(PX$mu0) == 1L) PX$mu0 <- rep.int(PX$mu0, q0)
      if (is.numeric(PX$Q0)) {
        if (length(PX$Q0) == 1L) PX$Q0 <- rep.int(PX$Q0, q0)
        PX$Q0 <- Cdiag(PX$Q0)
      }
      PX$Q0 <- economizeMatrix(PX$Q0, sparse=PX$sparse, symmetric=TRUE)
      base_tcrossprod <- base::tcrossprod  # faster access (Matrix promotes tcrossprod to S4 generic)
    }
    if (length(PX$mu0) != PX$dim) stop("wrong length for 'PX$mu0'")
    PX$Qmu0 <- PX$Q0 %m*v% PX$mu0
    name_xi <- paste0(name, "_xi_")  # trailing '_' so that it is not stored by MCMCsim even if store.all=TRUE
  }

  if (Leroux_type != "none") {
    if (isUnitDiag(QA)) warning("Leroux extension for iid random effects has no effect")
    if (Leroux_type == "fixed") {
      # compute the weighted average once and for all
      QA <- economizeMatrix(Leroux * QA + (1 - Leroux) * CdiagU(l), symmetric=TRUE, sparse=if (in_block) TRUE else NULL)
    } else {
      # use a sparse (sum) template
      idL <- CdiagU(l)
      mat_sum_Leroux <- make_mat_sum(M1=QA, M2=idL)
      # base the determinant template on 0.5*(QA + I) to ensure it is non-singular
      det <- make_det(mat_sum_Leroux(QA, idL, 0.5, 0.5))
    }
  }

  # construct equality restriction matrix
  # TODO if R is supplied, it is assumed to be the full constraint matrix...
  #      in that case the derivation from R0 and RA can be skipped
  if (!is.null(R0)) {
    if (var != "scalar") stop("R0 constraint matrix only allowed if var='scalar'")
    if (nrow(R0) != q0) stop("incompatible matrix R0")
  }

  # if a matrix 'RA' is provided by the user, use these constraints in addition to possible GMRF constraints
  if (is.null(RA)) {
    RA <- GMRFmats$R
  } else {
    if (nrow(RA) != l) stop("incompatible matrix RA")
    if (!is.null(GMRFmats$R)) {
      RA <- economizeMatrix(cbind(RA, GMRFmats$R), allow.tabMatrix=FALSE)
      RA <- remove_redundancy(RA)  # combination of user-supplied and IGMRF constraints may contain redundant columns
    }
  }

  R <- switch(paste0(is.null(RA), is.null(R0)),
    TRUETRUE = NULL,
    TRUEFALSE = kronecker(CdiagU(l), R0),
    FALSETRUE = kronecker(RA, CdiagU(q0)),
    FALSEFALSE = cbind(kronecker(CdiagU(l), R0), kronecker(RA, CdiagU(q0)))
  )
  if (!is.null(R)) {
    # tabMatrix doesn't work yet because of lack of solve method
    R <- economizeMatrix(R, allow.tabMatrix=FALSE)
  }
  if (!is.null(R0) && !is.null(RA)) R <- remove_redundancy(R)
  rm(R0, GMRFmats)

  # construct inequality restriction matrix (TODO S instead of S0, SA)
  if (!is.null(S0)) {
    if (nrow(S0) != q0) stop("incompatible matrix S0")
  }
  if (!is.null(SA)) {
    if (nrow(SA) != l) stop("incompatible matrix SA")
  }
  S <- switch(paste0(is.null(SA), is.null(S0)),
    TRUETRUE = NULL,
    TRUEFALSE = kronecker(CdiagU(l), S0),
    FALSETRUE = kronecker(SA, CdiagU(q0)),
    FALSEFALSE = cbind(kronecker(CdiagU(l), S0), kronecker(SA, CdiagU(q0)))
  )
  rm(S0, SA)

  if (is.null(prior)) {
    prior <- switch(var,
      unstructured = pr_invwishart(df=q0+1, scale=diag(q0)),
      diagonal=, scalar = pr_invchisq(if (usePX) 1 else 1e-3, 1)
    )
  }
  if (!(prior$type %in% if (var == "unstructured") "invwishart" else c("invchisq", "exp"))) stop("unsupported prior")
  if (is.list(prior$df)) stop("not supported: modeled df for random effect variance prior")
  switch(prior$type,
    invwishart = prior <- pr_invwishart(prior$df, prior$scale, q0),
    invchisq = {
      if (var == "scalar") {
        prior <- pr_invchisq(prior$df, prior$scale, 1L, !e$prior.only)
        df.data <- if (is.null(R)) q else q - ncol(R)
      } else {  # var "diagonal"
        prior <- pr_invchisq(prior$df, prior$scale, q0, !e$prior.only)
        # df based on number of unconstrained coefficients
        df.data <- if (is.null(RA)) l else l - ncol(RA)
      }
    },
    exp = {
      if (var == "scalar") {
        prior <- pr_exp(prior$scale, 1L, !e$prior.only)
        ppar <- if (is.null(R)) q else q - ncol(R)
      } else {  # var "diagonal"
        prior <- pr_exp(prior$scale, q0, !e$prior.only)
        ppar <- if (is.null(RA)) l else l - ncol(RA)
      }
      ppar <- 1 - ppar/2
    }
  )

  if (var == "scalar") {
    # also check given precision matrix Q0, if any
    if (!is.null(Q0)) {
      Q0 <- economizeMatrix(Q0, symmetric=TRUE, vec.diag=TRUE)
      if (is.vector(Q0)) {
        if (!(length(Q0) %in% c(1L, q0))) stop("incompatible 'Q0'")
      } else {
        if (!identical(dim(Q0), c(q0, q0))) stop("incompatible 'Q0'")
      }
      if (is.vector(Q0) && all(Q0 == 1)) Q0 <- NULL  # identity Q0, can be ignored
    }
  } else {
    if (var == "unstructured") {
      # TODO add rprior and draw methods to pr_invwishart
      df <- prior$df + l
      if (!is.null(RA)) df <- df - ncol(RA)
    }
    if (!is.null(Q0)) warning("'Q0' ignored")
  }
  rm(RA)

  gl <- !is.null(formula.gl)  # group-level predictors?
  if (gl) {  # group-level covariates
    if (!inherits(formula.gl, "formula")) stop("'formula.gl' must be a formula")
    vars <- as.list(attr(terms(formula.gl), "variables"))[-1L]
    if (length(vars) != 1L || as.character(vars[[1L]][[1L]]) != "glreg") stop("'formula.gl' should be a formula with a single 'glreg' component")
    name_gl <- paste0(name, "_gl")
    glp <- vars[[1L]]
    glp$name <- name_gl
    glp$e <- self
    gl <- TRUE
    rm(vars)
  }

  if (e$modeled.Q) {
    if (gl) {
      # ensure that XX is always sparse because we need a sparse block diagonal template in this case
      X <- economizeMatrix(X, sparse=TRUE)  # --> crossprod_sym(X, Q) also sparse
    }
    XX <- crossprod_sym(X, crossprod_sym(Diagonal(x=runif(e$n, 0.9, 1.1)), e$Q0))
  } else {
    XX <- economizeMatrix(crossprod_sym(X, e$Q0), symmetric=TRUE)
  }

  # for both memory and speed efficiency define unit_Q case (random intercept and random slope with scalar var components, no constraints)
  unit_Q <- is.null(priorA) && isUnitDiag(QA) && (var == "scalar" && is.null(Q0)) && !constr

  name_sigma <- paste0(name, "_sigma")
  if (var == "unstructured") name_rho <- paste0(name, "_rho")

  if (q0 > 1L) {  # add labels for sigma, rho, xi to e$coef.names
    varnames <- try(colnames(model_Matrix(formula, e$data[1L, , drop=FALSE], remove.redundant=remove.redundant, sparse=sparse)), silent=TRUE)
    if (inherits(varnames, "try-error")) {
      # something wrong, which can happen for formulas like ~as.factor(area), as well as for character variables --> use full dataset
      # TODO function that uses full dataset but only returns colnames
      varnames <- colnames(model_Matrix(formula, e$data, remove.redundant=remove.redundant, sparse=sparse))
    }
    if (var %in% c("unstructured", "diagonal") || (usePX && PX$vector)) {
      e$coef.names[[name_sigma]] <- varnames
      if (var == "unstructured")
        e$coef.names[[name_rho]] <- outer(varnames, varnames, FUN=paste, sep=":")[upper.tri(diag(q0))]
      if (usePX && PX$vector)
        e$coef.names[[name_xi]] <- varnames
    }
    if (gl) {
      if (is.null(e$coef.names[[name_gl]]))
        e$coef.names[[name_gl]] <- varnames
      else
        e$coef.names[[name_gl]] <- as.vector(outer(varnames, e$coef.names[[name_gl]], FUN=paste, sep=":"))
    }
  }

  linpred <- function(p) X %m*v% p[[name]]

  if (gl) {  # group-level covariates
    glp <- eval(glp)
    sparse_template(self, modeled.Q=e$modeled.Q, prior.only=e$prior.only)
    if (!in_block && !e$prior.only) glp$R <- NULL
  } else {
    sparse_template(self, modeled.Q=e$modeled.Q, prior.only=e$prior.only)
  }

  # BEGIN rprior function
  # draw from prior; draws are independent and stored in p
  rprior <- function(p) {}
  if (usePX) {
    if (PX$data.scale && !e$sigma.fixed)
      rprior <- add(rprior, quote(xi <- PX$mu0 + drawMVN_Q(PX$Q0, sd=p[["sigma_"]])))
    else
      rprior <- add(rprior, quote(xi <- PX$mu0 + drawMVN_Q(PX$Q0)))
    rprior <- add(rprior, bquote(p[[.(name_xi)]] <- xi))
    rprior <- add(rprior, quote(inv_xi <- 1/xi))
  } else {
    rprior <- add(rprior, quote(inv_xi <- 1))
  }

  if (var == "unstructured") {
    if (is.list(prior$scale)) {
      psi0 <- prior$scale$df / prior$scale$scale
      # TODO C++ version of dense diagonal matrix creation
      if (prior$scale$common)
        rprior <- add(rprior, quote(psi0 <- diag(rep.int(rchisq_scaled(1L, prior$scale$df, psi=psi0), q0))))
      else
        rprior <- add(rprior, quote(psi0 <- diag(rchisq_scaled(q0, prior$scale$df, psi=psi0))))
    } else {
      psi0 <- prior$scale  # matrix
    }
    rprior <- add(rprior, quote(Qraw <- rWishart(1L, df=prior$df, Sigma=inverseSPD(psi0))[,,1L]))
    if (usePX && PX$vector)
      rprior <- add(rprior, quote(inv_xi_qform <- base_tcrossprod(inv_xi)))
    else
      rprior <- add(rprior, quote(inv_xi_qform <- inv_xi^2))
    rprior <- add(rprior, quote(Qv <- Qraw * inv_xi_qform))  # use Qv to distinguish from data-level precision Q
    # convert precision matrix to standard errors and correlations
    names_se_cor <- c(name_sigma, name_rho)
    rprior <- add(rprior, quote(p[names_se_cor] <- prec2se_cor(Qv)))
  } else {
    rprior <- add(rprior, quote(Qraw <- 1 / prior$rprior()))
    rprior <- add(rprior, quote(Qv <- Qraw * inv_xi^2))
    rprior <- add(rprior, bquote(p[[.(name_sigma)]] <- sqrt(1/Qv)))
  }
  if (Leroux_update) {
    name_Leroux <- paste0(name, "_Leroux")
    if (Leroux_type == "default")
      rprior <- add(rprior, bquote(p[[.(name_Leroux)]] <- runif(1L)))
    else
      rprior <- add(rprior, bquote(p[[.(name_Leroux)]] <- rbeta(1L, Leroux$a, Leroux$b)))
    rprior <- add(rprior, bquote(QA <- mat_sum_Leroux(QA, idL, p[[.(name_Leroux)]], 1 - p[[.(name_Leroux)]])))
  }

  if (gl) {  # draw group-level predictor coefficients
    rprior <- add(rprior, bquote(p[[.(name_gl)]] <- drawMVN_Q(glp$Q0, sd=.(if (e$sigma.fixed) 1 else quote(p[["sigma_"]])))))
  }

  if (!is.null(priorA)) {
    switch(priorA$type,
      invchisq = {
        if (is.list(priorA$df)) {
          name_df <- paste0(name, "_df")
          rprior <- add(rprior, bquote(p[[.(name_df)]] <- priorA$rprior_df()))
          rprior <- add(rprior, bquote(p[[.(name_omega)]] <- priorA$rprior(p[[.(name_df)]])))
        } else {
          rprior <- add(rprior, bquote(p[[.(name_omega)]] <- priorA$rprior()))
        }
      },
      exp = {
        rprior <- add(rprior, bquote(p[[.(name_omega)]] <- priorA$rprior()))
      }
    )
    rprior <- add(rprior, bquote(QA <- crossprod_sym(DA, 1/p[[.(name_omega)]])))
  }

  # draw coefficient from prior
  if (unit_Q) {  # slightly faster alternative for random intercept models (or random slope with scalar variance)
    rprior <- add(rprior, bquote(p[[.(name)]] <- rnorm(.(q), sd=p[[.(name_sigma)]])))
  } else {
    if (fastGMRFprior) {
      rGMRFprior <- NULL  # to be created at the first occasion rprior is called
      rprior <- add(rprior, quote(if (is.null(rGMRFprior)) setup_priorGMRFsampler(self, Qv)))
      rprior <- add(rprior, bquote(p[[.(name)]] <- rGMRFprior(Qv)))
    } else {
      if (q0 == 1L)
        rprior <- add(rprior, quote(Q <- Qv * QA))
      else
        rprior <- add(rprior, quote(Q <- kron_prod(QA, Qv)))
      rMVNprior <- NULL  # to be created at the first occasion rprior is called
      rprior <- add(rprior, quote(if (is.null(rMVNprior)) setup_priorMVNsampler(self, Q)))
      rprior <- add(rprior, bquote(p[[.(name)]] <- rMVNprior(p, Q)))
    }
  }
  rprior <- add(rprior, quote(p))
  # END rprior function

  if (e$prior.only) return(self)  # no need for posterior draw function in this case

  if (!is.null(priorA) || is.list(prior$scale) || Leroux_update) {
    # store Qraw --> no need to recompute at next MCMC iteration
    name_Qraw <- paste0(name, "_Qraw_")
  }

  # BEGIN draw function
  draw <- function(p) {}
  if (debug) draw <- add(draw, quote(browser()))
  if (e$single.block && length(e$mod) == 1L) {
    # optimalization in case of a single regression term, only in case of single mc (due to PX)
    draw <- add(draw, quote(p$e_ <- e$y_eff()))
  } else {
    if (e$e.is.res)
      draw <- add(draw, bquote(p$e_ <- p[["e_"]] + X %m*v% p[[.(name)]]))
    else
      draw <- add(draw, bquote(p$e_ <- p[["e_"]] - X %m*v% p[[.(name)]]))
  }

  if (!e$e.is.res && e$single.block && !e$use.offset && length(e$mod) != 1L) {
    # need to correct here Q_e function; this should not be necessary if we draw all xi's in a single block too !!
    if (e$family$link == "probit")
      draw <- add(draw, quote(Xy <- crossprod_mv(X, e$Q_e(p) - p[["e_"]])))
    else
      draw <- add(draw, quote(Xy <- crossprod_mv(X, e$Q_e(p) - p[["Q_"]] * p[["e_"]])))
  } else {
    draw <- add(draw, quote(Xy <- crossprod_mv(X, e$Q_e(p))))
  }
  if (e$modeled.Q) {
    if (e$Q0.type == "symm")
      draw <- add(draw, quote(XX <- crossprod_sym(X, p[["QM_"]])))
    else
      draw <- add(draw, quote(XX <- crossprod_sym(X, p[["Q_"]])))
  }

  if (usePX)
    draw <- add(draw, bquote(inv_xi <- 1 / p[[.(name_xi)]]))
  else
    draw <- add(draw, quote(inv_xi <- 1))
  draw <- add(draw, bquote(coef_raw <- inv_xi * p[[.(name)]]))
  if (q0 == 1L)
    draw <- add(draw, quote(M_coef_raw <- coef_raw))
  else
    draw <- add(draw, bquote(M_coef_raw <- matrix(coef_raw, nrow=.(l), byrow=TRUE)))

  draw <- add(draw, bquote(tau <- .(if (e$sigma.fixed) 1 else quote(1 / p[["sigma_"]]^2))))

  if (usePX) {
    if (PX$vector) {
      if (PX$sparse) {  # sparse M_ind, Dv
        M_ind <- as(as(kronecker(rep.int(1, l), CdiagU(q0)), "CsparseMatrix"), "dgCMatrix")
        mat_sum_xi <- make_mat_sum(M0=PX$Q0, M1=crossprod_sym(M_ind, XX))
        chol_xi <- build_chol(mat_sum_xi(crossprod_sym(M_ind, XX)))
        draw <- add(draw, quote(Dv <- M_ind))
        draw <- add(draw, quote(attr(Dv, "x") <- as.numeric(M_coef_raw)))
        if (PX$data.scale)
          draw <- add(draw, quote(chol_xi$update(mat_sum_xi(crossprod_sym(Dv, XX)))))
        else
          draw <- add(draw, quote(chol_xi$update(mat_sum_xi(crossprod_sym(Dv, XX), w1=tau))))
      } else {  # matrix M_ind, Dv
        M_ind <- kronecker(rep.int(1, l), diag(q0))
        chol_xi <- build_chol(crossprod_sym(M_ind, XX) + PX$Q0)
        draw <- add(draw, quote(Dv <- coef_raw * M_ind))
        if (PX$data.scale)
          draw <- add(draw, quote(chol_xi$update(crossprod_sym(Dv, XX) + PX$Q0)))
        else
          draw <- add(draw, quote(chol_xi$update(crossprod_sym(Dv, XX) * tau + PX$Q0)))
      }
      if (PX$data.scale)
        draw <- add(draw, bquote(xi <- drawMVN_cholQ(chol_xi, crossprod_mv(Dv, Xy) + PX$Qmu0, sd=.(if (e$sigma.fixed) 1 else quote(p[["sigma_"]])))))
      else
        draw <- add(draw, quote(xi <- drawMVN_cholQ(chol_xi, crossprod_mv(Dv, Xy) * tau + PX$Qmu0)))
    } else {  # scalar xi
      if (PX$data.scale) {
        draw <- add(draw, quote(V <- 1 / (dotprodC(coef_raw, XX %m*v% coef_raw) + PX$Q0)))
        draw <- add(draw, quote(E <- V * (dotprodC(coef_raw, Xy) + PX$Qmu0)))
        draw <- add(draw, bquote(xi <- rnorm(1L, mean=E, sd=.(if (e$sigma.fixed) quote(sqrt(V)) else quote(p[["sigma_"]] * sqrt(V))))))
      } else {
        draw <- add(draw, quote(V <- 1 / (dotprodC(coef_raw, XX %m*v% coef_raw) * tau + PX$Q0)))
        draw <- add(draw, quote(E <- V * (dotprodC(coef_raw, Xy) * tau + PX$Qmu0)))
        draw <- add(draw, quote(xi <- rnorm(1L, mean=E, sd=sqrt(V))))
      }
    }
    if (!is.null(S)) {
      stop("TBI: inequality constrained coefficients with parameter expansion")  # --> f.c. of xi also TMVN(?)
      draw <- add(draw, bquote(p[[.(name)]] <- xi * coef_raw))  # need updated coefficient as input for TMVN sampler
    }
  } else {
    xi <- 1  # no parameter expansion
  }

  if (e$modeled.Q) rm(XX)  # XX recomputed at each iteration

  if (gl) {  # group-level predictors
    if (usePX) {
      # MH step
      if (PX$vector) {
        draw <- add(draw, bquote(logr <- .(glp$p0) * (sum(log(abs(xi))) - sum(log(abs(p[[.(name_xi)]]))))))
      } else {
        draw <- add(draw, bquote(logr <- .(q0 * glp$p0) * (sum(log(abs(xi))) - sum(log(abs(p[[.(name_xi)]]))))))
      }
      draw <- add(draw, bquote(delta <- (xi * inv_xi) * p[[.(name_gl)]]))
      draw <- add(draw, quote(logr <- logr - 0.5 * dotprodC(delta, glp$Q0 %m*v% delta)))
      draw <- add(draw, bquote(delta <- p[[.(name_gl)]]))
      draw <- add(draw, quote(logr <- logr + 0.5 * dotprodC(delta, glp$Q0 %m*v% delta)))
      draw <- add(draw, bquote(if (logr < log(runif(1L))) xi <- p[[.(name_xi)]]))  # reject, otherwise accept
    }
    # NB use old value of xi to get 'raw' group-level coefficients
    if (q0 == 1L) {
      draw <- add(draw, bquote(M_coef_raw <- M_coef_raw - glp$X %m*v% (inv_xi * p[[.(name_gl)]])))
    } else {
      draw <- add(draw, bquote(M_coef_raw <- M_coef_raw - glp$X %*% matrix(inv_xi * p[[.(name_gl)]], nrow=.(glp$p0), byrow=TRUE)))
    }
  }
  if (usePX) {
    draw <- add(draw, bquote(p[[.(name_xi)]] <- xi))
    draw <- add(draw, quote(inv_xi <- 1 / xi))
  }

  if (!is.null(priorA)) {
    if (q0 == 1L)  # DAM is lD vector
      draw <- add(draw, quote(DAM <- DA %m*v% M_coef_raw))
    else  # DAM is of type matrix (lD x q0) as M_coef_raw is
      draw <- add(draw, quote(DAM <- DA %*% M_coef_raw))
    switch(var,
      unstructured = {
        draw <- add(draw, bquote(SSR <- .rowSums((DAM %*% p[[.(name_Qraw)]]) * DAM, .(lD), .(q0))))
      },
      diagonal = {
        #draw <- add(draw, bquote(SSR <- .rowSums((DAM %*% Cdiag(p[[.(name_temp)]]$Qraw)) * DAM, .(lD), .(q0))))
        draw <- add(draw, bquote(SSR <- .rowSums(rep_each(p[[.(name_Qraw)]], .(lD)) * DAM^2, .(lD), .(q0))))
      },
      scalar = {
        if (q0 == 1L) {
          if (is.null(Q0)) {
            draw <- add(draw, bquote(SSR <- p[[.(name_Qraw)]] * DAM^2))
          } else {
            draw <- add(draw, bquote(SSR <- Q0 * p[[.(name_Qraw)]] * DAM^2))
          }
        } else {
          if (is.null(Q0)) {
            draw <- add(draw, bquote(SSR <- p[[.(name_Qraw)]] * .rowSums(DAM * DAM, .(lD), .(q0))))
          } else {
            draw <- add(draw, bquote(SSR <- p[[.(name_Qraw)]] * .rowSums((DAM %m*m% Q0) * DAM, .(lD), .(q0))))
          }
        }
      }
    )
    switch(priorA$type,
      invchisq = {
        if (is.list(priorA$df)) {
          draw <- add(draw, bquote(p[[.(name_df)]] <- priorA$draw_df(p[[.(name_df)]], 1 / p[[.(name_omega)]])))
          if (is.list(priorA$scale)) {
            draw <- add(draw, bquote(p[[.(name_omega)]] <- priorA$draw(p[[.(name_df)]], df.data.omega, SSR, 1 / p[[.(name_omega)]])))
          } else {
            draw <- add(draw, bquote(p[[.(name_omega)]] <- priorA$draw(p[[.(name_df)]], df.data.omega, SSR)))
          }
        } else {
          if (is.list(priorA$scale)) {
            draw <- add(draw, bquote(p[[.(name_omega)]] <- priorA$draw(df.data.omega, SSR, 1 / p[[.(name_omega)]])))
          } else {
            draw <- add(draw, bquote(p[[.(name_omega)]] <- priorA$draw(df.data.omega, SSR)))
          }
        }
      },
      exp = {
        aparA <- 2 / priorA$scale
        draw <- add(draw, bquote(p[[.(name_omega)]] <- priorA$draw(1 - q0/2, aparA, SSR)))
      }
    )
    # then update QA
    draw <- add(draw, bquote(QA <- crossprod_sym(DA, 1 / p[[.(name_omega)]])))
  }

  if (Leroux_update) {
    # draw Leroux parameter
    draw <- add(draw, quote(p <- draw_Leroux(p, self, M_coef_raw)))
    # compute QA using current value of Leroux parameter
    # alternatively, store it in p[[name_temp]]$QA_L and recompute only if a new Leroux parameter has been accepted
    draw <- add(draw, bquote(QA <- mat_sum_Leroux(QA, idL, p[[.(name_Leroux)]], 1 - p[[.(name_Leroux)]])))
  }

  # draw V
  if (var == "unstructured") {
    if (is.list(prior$scale)) {
      # Qraw stored in p; NB Qv has changed because of updated xi (PX), but Qraw has not
      if (prior$scale$common) {
        draw <- add(draw, bquote(psi0 <- psi0 + sum(diagC(p[[.(name_Qraw)]]))))
        draw <- add(draw, quote(psi0 <- rchisq_scaled(1L, q0 * prior$df + prior$scale$df, psi = psi0)))
      } else {
        draw <- add(draw, bquote(psi0 <- psi0 + diagC(p[[.(name_Qraw)]])))
        draw <- add(draw, bquote(psi0 <- rchisq_scaled(.(q0), prior$df + prior$scale$df, psi=psi0)))
      }
      draw <- add(draw, quote(Qraw <- rWishart(1L, df=df, Sigma = inverseSPD(add_diagC(crossprod_sym(M_coef_raw, QA), psi0)))[,,1L]))
    } else {
      # TODO draw V_raw using invWishart and form Qv by solving --> 1 instead of 2 solves
      draw <- add(draw, quote(Qraw <- rWishart(1L, df=df, Sigma = inverseSPD(psi0 + crossprod_sym(M_coef_raw, QA)))[,,1L]))
    }
    if (usePX && PX$vector)
      draw <- add(draw, quote(inv_xi_qform <- base_tcrossprod(inv_xi)))
    else
      draw <- add(draw, quote(inv_xi_qform <- inv_xi^2))
    draw <- add(draw, quote(Qv <- Qraw * inv_xi_qform))
    # convert precision matrix to standard errors and correlations
    draw <- add(draw, quote(p[names_se_cor] <- prec2se_cor(Qv)))
  } else {
    if (var == "diagonal") {
      draw <- add(draw, bquote(SSR <- .colSums(M_coef_raw * (QA %*% M_coef_raw), .(l), .(q0))))
    } else {  # scalar variance
      if (q0 == 1L) {
        if (is.null(Q0))
          draw <- add(draw, quote(SSR <- dotprodC(M_coef_raw, QA %m*v% M_coef_raw)))
        else
          draw <- add(draw, quote(SSR <- Q0 * dotprodC(M_coef_raw, QA %m*v% M_coef_raw)))
      } else {
        if (is.null(Q0))
          draw <- add(draw, quote(SSR <- sum(M_coef_raw * (QA %*% M_coef_raw))))
        else
          draw <- add(draw, quote(SSR <- sum(M_coef_raw * (QA %*% M_coef_raw %m*m% Q0))))
      }
    }
    switch(prior$type,
      invchisq = {
        if (is.list(prior$scale))
          draw <- add(draw, bquote(Qraw <- 1 / prior$draw(df.data, SSR, p[[.(name_Qraw)]])))
        else
          draw <- add(draw, quote(Qraw <- 1 / prior$draw(df.data, SSR)))
      },
      exp = {
        apar <- 2 / prior$scale
        draw <- add(draw, quote(Qraw <- 1 / prior$draw(ppar, apar, SSR)))
      }
    )
    # NB if xi is a vector so is Qv, even for scalar Qraw
    draw <- add(draw, quote(Qv <- Qraw * inv_xi^2))
    if (is.null(Q0)) {
      draw <- add(draw, bquote(p[[.(name_sigma)]] <- sqrt(1/Qv)))
    } else {
      draw <- add(draw, quote(sqrtQv <- sqrt(Qv)))
      draw <- add(draw, quote(Qv <- scale_mat(Q0, sqrtQv)))
      draw <- add(draw, bquote(p[[.(name_sigma)]] <- 1/sqrtQv))
    }
  }
  if (!is.null(priorA) || is.list(prior$scale) || Leroux_update) {
    # store Qraw for next iteration -> saves reconstructing it from sigma, rho, xi
    draw <- add(draw, bquote(p[[.(name_Qraw)]] <- Qraw))
  }

  # draw coefficients
  if (gl) {
    i.v <- seq_len(q)  # indices for random effect vector in u=(v, alpha)
    i.alpha <- (q + 1L):(q + glp$p0 * q0)  # indices for group-level effect vector in u=(v, alpha)
  }
  if (in_block) {
    name_Q <- paste0(name, "_Q_")  # trailing "_" --> only temporary storage
    if (gl) {
      draw <- add(draw, bquote(p[[.(name_Q)]] <- kron_prod(glp$QA.ext, Qv, values.only=TRUE)))
    } else {
      draw <- add(draw, bquote(p[[.(name_Q)]] <- kron_prod(QA, Qv, values.only=TRUE)))
    }
    draw <- add(draw, bquote(p[[.(name)]] <- xi * coef_raw))
  } else {
    if (gl) {  # block sampling of (coef, glp)
      if (e$modeled.Q)
        draw <- add(draw, quote(attr(glp$XX.ext, "x") <- c(XX@x, glp$Q0@x)))
      if (Leroux_update)
        draw <- add(draw, quote(glp$QA.ext <- crossprod_sym(glp$IU0, QA)))
      draw <- add(draw, quote(Qlist <- update(glp$XX.ext, glp$QA.ext, Qv, 1/tau)))
      draw <- add(draw, bquote(coef <- MVNsampler$draw(p, .(if (e$sigma.fixed) 1 else quote(p[["sigma_"]])), Q=Qlist$Q, Imult=Qlist$Imult, Xy=c(Xy, glp$Q0b0))[[.(name)]]))
      draw <- add(draw, bquote(p[[.(name)]] <- coef[i.v]))
      draw <- add(draw, bquote(p[[.(name_gl)]] <- coef[i.alpha]))
    } else {  # no blocking, no group-level component
      draw <- add(draw, quote(Qlist <- update(XX, QA, Qv, 1/tau)))
      draw <- add(draw, bquote(p[[.(name)]] <- MVNsampler$draw(p, .(if (e$sigma.fixed) 1 else quote(p[["sigma_"]])), Q=Qlist$Q, Imult=Qlist$Imult, Xy=Xy)[[.(name)]]))
    }
  }

  if (e$e.is.res)
    draw <- add(draw, bquote(p$e_ <- p[["e_"]] - X %m*v% p[[.(name)]]))
  else
    draw <- add(draw, bquote(p$e_ <- p[["e_"]] + X %m*v% p[[.(name)]]))
  draw <- add(draw, quote(p))
  # END draw function


  start <- function(p) {}

  # TODO check formats of user-provided start values
  if (!in_block) {
    if (gl) {
      # TODO conditional sampling if one of p[[.(name)]] and p[[.(name_gl)]] is provided
      start <- add(start, bquote(coef <- MVNsampler$start(p, .(e$scale_y))[[.(name)]]))
      #if (!is.null(glp$R))
      #  start <- add(start, quote(coef <- constrain_cholQ(coef, ops$ch, glp$R)))
      start <- add(start, bquote(if (is.null(p[[.(name)]])) p[[.(name)]] <- coef[i.v]))
      start <- add(start, bquote(if (is.null(p[[.(name_gl)]])) p[[.(name_gl)]] <- coef[i.alpha]))
    } else {
      start <- add(start, bquote(if (is.null(p[[.(name)]])) p[[.(name)]] <- MVNsampler$start(p, .(e$scale_y))[[.(name)]]))
    }
  }

  if (usePX) {
    if (PX$vector)
      start <- add(start, bquote(if (is.null(p[[.(name_xi)]])) p[[.(name_xi)]] <- rep.int(1, .(q0))))
    else
      start <- add(start, bquote(if (is.null(p[[.(name_xi)]])) p[[.(name_xi)]] <- 1))
  }

  if (!is.null(priorA) || is.list(prior$scale) || Leroux_update) {
    switch(var,
      unstructured = {
        start <- add(start, bquote(if (is.null(p[[.(name_Qraw)]])) p[[.(name_Qraw)]] <- rWishart(1L, prior$df, if (is.list(prior$scale)) (1/psi0)*diag(q0) else inverseSPD(psi0))[,,1L]))
      },
      diagonal=, scalar = {
        if (prior$type == "exp")
          start <- add(start, bquote(if (is.null(p[[.(name_Qraw)]])) p[[.(name_Qraw)]] <- rexp(.(if (var == "diagonal") q0 else 1L), rate=1/prior$scale)))
        else
          start <- add(start, bquote(if (is.null(p[[.(name_Qraw)]])) p[[.(name_Qraw)]] <- rchisq_scaled(.(if (var == "diagonal") q0 else 1L), prior$df, psi=prior$psi0)))
      }
    )
    if (Leroux_update) {
      start <- add(start, bquote(if (is.null(p[[.(name_Leroux)]])) p[[.(name_Leroux)]] <- runif(1L)))
      name_detQA <- paste0(name, "_detQA_")
      start <- add(start, bquote(if (is.null(p[[.(name_detQA)]])) p[[.(name_detQA)]] <- det(2*p[[.(name_Leroux)]], 1 - 2*p[[.(name_Leroux)]])))
    }
    if (!is.null(priorA)) {
      if (is.list(priorA$df))
        start <- add(start, bquote(if (is.null(p[[.(name_df)]])) p[[.(name_df)]] <- runif(1L, 1, 25)))
      if (is.list(priorA$df) || is.list(priorA$scale))
        start <- add(start, bquote(if (is.null(p[[.(name_omega)]])) p[[.(name_omega)]] <- runif(.(lD), 0.75, 1.25)))
    }
  }

  start <- add(start, quote(p))

  # TODO rm more components no longer needed (Leroux_type, Leroux, other switches such as gl, ...)
  self
}


setup_priorGMRFsampler <- function(mc, Qv) {
  # TODO
  # - support local scale parameters DA --> diag(omega_i)^-1/2 DA
  # - in some cases perm=TRUE may be more efficient
  mc$cholDD <- build_chol(tcrossprod(mc$DA), perm=FALSE)
  mc$cholQv <- build_chol(Qv, perm=FALSE)
  rGMRF <- function(Qv) {}
  rGMRF <- add(rGMRF, quote(cholQv$update(Qv)))
  if (mc$q0 == 1L) {
    rGMRF <- add(rGMRF, bquote(Z <- Crnorm(.(nrow(mc$DA)))))
    rGMRF <- add(rGMRF, quote(Z <- cholQv$solve(Z, system="Lt", systemP=TRUE)))
    rGMRF <- add(rGMRF, quote(crossprod_mv(DA, cholDD$solve(Z))))
  } else {
    rGMRF <- add(rGMRF, bquote(Z <- matrix(Crnorm(.(mc$q0 * nrow(mc$DA))), nrow=.(mc$q0))))
    rGMRF <- add(rGMRF, quote(Z <- cholQv$solve(Z, system="Lt", systemP=TRUE)))
    rGMRF <- add(rGMRF, quote(coef <- crossprod(DA, cholDD$solve(t.default(Z)))))
    rGMRF <- add(rGMRF, quote(as.numeric(t.default(coef))))
  }
  mc$rGMRFprior <- rGMRF
  environment(mc$rGMRFprior) <- mc
}

setup_priorMVNsampler <- function(mc, Q) {
  rMVN <- function(p, Q) {}
  if (is_proper_GMRF(mc)) {
    mc$priorMVNsampler <- create_TMVN_sampler(Q, update.Q=TRUE, update.mu=FALSE, name="coef", R=mc$R)
  } else {
    if (is.null(mc$R)) stop("cannot sample from improper GMRF without constraints")
    # for IGMRF add tcrossprod(mc$R) to Q to make it non-singular (--> inefficient)
    # TODO maybe mc$R is removed and stored in (posterior) MVNsampler --> keep a reference(!) to R in mc
    mc$mat_sum_prior <- make_mat_sum(M0=economizeMatrix(tcrossprod(mc$R), symmetric=TRUE), M1=Q)
    mc$priorMVNsampler <- create_TMVN_sampler(mc$mat_sum_prior(Q), update.Q=TRUE, update.mu=FALSE, name="coef", R=mc$R)
    rMVN <- add(rMVN, quote(Q <- mat_sum_prior(Q)))
  }
  rMVN <- add(rMVN, quote(priorMVNsampler$draw(p, Q=Q)[["coef"]]))
  mc$rMVNprior <- rMVN
  environment(mc$rMVNprior) <- mc
}
