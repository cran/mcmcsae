#' Create a model object for group-level regression effects within a
#' generic random effects component.
#'
#' This function is intended to be used to specify the \code{formula.gl} argument to
#' the \code{\link{gen}} model component specification function.
#' Group-level predictors and hierarchical centering are
#' not used by default, and they currently cannot be used in a model component that is sampled
#' together with another model component in the same Gibbs block.
#'
#' @export
#' @param formula a formula specifying the group-level predictors to be used within a model
#'  component. If no \code{data} is supplied the group-level predictors are derived as
#'  group-level means from the unit-level data passed as \code{data} argument to
#'  \code{\link{create_sampler}} or \code{\link{generate_data}}.
#' @param remove.redundant whether redundant columns should be removed from the design matrix.
#'  Default is \code{FALSE}.
#' @param Q0 prior precision matrix for the group-level effects. The default is a
#'  zero matrix corresponding to a noninformative improper prior.
#  NB here we do not support non-zero prior mean b0
#' @param data group-level data frame in which the group-level variables specified in
#'  \code{formula} are looked up.
#' @param name the name of the model component. This name is used in the output of the
#'  MCMC simulation function \code{\link{MCMCsim}}. By default this name will be
#'  the name of the corresponding generic random effects component appended by '_gl'.
#' @param e for internal use only.
#' @return An object with precomputed quantities for sampling from
#'  prior or conditional posterior distributions for this model component. Only intended
#'  for internal use by other package functions.
glreg <- function(formula=NULL, remove.redundant=FALSE, Q0=NULL,
                  data=NULL, name="", e=parent.frame()) {

  type <- "reg"
  if (name == "") stop("missing model component name")

  if (is.null(data)) {
    # derive group-level design matrix as group means of unit-level data
    if (is.null(formula)) stop("no group-level formula: cannot derive group-level design matrix")
    if (intercept_only(formula)) {
      X <- matrix(rep.int(1, e$l), ncol=1L, dimnames=list(NULL, "(Intercept)"))
    } else {
      #X <- compute_X(formula, data=e$e$data, remove.redundant=remove.redundant, sparse=FALSE)
      X <- model_matrix(formula, e$e$data, sparse=FALSE)
      if (remove.redundant) X <- remove_redundancy(X)
      if (is.null(e$factor)) stop("cannot derive group-level design matrix")
      factor.info <- get_factor_info(e$factor, e$e$data)
      if ("spline" %in% factor.info$types) stop("unsupported combination: splines and group-level covariates")
      fac <- combine_factors(factor.info$variables, e$e$data)
      X <- economizeMatrix(crossprod(aggrMatrix(fac, mean=TRUE), X), strip.names=FALSE)
      rm(fac)
    }
  } else {  # formula + group-level data provided
    #X <- compute_X(formula, data=data, remove.redundant=remove.redundant, sparse=FALSE)
    X <- model_matrix(formula, data, sparse=FALSE)
    if (remove.redundant) X <- remove_redundancy(X)
    # TODO if formula is used, match levels of factor to glp$data
  }
  if (nrow(X) != e$l) stop("wrong number of rows of group-level design matrix")
  if (!is.null(colnames(X))) {
    e$e$coef.names[[name]] <- colnames(X)
  }
  X <- unname(X)  # l x p0

  p0 <- ncol(X)
  # if p0 > 1, each varying effect has its own coefficient and the same design matrix glp$X
  q <- p0 * e$q0  # number of gl coefficients

  # if modeled.Q need sparse Q0 and XX in order to combine them easily using a block-diagonal sparse template
  Q0 <- set_prior_precision(Q0, q, sparse=e$e$modeled.Q)
  informative.prior <- !is_zero_matrix(Q0)
  Q0b0 <- rep.int(0, q)  # need this in draw function, even though only b0=0 is supported

  # for modeled Q, the XX block of XX.ext is updated -> choose sparse
  XX.ext <- economizeMatrix(bdiag(e$XX, Q0), symmetric=TRUE, sparse=if (e$e$modeled.Q) TRUE else NULL)
  IU0 <- economizeMatrix(cbind(Diagonal(e$l), -X), drop.zeros=TRUE)
  if (e$Leroux_update) {
    # TODO crossprod_sym may introduce 0 fill-in --> kronecker template may fail for "unstructured" or "diagonal" var
    #      may drop (prune) 0s, but need to do that in each next crossprod_sym call too! (unavoidable?)
    w <- runif(1L, 0.25, 0.75)
    QA.ext <- crossprod_sym(IU0, e$mat_sum_Leroux(e$QA, e$idL, w, 1-w))
    rm(w)
  } else {
    QA.ext <- economizeMatrix(crossprod_sym(IU0, e$QA), symmetric=TRUE, drop.zeros=TRUE)
  }
  if (!is.null(e$R)) {
    R <- economizeMatrix(crossprod(kronecker(IU0, Diagonal(e$q0)), e$R), allow.tabMatrix=FALSE)
  }

  rm(data, e)
  environment()
}
