#' Create a model component object for a BART (Bayesian Additive Regression Trees)
#' component in the linear predictor
#'
#' This function is intended to be used on the right hand side of the
#' \code{formula} argument to \code{\link{create_sampler}} or
#' \code{\link{generate_data}}. It creates a BART term in the
#' model's linear predictor. To use this model component one needs
#' to have R package \pkg{dbarts} installed.
#'
#' @export
#' @param formula a formula specifying the predictors to be used in the BART
#'  model component. Variable names are looked up in the data frame
#'  passed as \code{data} argument to \code{\link{create_sampler}} or
#'  \code{\link{generate_data}}, or in \code{environment(formula)}.
#' @param X a design matrix can be specified directly, as an alternative
#'  to the creation of one based on \code{formula}. If \code{X} is
#'  specified \code{formula} is ignored.
#' @param name the name of the model component. This name is used in the output of the
#'  MCMC simulation function \code{\link{MCMCsim}}. By default the name will be 'bart'
#'  with the number of the model term attached.
#' @param debug if \code{TRUE} a breakpoint is set at the beginning of the posterior
#'  draw function associated with this model component. Mainly intended for developers.
#' @param e for internal use only.
#' @param ... parameters passed to \code{\link[dbarts]{dbarts}}.
#' @return an object with precomputed quantities and functions for sampling from
#'  prior or conditional posterior distributions for this model component. Intended
#'  for internal use by other package functions.
#' @references
#'  H.A. Chipman, E.I. Georgea and R.E. McCulloch (2010).
#'    BART: Bayesian additive regression trees.
#'    The Annals of Applied Statistics 4(1), 266-298.
bart <- function(formula, X=NULL, name="",
                 debug=FALSE, e=parent.frame(), ...) {

  type <- "bart"
  if (name == "") stop("missing model component name")

  if (!requireNamespace("dbarts", quietly=TRUE)) stop("Package dbarts required for a model including Bayesian Additive Regression Trees")

  if (e$modeled.Q && e$Q0.type == "symm") stop("BART component not compatible with non-diagonal residual variance matrix")

  if (is.null(X)) {
    # TODO maybe use dbarts::makeModelMatrixFromDataFrame instead
    X <- model_matrix(formula, e$data, sparse=FALSE)
  }
  X <- economizeMatrix(X, sparse=FALSE, strip.names=FALSE)

  control <- dbarts::dbartsControl(n.chains=1L, updateState=FALSE)

  in_block <- FALSE

  name_sampler <- paste0(name, "_sampler_")  # trailing '_' --> not stored by MCMCsim even if store.all=TRUE
  linpred <- function(p) p[[name]]

  make_predict <- function(newdata=NULL, Xnew, verbose=TRUE) {
    stop("out-of-sample prediction for BART model component not yet supported")
  }

  # TODO if sigma.fixed set prior to (almost) fix sigma to 1
  create_BART_sampler <- function() {
    if (e$Q0.type == "diag") {
      dbarts::dbarts(formula=X, data=e$y_eff(),
        weights = e$Q0@x, sigma = if (e$sigma.fixed) 1 else NA_real_,
        control=control, ...
      )
    } else {
      dbarts::dbarts(formula=X, data=e$y_eff(),
        sigma = if (e$sigma.fixed) 1 else NA_real_,
        control=control, ...
      )
    }
  }

  start <- function(p) {
    # TODO multinomial family
    if (is.null(p[[name_sampler]])) {
      p[[name_sampler]] <- create_BART_sampler()
    }
    p[[name]] <- p[[name_sampler]]$run(0L, 1L, FALSE)$train[, 1L]
    p
  }

  rprior <- function(p) {
    if (is.null(p[[name_sampler]])) {
      p[[name_sampler]] <- create_BART_sampler()
    }
    p[[name_sampler]]$sampleTreesFromPrior()
    p[[name_sampler]]$sampleNodeParametersFromPrior()
    p[[name]] <- p[[name_sampler]]$predict(X)
    p
  }

  if (!e$prior.only) {
    draw <- function(p) {}
    if (debug) draw <- add(draw, quote(browser()))
    if (e$e.is.res)
      draw <- add(draw, bquote(p$e_ <- p[["e_"]] + p[[.(name)]]))
    else
      draw <- add(draw, bquote(p$e_ <- p[["e_"]] - p[[.(name)]]))
    draw <- add(draw, bquote(p[[.(name_sampler)]]$setResponse(p[["e_"]])))
    if (e$modeled.Q) {
      draw <- add(draw, bquote(p[[.(name_sampler)]]$setWeights(p[["Q_"]])))
    }
    draw <- add(draw, bquote(p[[.(name_sampler)]]$setSigma(.(if (e$sigma.fixed) 1 else quote(p[["sigma_"]])))))
    draw <- add(draw, bquote(p[[.(name)]] <- p[[.(name_sampler)]]$run(0L, 1L)$train[, 1L]))
    if (e$e.is.res)
      draw <- add(draw, bquote(p$e_ <- p[["e_"]] - p[[.(name)]]))
    else
      draw <- add(draw, bquote(p$e_ <- p[["e_"]] + p[[.(name)]]))
    draw <- add(draw, quote(p))
  }

  environment()
}
