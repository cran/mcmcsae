#' Generate draws from the predictive distribution
#'
#' @examples
#' \donttest{
#' n <- 250
#' dat <- data.frame(x=runif(n))
#' dat$y <- 1 + dat$x + rnorm(n)
#' sampler <- create_sampler(y ~ x, data=dat)
#' sim <- MCMCsim(sampler)
#' summary(sim)
#' # in-sample prediction
#' pred <- predict(sim, ppcheck=TRUE)
#' hist(attr(pred, "ppp"))
#' # out-of-sample prediction
#' pred <- predict(sim, newdata=data.frame(x=seq(0, 1, by=0.1)))
#' summary(pred)
#' }
#'
#' @export
## @method predict draws
#' @param object draws object.
#' @param newdata data frame with auxiliary information to be used for prediction.
#  use X. instead of X to avoid argument name clash in parallel version par_predict with parLapply
#' @param X. a list of design matrices; alternatively, \code{X.} equals 'in-sample' or 'linpred'.
#'  If 'in-sample' (the default if newdata is not supplied), the design matrices for in-sample
#'  prediction are used. If 'linpred' the 'linpred_' component of \code{object} is used.
#' @param type the type of predictions. The default is \code{"data"}, meaning that
#'  new data is generated according to the predictive distribution.
#'  If \code{type="link"} only the linear predictor for the mean is generated, and
#'  in case \code{type="response"} the linear predictor is transformed to the response scale.
#'  For Gaussian models \code{type="link"} and \code{type="response"} are equivalent.
#'  For binomial and negative binomial models \code{type="response"} returns the simulations
#'  of the latent probabilities.
#' @param var variance(s) used for out-of-sample prediction. By default 1.
#' @param ny number of trials for used for out-of-sample prediction in case of a binomial model. By default 1.
# use fun. instead of fun to avoid argument name clash with parLapply
#' @param fun. function applied to the vector of posterior predictions to compute one or multiple summaries
#'  or test statistics. The function can have one or two arguments. The first argument is always the vector
#'  of posterior predictions. The optional second argument represents a list of model parameters, needed only
#'  when a test statistic depends on them.
#' @param labels optional names for the output object. Must be a vector of the same length as the result of \code{fun.}.
#' @param ppcheck if \code{TRUE}, function \code{fun.} is also applied to the observed data and
#'  an MCMC approximation is computed of the posterior predictive probability that the test statistic for
#'  predicted data is greater than the test statistic for the observed data.
#' @param iters iterations in \code{object} to use for prediction.
#'  Default \code{NULL} means that all draws from \code{object} are used.
#' @param to.file if \code{TRUE} the predictions are streamed to file.
#' @param filename name of the file to write predictions to in case \code{to.file=TRUE}.
#' @param write.single.prec Whether to write to file in single precision. Default is \code{FALSE}.
#' @param show.progress whether to show a progress bar.
#' @param n.cores the number of cpu cores to use. Default is one, i.e. no parallel computation.
#'  If an existing cluster \code{cl} is provided, \code{n.cores} will be set to the number
#'  of workers in that cluster.
#' @param cl an existing cluster can be passed for parallel computation. If \code{NULL} and
#'  \code{n.cores > 1}, a new cluster is created.
#' @param seed a random seed (integer). For parallel computation it is used to independently
#'  seed RNG streams for all workers.
#' @param export a character vector with names of objects to export to the workers. This may
#'  be needed for parallel execution if expressions in \code{fun.} depend on global variables.
#' @param ... currently not used.
#' @return An object of class \code{dc}, containing draws from the posterior (or prior) predictive distribution.
#'  If \code{ppcheck=TRUE} posterior predictive p-values are returned as an additional attribute.
#'  In case \code{to.file=TRUE} the file name used is returned.
# TODO: pp check only, i.e. without storing predictions either in an object or a file
predict.draws <- function(object, newdata=NULL, X.=if (is.null(newdata)) "in-sample" else NULL,
                          type=c("data", "link", "response"),
                          var=NULL, ny=NULL,
                          fun.=identity, labels=NULL, ppcheck=FALSE, iters=NULL,
                          to.file=FALSE, filename, write.single.prec=FALSE,
                          show.progress=TRUE,
                          n.cores=1L, cl=NULL, seed=NULL, export=NULL, ...) {
  type <- match.arg(type)
  if (ppcheck && type != "data") stop("posterior predictive checks only possible with type = 'data'")
  model <- object[["_model"]]
  fam <- model$family
  par.names <- par_names(object)
  if (is.null(iters))
    iters <- seq_len(object[["_info"]]$n.draw)
  else
    if (!all(iters %in% seq_len(object[["_info"]]$n.draw))) stop("non-existing iterations selected")
  n.chain <- nchains(object)
  cholQ <- NULL
  V <- NULL
  if (is.null(newdata)) {
    if (is.null(X.)) stop("one of 'newdata' and 'X.' must be supplied")
    if (identical(X., "in-sample")) {
      use.linpred_ <- "linpred_" %in% par.names && model$do.linpred && is.null(model$linpred)
      if (!use.linpred_ && !all(names(model$mod) %in% par.names)) stop("for prediction all coefficients must be stored in 'object' (use 'store.all=TRUE' in MCMCsim)")
      X. <- NULL
      n <- model$n
      if (!is.null(var)) {
        warning("argument 'var' ignored for in-sample prediction", immediate.=TRUE)
        var <- NULL
      }
      if (model$modeled.Q) {
        cholQ <- build_chol(
          switch(model$Q0.type,
            unit = rep.int(1, model$n),
            diag = model$Q0@x,
            symm = model$Q0
          )
        )
      } else {
        cholQ <- build_chol(model$Q0)
      }
      if (!is.null(ny)) {
        warning("argument 'ny' ignored for in-sample prediction", immediate.=TRUE)
        ny <- NULL
      }
      # BayesLogit::rpg flags ny=0, so we have set ny to a tiny value > 0 in sampler; need to undo that here as rbinom yields NA for non-integral ny
      if (fam$family == "binomial") ny <- as.integer(round(model$ny))
    } else {
      if (identical(X., "linpred")) {
        if (!("linpred_" %in% par.names))
          stop("'linpred_' not found in object. Use linpred='fitted' in create_sampler.")
        n <- nvars(object[["linpred_"]])
        use.linpred_ <- TRUE
        X. <- NULL
      } else {
        if (!is.list(X.)) stop("'X.' must be a list")
        if (!all(names(X.) %in% names(model$mod)))
          stop("X. names not corresponding to a model component: ", paste(setdiff(names(X.), names(model$mod)), collapse=", "))
        for (i in seq_along(X.)) {
          X.[[i]] <- economizeMatrix(X.[[i]])
          d <- dim(X.[[i]])
          if (i == 1L) 
            n <- d[1L]
          else
            if (d[1L] != n) stop("not all matrices in 'X.' have the same number of rows")
          if (d[2L] != model$mod[[names(X.)[i]]]$q) stop("wrong number of columns for X. component ", names(X.)[i])
        }
        use.linpred_ <- FALSE
      }
      # by default use unit variances for new (oos) predictions, but warn if in-sample variance matrix is not unit diagonal
      if (is.null(var)) {
        if (type == "data" && model$Q0.type != "unit") warning("no 'var' specified, so unit variances are assumed for prediction", immediate.=TRUE)
        var <- 1
      }
      if (!(length(var) %in% c(1L, n))) stop("'var' should be either a scalar value or a vector of length 'nrow(newdata)'")
      if (model$modeled.Q && fam$family == "gaussian") stop("Please use argument 'newdata' instead of 'X.' to also supply information about the factors determining the (modeled) sampling variances.")
      if (fam$family == "binomial") {
        if (is.null(ny)) ny <- model$ny.input
        ny <- check_ny(ny, n, newdata)
      }
    }
    if (!use.linpred_) {
      # check that for mec components name_X is stored
      name_Xs <- unlist(lapply(model$mod, `[[`, "name_X"))
      if (!is.null(X.)) name_Xs <- name_Xs[names(X.)]
      if (!all(name_Xs %in% par.names))
        stop("(in-sample) prediction for a model with measurement error components requires that ",
          "the samples of the corresponding covariates are stored (use 'store.all=TRUE' in MCMCsim)"
        )
    }
  } else {
    if (!is.null(X.)) warning("argument 'X.' is ignored", immediate.=TRUE)
    use.linpred_ <- FALSE
    X. <- list()
    for (mc in model$mod) {
      if (is.null(mc$formula)) stop("model component without 'formula' cannot be used for prediction")
      if (mc$type == "mec") {
        X.[[mc$name]] <- mc$make_predict(newdata)
      } else {
        # for prediction do not (automatically) remove redundant columns!
        X.[[mc$name]] <- compute_X(mc$formula, mc$factor, remove.redundant=FALSE,
                                   drop.empty.levels=FALSE, sparse=mc$sparse, data=newdata)
        if (mc$remove.redundant)
          X.[[mc$name]] <- X.[[mc$name]][, model$coef.names[[mc$name]], drop=FALSE]
        # check that X has the right number of columns
        if (ncol(X.[[mc$name]]) != mc$q) stop("'newdata' yields ", ncol(X.[[mc$name]]), " predictor columns for model term '", mc$name, "' versus ", mc$q, " in the MCMC draws object")
        # TODO check the category names as well; allow oos categories; insert missing categories in newdata
        X.[[mc$name]] <- unname(X.[[mc$name]])
      }
    }
    if (!is.null(X.) && !all(names(X.) %in% par.names)) stop("for prediction all coefficients must be stored in 'object' (use 'store.all=TRUE' in MCMCsim)")
    n <- nrow(newdata)
    if (is.null(var)) {
      if (type == "data" && model$Q0.type != "unit") warning("no 'var' specified, so unit variances are assumed for prediction", immediate.=TRUE)
      var <- 1
    }
    if (!(length(var) %in% c(1L, n))) stop("'var' should be either a scalar value or a vector of length 'nrow(newdata)'")
    if (model$modeled.Q && fam$family == "gaussian") {
      V <- list()
      for (Vmc in model$Vmod) V[[Vmc$name]] <- Vmc$make_predict_Vfactor(newdata)
    }
    if (fam$family == "binomial") {
      if (is.null(ny)) ny <- model$ny.input
      ny <- check_ny(ny, n, newdata)
    }
  }
  fun. <- match.fun(fun.)
  arg1 <- length(formals(fun.))
  if (arg1 < 1L || arg1 > 2L) stop("'fun.' must be a funtion of 1 or 2 arguments")
  arg1 <- arg1 == 1L  # flag for single argument function, i.e. function of prediction x only
  # determine dimension of test statistic
  if (arg1)
    d <- length(fun.(rep.int(0, n)))
  else
    d <- length(fun.(rep.int(0, n), get_draw(object, 1L, 1L)))
  if (!is.null(labels) && length(labels) != d) stop("incompatible 'labels' vector")

  if (is.null(cl))
    n.cores <- min(as.integer(n.cores)[1L], n.chain)
  else
    n.cores <- length(cl)
  if (n.cores > 1L) {
    if (n.cores > parallel::detectCores()) stop("too many cores")
    reserved.export.names <- c("X.", "cholQ", "var", "V", "ny", "d", "ppcheck")
    if (is.null(cl)) {
      cl <- setup_cluster(n.cores, seed, export)
      on.exit(stop_cluster(cl))
    } else {
      if (!is.null(export)) {
        if (!is.character(export)) stop("'export' must be a character vector")
        if (any(export %in% reserved.export.names))
          stop("please use variable names other than ", paste(reserved.export.names, collapse=", "),
               " in the definition of 'fun.'")
        parallel::clusterExport(cl, export)
      }
      if (!is.null(seed)) parallel::clusterSetRNGStream(cl, seed)
    }
    parallel::clusterExport(cl, reserved.export.names, environment())
    message(length(iters), " draws distributed over ", n.cores, " cores")
    sim_list <- split_iters(object, iters, parts=n.cores)
    predict_obj <- function(obj) {
      chains <- seq_len(nchains(obj))
      iters <- seq_len(ndraws(obj))
      out <- rep.int(list(matrix(NA_real_, length(iters), d)), length(chains))
      if (ppcheck) {
        ppp <- rep.int(0, d)
        if (arg1)
          fy <- fun.(model$y)
        else
          y <- obj[["_model"]]$y
      }
      for (i in iters) {
        for (ch in chains) {
          p <- get_draw(obj, i, ch)
          if (use.linpred_)
            ystar <- p[["linpred_"]]
          else
            ystar <- model$lin_predict(p, X.)
          if (type == "data")
            ystar <- model$rpredictive(p, ystar, cholQ, var, V, ny)
          else if (type == "response")
            ystar <- fam$linkinv(ystar)
          out[[ch]][i, ] <- if (arg1) fun.(ystar) else fun.(ystar, p)
          if (ppcheck) ppp <- ppp + (out[[ch]][i, ] >= if (arg1) fy else fun.(y, p))
        }
      }
      out
    }
    message("running predict")
    out <- combine_iters_dc(parallel::parLapply(cl, X=sim_list, fun=predict_obj))
  } else {
    if (!is.null(seed)) set.seed(seed)
    n.it <- length(iters)
    if (to.file) {
      if (missing(filename)) filename <- "MCMCdraws_pred.dat"
      outfile <- file(filename, "wb")
      write_header(outfile, n.it, n.chain, d, labels, write.single.prec)
      write.size <- if (write.single.prec) 4L else NA_integer_
      ppcheck <- FALSE  # TODO allow ppcheck in the case that predictions are written to file
    } else {
      out <- rep.int(list(matrix(NA_real_, n.it, d)), n.chain)
    }
    if (ppcheck) {
      ppp <- rep.int(0, d)
      if (arg1)
        fy <- fun.(model$y)
      else
        y <- model$y
    }
    r <- 1L
    show.progress <- show.progress && n.it > 1L
    if (show.progress) pb <- txtProgressBar(min=1L, max=n.it, style=3L)
    for (i in iters) {
      for (ch in seq_len(n.chain)) {
        p <- get_draw(object, i, ch)
        if (use.linpred_)
          ystar <- p[["linpred_"]]
        else
          ystar <- model$lin_predict(p, X.)
        if (type == "data")
          ystar <- model$rpredictive(p, ystar, cholQ, var, V, ny)
        else if (type == "response")
          ystar <- fam$linkinv(ystar)
        if (to.file)
          writeBin(if (arg1) fun.(ystar) else fun.(ystar, p), con=outfile, size=write.size)
        else
          out[[ch]][r, ] <- if (arg1) fun.(ystar) else fun.(ystar, p)
        if (ppcheck) ppp <- ppp + (out[[ch]][r, ] >= if (arg1) fy else fun.(y, p))
      }
      if (show.progress) setTxtProgressBar(pb, r)
      r <- r + 1L
    }
    if (show.progress) close(pb)
    if (to.file) {
      close(outfile)
      return(filename)
    }
  }
  if (ppcheck) attr(out, "ppp") <- ppp/(n.chain * length(iters))
  if (!is.null(labels)) attr(out, "labels") <- labels
  class(out) <- "dc"
  out
}

#' Generate a data vector according to a model
#'
#' This function generates draws from the prior predictive distribution.
#' Parameter values are drawn from their priors, and consequently data is
#' generated from the sampling distribution given these parameter values.
#'
#' @examples
#' \donttest{
#' n <- 250
#' dat <- data.frame(
#'   x = rnorm(n),
#'   g = factor(sample(1:10, n, replace=TRUE)),
#'   ny = 10
#' )
#' gd <- generate_data(
#'   ~ reg(~ 1 + x, Q0=10, b0=c(0, 1), name="beta") + gen(factor = ~ g, name="v"),
#'   family="binomial", ny="ny", data=dat
#' )
#' gd
#' plot(dat$x, gd$y)
#' }
#'
#' @export
#' @param formula A model formula, see \code{create_sampler}.
#'  Any left-hand-side of the formula is ignored.
#' @param data see \code{create_sampler}.
#' @param family see \code{create_sampler}.
#' @param ny see \code{create_sampler}.
#' @param ry see \code{create_sampler}.
#' @param r.mod see \code{create_sampler}.
#' @param sigma.fixed see \code{create_sampler}.
#' @param sigma.mod see \code{create_sampler}.
#' @param Q0 see \code{create_sampler}.
#' @param formula.V see \code{create_sampler}.
#' @param linpred see \code{create_sampler}.
#' @return A list with a generated data vector and a list of prior means of the
#'  parameters. The parameters are drawn from their priors.
generate_data <- function(formula, data=NULL, family="gaussian",
                          ny=NULL, ry, r.mod,
                          sigma.fixed=(family != "gaussian"),
                          sigma.mod=NULL, Q0=NULL, formula.V=NULL,
                          linpred=NULL) {
  if (!sigma.fixed && is.null(sigma.mod))
    sigma.mod <- pr_invchisq(df=1, scale=1, n=1L)  # use a proper prior for sigma by default
  sampler <- create_sampler(formula=formula, data=data, family=family, ny=ny, ry=ry, r.mod=r.mod,
                            sigma.fixed=sigma.fixed, sigma.mod=sigma.mod, Q0=Q0, formula.V=formula.V,
                            linpred=linpred, prior.only=TRUE)
  sim <- MCMCsim(sampler, n.iter=1L, n.chain=1L, store.all=TRUE, from.prior=TRUE, verbose=FALSE)
  pars <- lapply(sim[par_names(sim)], function(x) setNames(x[[1L]][1L, ], attr(x, "labels")))
  list(y=as.vector(as.matrix.dc(predict(sim), colnames=FALSE)), pars=pars)
}
