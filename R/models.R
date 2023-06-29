#' Check names of model components
#' 
#' @noRd
#' @param x vector of names of model components.
#' @return \code{TRUE} if all names are OK; throws an error otherwise.
check_mod_names <- function(x) {
  if (any(duplicated(x))) stop("duplicate model component name '", x[duplicated(x)], "'")
  firstlast <- c(substring(x, 1L, 1L), substring(x, nchar(x), nchar(x)))
  if (any("_" == firstlast)) stop("\'_\' as first or last character of a model component name is reserved for internal use")
  # check for names ending with extensions like _sigma, _rho, _xi etc since they might clash
  reserved.exts <- c("_df", "_gl", "_Leroux", "_omega", "_rho", "_sigma", "_xi")
  if (any(grepl(paste0("(", paste(reserved.exts, collapse="|"), ")$"), x)))
    stop("model component names ending with any of (", paste(reserved.exts, collapse=", "), ") are not allowed")
  TRUE
}

#' Compute a list of design matrices for all terms in a model formula,
#' or based on a sampler environment
#'
#' If \code{sampler} is provided instead of \code{formula}, the design matrices
#' are based on the model used to create the sampler environment. In that case, if
#' \code{data} is \code{NULL}, the design matrices stored in \code{sampler} are returned,
#' otherwise the design matrices are computed for the provided data based on the sampler's model.
#' The output is a list of dense or sparse design matrices for the model components
#' with respect to \code{data}.
#'
#' @examples
#' n <- 1000
#' dat <- data.frame(
#'   x = rnorm(n),
#'   f = factor(sample(1:50, n, replace=TRUE))
#' )
#' str(computeDesignMatrix(~ x, dat)[[1]])
#' model <- ~ reg(~x, name="beta") + gen(~x, factor=~f, name="v")
#' X <- computeDesignMatrix(model, dat)
#' str(X)
#'
#' @export
#' @param formula model formula.
#' @param data data frame to be used in deriving the design matrices.
#' @param labels if \code{TRUE}, column names are assigned.
#' @return A list of design matrices.
computeDesignMatrix <- function(formula=NULL, data=NULL, labels=TRUE) {
  if (!inherits(formula, "formula")) stop("'formula' must be a formula")
  formula <- standardize_formula(formula, data=data)
  out <- to_mclist(formula)
  for (m in seq_along(out)) {
    mc <- as.list(out[[m]])[-1L]
    if (is.null(mc$formula))
      mc$formula <- ~ 1
    else
      mc$formula <- as.formula(eval(mc$formula))
    environment(mc$formula) <- environment(formula)
    if (any("|" == all.names(mc$formula))) {
      # assume this is a mec component; we use as design matrix the covariate matrix subject to error
      vs <- as.list(attr(terms(mc$formula), "variables")[-1L])
      formula.X <- as.formula(paste0("~ 0 + ", paste(sapply(vs, function(x) deparse(x[[2L]])), collapse=" + ")), env=environment(formula))
      out[[m]] <- model_matrix(formula.X, data, sparse=mc$sparse)
    } else {
      if (!is.null(mc$factor))
        mc$factor <- as.formula(mc$factor, env=environment(formula))
      if (is.null(mc$remove.redundant)) mc$remove.redundant <- FALSE
      if (is.null(mc$drop.empty.levels)) mc$drop.empty.levels <- FALSE
      out[[m]] <- compute_X(mc$formula, mc$factor, mc$remove.redundant,
          mc$drop.empty.levels, mc$sparse, data)
    }        
    if (!labels) colnames(out[[m]]) <- NULL
  }
  out
}

#' Compute design matrix for a model component
#'
#' The output is a dense or sparse design matrix for the specified model component with respect to
#' the specified data frame. The design matrix contains appropriate column names.
#'
#' @noRd
#' @param formula the formula part of a model component.
#' @param factor the factor part of a model component (also a formula object).
#' @param remove.redundant whether to remove redundant columns.
#' @param drop.empty.levels whether to remove factor levels without observations.
#' @param sparse whether the design matrix based on \code{formula} alone should be in a sparse matrix format.
#' @param data data frame to be used in deriving the design matrix.
#' @return A design matrix including column names.
compute_X <- function(formula=~1, factor=NULL,
                      remove.redundant=TRUE, drop.empty.levels=FALSE,
                      sparse=NULL, data=NULL) {
  X0 <- model_matrix(formula, data, sparse=sparse)
  if (remove.redundant) X0 <- remove_redundancy(X0)
  if (is.null(factor) || (inherits(factor, "formula") && intercept_only(factor)))
    return(X0)
  XA <- compute_XA(factor, data)
  if (drop.empty.levels) {
    cols2remove <- which(zero_col(XA))
    if (length(cols2remove))
      XA <- XA[, -cols2remove]
  } else
    cols2remove <- NULL
  X <- combine_X0_XA(X0, XA)
  if (length(cols2remove))
    attr(X, "factor.cols.removed") <- cols2remove
  if (ncol(X0) > 1L)
    attr(X, "formula.colnames") <- colnames(X0)
  X
}

combine_X0_XA <- function(X0, XA) {
  n <- nrow(X0)
  q0 <- ncol(X0)
  if ((isUnitDiag(XA) || class(XA)[1L] == "tabMatrix") && !large_and_sparse(X0)) {
    fac <- if (isUnitDiag(XA)) seq_len(dim(XA)[1L]) else XA@perm + 1L
    x <- X0
    if (class(XA)[1L] == "tabMatrix") {
      if (XA@num) x <- x * XA@x
      if (XA@reduced) x <- x * as.numeric(XA@perm != -1L)
    }
    out <- drop0(sparseMatrix(
        i = rep.int(seq_len(n), q0),
        j = rep.int(q0 * (fac - 1L) + 1L, q0) + rep_each(0L:(q0 - 1L), n),
        x = as.numeric(x),  # inefficient for large sparse X0
        dims = c(n, q0 * ncol(XA)), check=FALSE
      ), is.Csparse=TRUE
    )
  } else {
    # use Matrix::KhatriRao as it is more memory-efficient
    if (isUnitDiag(XA)) XA <- expandUnitDiag(XA)  # KhatriRao currently fails for unit-ddiMatrix
    if (isUnitDiag(X0)) X0 <- expandUnitDiag(X0)
    out <- t(KhatriRao(t(XA), t(X0)))  # a dgCMatrix
    attr(out@i, "names") <- NULL  # KhatriRao result somehow has names on i and x slots
    attr(out@x, "names") <- NULL
  }
  if (q0 > 1L) {
    if (ncol(XA) * q0 > 1e6L) {
      labs <- NULL
    } else {  # this may be too much; forget about column names
      labs <- paste(rep.int(colnames(X0), ncol(XA)), rep_each(colnames(XA), q0), sep=":")
    }
  } else {
    labs <- colnames(XA)
  }
  attr(out, "Dimnames") <- list(NULL, labs)
  economizeMatrix(out, strip.names=FALSE)
}

compute_XA <- function(factor=NULL, data=NULL, enclos=.GlobalEnv) {
  if (is.null(factor) || intercept_only(factor)) return(NULL)
  factor.info <- get_factor_info(factor, data, enclos=enclos)
  n <- n_row(data)
  if (any("spline" == factor.info$types)) {  # B-spline design matrix components
    fs <- as.list(attr(terms(factor), "variables"))[-1L]
    out <- matrix(1, nrow=1L, ncol=n)
    labs <- ""
    for (f in seq_len(nrow(factor.info))) {
      variable <- eval_in(factor.info$variables[f], data, enclos)
      switch(factor.info$types[f],
        spline = {
          # P-splines, i.e. penalized B-splines
          # RW1/2(t) are special cases with degree 1, knots (tmin+1):(tmax-1), and the same precision matrix
          if (is.null(fs[[f]]$knots)) fs[[f]]$knots <- min(variable) + ((-4:35)/30) * (max(variable) - min(variable))
          if (is.null(fs[[f]]$degree)) fs[[f]]$degree <- 3L
          if (length(fs[[f]]$knots) == 1L) {  # assumed to be the number of knots
            if (fs[[f]]$knots <= 2L*(fs[[f]]$degree + 2L)) stop("spline: too few knots")
            fs[[f]]$knots <- min(variable) + ((seq_len(fs[[f]]$knots) - (fs[[f]]$degree + 1L))/(fs[[f]]$knots - 2L*(fs[[f]]$degree + 2L))) * (max(variable) - min(variable))
          }
          Xf <- drop0(splines::splineDesign(fs[[f]]$knots, variable, fs[[f]]$degree + 1L, sparse=TRUE), tol=sqrt(.Machine$double.eps), is.Csparse=TRUE)
          labs <- as.vector(outer(labs, paste0("bs", seq_len(ncol(Xf))), FUN=paste, sep=if (identical(labs, "")) "" else ":"))
        },
        {
          Xf <- aggrMatrix(variable, facnames=TRUE)
          labs <- as.vector(outer(labs, colnames(Xf), FUN=paste, sep=if (identical(labs, "")) "" else ":"))
        }
      )
      out <- KhatriRao(t(Xf), out)
    }
    out <- economizeMatrix(t(out), strip.names=FALSE)  # names of @i and @x removed by t()
    dimnames(out) <- list(NULL, labs)
    out
  } else {
    fac <- combine_factors(factor.info$variables, data, enclos=enclos)
    if (anyNA(fac)) stop("NA's in 'factor' not allowed")
    aggrMatrix(fac, facnames=TRUE)
  }
}

# extract information from the factor formula component of a model component list
get_factor_info <- function(formula, data, enclos=.GlobalEnv) {
  fs <- as.list(attr(terms(formula), "variables"))[-1L]
  variables <- sapply(fs, function(x) if (is.symbol(x)) deparse(x) else deparse(x[[2L]]))
  types <- sapply(fs, function(x) if (inherits(x, "name")) "iid" else as.character(x[[1L]]))
  ncols <- integer(length(fs))
  for (f in seq_along(ncols)) {
    if (types[f] == "spline") {
      if (is.null(fs[[f]]$knots)) {
        n.knots <- 40L
      } else {
        if (length(fs[[f]]$knots) == 1L)
          n.knots <- fs[[f]]$knots
        else
          n.knots <- length(eval(fs[[f]]$knots))
      }
      if (is.null(fs[[f]]$degree)) fs[[f]]$degree <- 3L
      ncols[f] <- n.knots - fs[[f]]$degree - 1L
    } else {
      if (variables[[f]] == "local_")
        ncols[f] <- n_row(data)
      else
        ncols[f] <- nlevels(as.factor(eval_in(variables[f], data, enclos)))
    }
  }
  data.frame(types=types, variables=variables, n=ncols, stringsAsFactors=FALSE)
}

eval_in <- function(text, data, enclos=.GlobalEnv) {
  if (text == "local_" && all("local_" != colnames(data))) return(seq_len(n_row(data)))
  if (text == "global_" && all("global_" != colnames(data))) return(rep.int(1, n_row(data)))
  if (is.environment(data)) {
    # if data is not a (pair)list/data.frame enclos is not searched, so do it manually if needed
    tryCatch(
      eval(substitute(evalq(expr, data, enclos), list(expr=str2lang_(text)))),
      error = function(e) eval(substitute(evalq(expr, enclos), list(expr=str2lang_(text))))
    )
  } else
    eval(substitute(evalq(expr, data, enclos), list(expr=str2lang_(text))))
}

#' Compute (I)GMRF incidence, precision and restriction matrices corresponding to a generic model component
#'
#' This function computes incidence, precision and restriction matrices, or
#' a subset thereof, for a Gaussian Markov Random Field (GMRF).
#' A GMRF is specified by a formula passed to the \code{factor} argument,
#' in the same way as for the \code{factor} argument of \code{\link{gen}}.
#'
#' @examples
#' n <- 1000
#' dat <- data.frame(
#'   x = rnorm(n),
#'   f1 = factor(sample(1:50, n, replace=TRUE)),
#'   f2 = factor(sample(1:10, n, replace=TRUE))
#' )
#' mats <- compute_GMRF_matrices(~ f1 * RW1(f2), dat)
#' str(mats)
#'
#' @export
#' @param factor factor formula of a generic model component,
#'  see \code{\link{gen}}.
#' @param data data frame to be used in deriving the matrices.
#' @param D if \code{TRUE} compute the incidence matrix.
#' @param Q if \code{TRUE} compute the precision matrix.
#' @param R if \code{TRUE} compute the restriction matrix.
#' @param cols2remove if an integer vector is passed, the dimensions (columns of D,
#'  rows and columns of Q and rows of R) that are removed. This can be useful in the
#'  case of empty domains.
#' @param remove.redundant.R.cols whether to test for and remove redundant restrictions from restriction matrix R
#' @param enclos enclosure to look for objects not found in \code{data}.
#' @param n.parent for internal use; in case of custom factor, the number of frames up
#'   the calling stack in which to evaluate any custom matrices
#' @param ... further arguments passed to \code{economizeMatrix}.
#' @return A list containing some or all of the components \code{D} (incidence matrix),
#'  \code{Q} (precision matrix) and \code{R} (restriction matrix).
compute_GMRF_matrices <- function(factor, data, D=TRUE, Q=TRUE, R=TRUE, cols2remove=NULL,
                                  remove.redundant.R.cols=TRUE, enclos=.GlobalEnv, n.parent=1L, ...) {
  out <- list()
  if (!D && !Q && !R) return(out)
  if (D || Q) DQ <- CdiagU(1L)
  if (R) {
    out$R <- NULL
    nIGMRF <- 0L  # number of (singular) IGMRF factors
  }
  if (is.null(factor)) {
    if (D) out$D <- DQ
    if (Q) out$Q <- DQ
    return(out)
  }
  factor.info <- get_factor_info(factor, data, enclos=enclos)
  fs <- as.list(attr(terms(factor), "variables"))[-1L]
  for (f in seq_len(nrow(factor.info))) {
    fcall <- fs[[f]]
    if (inherits(fcall, "name")) {
      fcall <- call("iid", fcall)  # if no GMRF type name is used assume iid
    }
    fcall[[2L]] <- NULL  # remove first argument, which is assumed to be the factor variable's name
    fcall$n <- factor.info$n[f]
    if (D || Q) {
      DQcall <- fcall
      DQcall[[1L]] <- as.name(paste0(if (D) "D_" else "Q_", DQcall[[1L]]))
    }
    if (R) {
      Rcall <- fcall
      Rcall[[1L]] <- as.name(paste0("R_", Rcall[[1L]]))
    }
    # NB partial matching by `$`, safe enough?
    switch(factor.info$types[f],
      spline = {
        penalty <- fcall$penalty
        if (is.null(penalty)) penalty <- "RW2"
        if (all(penalty != c("iid", "RW1", "RW2"))) stop("unsupported spline penalty argument")
        if (D || Q) {
          DQcall[[1L]] <- as.name(paste0(if (D) "D_" else "Q_", penalty))
          DQcall$knots <- DQcall$degree <- DQcall$penalty <- NULL
          DQf <- eval(DQcall)
        }
        if (R) {
          if (penalty == "iid") {
            Rf <- NULL
          } else {
            Rcall[[1L]] <- as.name(paste0("R_", penalty))
            Rcall$knots <- Rcall$degree <- Rcall$penalty <- NULL
            Rf <- eval(Rcall)
          }
        }
      },
      spatial = {
        poly.df <- if (D || Q) eval(DQcall$poly.df) else eval(Rcall$poly.df)
        if (is.null(poly.df)) stop("'spatial()' requires argument 'poly.df'")
        if (inherits(poly.df, "SpatialPolygonsDataFrame")) {
          if (!requireNamespace("sf", quietly=TRUE)) stop("Package sf required to handle a SpatialPolygonsDataFrame. Please install it.")
          poly.df <- sf::st_as_sf(poly.df)
        }
        if (!inherits(poly.df, "sf")) stop("unexpected input for argument 'poly.df'")
        if (nrow(poly.df) != factor.info$n[f]) stop("number of rows of data slot of '", DQcall$poly.df, "' not equal to number of levels of variable'", factor.info$variables[f], "'")
        # TODO check that level order is the same: first try variable of the same name as factor, then try rownames of @data; if unsuccessful issue a warning
        if (D || Q) {
          DQcall$poly.df <- quote(poly.df)
          DQcall$derive.constraints <- NULL
          DQf <- eval(DQcall)
        }
        if (R) {
          Rcall$poly.df <- quote(poly.df)
          if (isTRUE(Rcall$derive.constraints) && (D || Q))
            Rf <- derive_constraints(if (Q) DQf else crossprod(DQf))
          else
            Rf <- eval(Rcall)
        }
      },
      custom = {
        if (D || Q) {
          if (D) {
            if (is.null(DQcall[["D"]])) stop("missing argument 'D' in custom factor")
            DQf <- economizeMatrix(eval(DQcall[["D"]], parent.frame(n.parent)), sparse=TRUE)
            if (!is.null(DQcall[["Q"]])) warn("argument 'Q' in custom() ignored")
          } else if (is.null(DQcall[["Q"]])) {
            if (is.null(DQcall[["D"]])) stop("custom factor must specify either incidence matrix 'D' or precision matrix 'Q'")
            DQf <- economizeMatrix(crossprod(eval(DQcall[["D"]], parent.frame(n.parent))), sparse=TRUE)
          } else {
            DQf <- economizeMatrix(eval(DQcall[["Q"]], parent.frame(n.parent)), sparse=TRUE, symmetric=TRUE)
          }
          if (ncol(DQf) != factor.info$n[f]) stop("custom factor matrix has unexpected number of columns")
        }
        if (R) {
          if (is.null(Rcall[["R"]])) {
            if (isTRUE(Rcall[["derive.constraints"]])) {
              if (!D && !Q) warn("cannot derive constraints as custom incidence or precision matrix unavailable")
              Rf <- derive_constraints(if (D) crossprod(DQf) else DQf)
            } else {
              Rf <- NULL
            }
          } else {
            Rf <- economizeMatrix(eval(Rcall[["R"]], parent.frame(n.parent)), allow.tabMatrix=FALSE)
            if (dim(Rf)[1L] != factor.info$n[f]) stop("wrong dimension of custom restriction matrix")
            if (isTRUE(Rcall[["derive.constraints"]])) warn("argument 'derive.constraints' ignored as restriction matrix 'R' is specified")
          }
        }
      },
      {
        if (D || Q) DQf <- eval(DQcall)
        if (R) {
          if (any(factor.info$types[f] == c("iid", "AR1")))
            Rf <- NULL
          else
            Rf <- eval(Rcall)
        }
      }
    )
    if (D || Q) DQ <- cross(DQ, DQf)
    if (R) {
      if (!is.null(out$R)) {
        out$R <- kronecker(CdiagU(factor.info$n[f]), out$R)
      }
      if (!is.null(Rf)) {
        nIGMRF <- nIGMRF + 1L
        if (is.null(out$R)) {
          out$R <- zeroMatrix(prod(factor.info$n[seq_len(f)]), 0L)
        }
        out$R <- cbind(out$R, kronecker(Rf, CdiagU(prod(factor.info$n[seq_len(f - 1L)]))))
      }
    }
  }  # END for (f in seq_len(nrow(factor.info)))
  if (D) {
    if (!is.null(cols2remove)) {
      DQ <- DQ[, -cols2remove]
      # then see which rows become all-zero and remove them too
      DQ <- DQ[-which(rowSums(DQ*DQ) == 0), ]
    }
    out$D <- economizeMatrix(DQ, ...)
  }
  if (Q) {
    if (D) {
      out$Q <- crossprod(DQ)
    } else {
      out$Q <- DQ
      if (!is.null(cols2remove)) out$Q <- out$Q[-cols2remove, -cols2remove]
    }
    out$Q <- economizeMatrix(out$Q, symmetric=TRUE, ...)
  }
  if (R && !is.null(out$R)) {
    if (!is.null(cols2remove)) out$R <- out$R[-cols2remove, ]
    if ((remove.redundant.R.cols && nIGMRF >= 2L) || !is.null(cols2remove)) {
      # for multiple IGMRF factors R as constructed has redundant columns,
      #   because of duplicate inclusion of cross-poducts of null-vectors
      out$R <- remove_redundancy(out$R)
    }
    out$R <- economizeMatrix(out$R, allow.tabMatrix=FALSE, ...)
  }
  out
}

#######
# functions to compute sparse precision matrices Q_, incidence matrices D_
#   and for IGMRFs associated constraint matrices R_

#' Correlation structures
#'
#' Element 'factor' of a model component created using function
#' \code{\link{gen}} is a formula composed of several possible terms described
#' below. It is used to derive a (sparse) precision matrix for a set of
#' coefficients, and possibly a matrix representing a set of linear constraints
#' to be imposed on the coefficient vector.
#' \describe{
#'   \item{iid(f)}{Independent effects corresponding to the levels of factor \code{f}.}
#'   \item{RW1(f, circular=FALSE, w=NULL)}{First-order random walk over the levels of factor \code{f}.
#'     The random walk can be made circular and different (fixed) weights can be attached to the innovations.
#'     If specified, \code{w} must be a positive numeric vector of length one less than the number of
#'     factor levels. For example, if the levels correspond to different times, it would often be
#'     reasonable to choose \code{w} proportional to the reciprocal time differences. For equidistant
#'     times there is generally no need to specify \code{w}.}
#'   \item{RW2(f)}{Second-order random walk.}
#'   \item{AR1(f, phi, w=NULL)}{First-order autoregressive correlation structure among
#'     the levels of \code{f}. Required argument is the (fixed) autoregressive parameter \code{phi}.
#'     For irregularly spaced AR(1) processes weights can be specified, in the same way as for
#'     \code{RW1}.}
#'   \item{season(f, period)}{Dummy seasonal with period \code{period}.}
#'   \item{spatial(f, poly.df, snap, queen, derive.constraints=FALSE)}{CAR spatial correlation.
#'     Argument \code{poly.df} can either be an object of (S4) class \code{SpatialPolygonsDataFrame}
#'     or an object of (S3) class \code{sf}. The latter can be obtained, e.g., from reading in a
#'     shape file using function \code{\link[sf]{st_read}}. Arguments \code{snap} and \code{queen}
#'     are passed to \code{\link[spdep]{poly2nb}}.
#'     If \code{derive.constraints=TRUE} the constraint matrix for an IGMRF model component
#'     is formed by computing the singular vectors of the precision matrix.}
#'   \item{spline(f, knots, degree)}{P-splines, i.e. penalized B-splines structure over
#'     the domain of a quantitative variable f. Arguments knots and degree are passed to
#'     \code{\link[splines]{splineDesign}}. If \code{knots} is a single value it is interpreted as
#'     the number of knots, otherwise as a vector of knot positions. By default 40 equally spaced
#'     knots are used, and a degree of 3.}
#'   \item{custom(f, D=NULL, Q=NULL, R=NULL, derive.constraints=NULL)}{Either a custom precision or incidence
#'     matrix associated with factor f can be passed to argument \code{Q} or \code{D}. Optionally a
#'     constraint matrix can be supplied as \code{R}, or constraints can be derived from the null space
#'     of the precision matrix by setting \code{derive.constraints=TRUE}.}
#' }
#' 
#' @examples
#' \donttest{
#' # example of CAR spatial random effects
#' if (requireNamespace("sf")) {
#'   # 1. load a shape file of counties in North Carolina
#'   nc <- sf::st_read(system.file("shape/nc.shp", package="sf"))
#'   # 2. generate some data according to a model with a few regression
#'   # effects, as well as spatial random effects
#'   gd <- generate_data(
#'     ~ reg(~ AREA + BIR74, Q0=1, name="beta") +
#'       gen(factor = ~ spatial(NAME, poly.df=nc), name="vs"),
#'     sigma.mod = pr_invchisq(df=10, scale=1),
#'     data = nc
#'   )
#'   # add the generated target variable and the spatial random effects to the
#'   # spatial dataframe object
#'   nc$y <- gd$y
#'   nc$vs_true <- gd$pars$vs
#'   # 3. fit a model to the generated data, and see to what extent the
#'   #    parameters used to generate the data, gd$pars, are reproduced
#'   sampler <- create_sampler(
#'     y ~ reg(~ AREA + BIR74, Q0=1, name="beta") +
#'     gen(factor = ~ spatial(NAME, poly.df=nc), name="vs"),
#'     block=TRUE, data=nc
#'   )
#'   sim <- MCMCsim(sampler, store.all=TRUE)
#'   (summ <- summary(sim))
#'   nc$vs <- summ$vs[, "Mean"]
#'   plot(nc[c("vs_true", "vs")])
#'   plot(gd$pars$vs, summ$vs[, "Mean"]); abline(0, 1, col="red")
#' }
#' }
#'
#' @name correlation
#' @references
#'  B. Allevius (2018).
#'    On the precision matrix of an irregularly sampled AR(1) process.
#'    arXiv:1801.03791.
#'
#'  H. Rue and L. Held (2005).
#'    Gaussian Markov Random Fields.
#'    Chapman & Hall/CRC.
NULL

D_iid <- Q_iid <- function(n) CdiagU(n)

# use weights w for irregular spacing, see Allevius 2018 research report
Q_AR1 <- function(n, phi, w=NULL) {
  if (n < 2L) stop("AR1 model needs a sequence of at least 2 periods")
  if (is.null(w)) {
    bandSparse(n, n, 0:1, list(c(1, rep.int(1 + phi^2, n-2L), 1), rep.int(-phi, n-1L)), symmetric=TRUE)
  } else {
    # w_i to be interpreted as 1/(t_{i+1} - t_{i})
    w <- as.numeric(w)
    if (any(w <= 0) || length(w) != n-1L) stop("w must be a positive vector of length n-1")
    phi2 <- phi^2
    iw <- 1/w
    d0 <- vector("numeric", n)
    d0[1L] <- (1 - phi2) / (1 - phi2^iw[1L])
    d0[2:(n-1L)] <- (1 - phi2)*(1 - phi2^(iw[1:(n-2L)] + iw[2:(n-1L)])) / ((1 - phi2^iw[1:(n-2L)])*(1 - phi2^iw[2:(n-1L)]))
    d0[n] <- (1 - phi2) / (1 - phi2^iw[n-1L])
    d1 <- - (1 - phi2) * phi^iw / (1 - phi2^iw)
    bandSparse(n, n, 0:1, list(d0, d1), symmetric=TRUE)
  }
}

D_AR1 <- function(n, phi, w=NULL) {
  if (n < 2L) stop("AR1 model needs a sequence of at least 2 periods")
  if (is.null(w)) {
    as(as(bandSparse(n, n, 0:1, list(c(rep.int(1, n-1L), sqrt(1 - phi^2)), rep.int(-phi, n-1L))), "generalMatrix"), "CsparseMatrix")
  } else {
    # w_i to be interpreted as 1/(t_{i+1} - t_{i})
    w <- as.numeric(w)
    if (any(w <= 0) || length(w) != n-1L) stop("w must be a positive vector of length n-1")
    phi2 <- phi^2
    iw <- 1/w
    d0 <- c(sqrt((1 - phi2)/(1 - phi2^iw)), sqrt(1 - phi2))
    d1 <- - sqrt((1 - phi2)/(1 - phi2^iw)) * phi^iw
    as(as(bandSparse(n, n, 0:1, list(d0, d1)), "generalMatrix"), "CsparseMatrix")
  }
}

Q_RW1 <- function(n, circular=FALSE, w=NULL) {
  if (n < 2L) stop("RW1 model needs a sequence of at least 2 periods")
  if (is.null(w)) {
    if (circular)
      bandSparse(n, n, c(0L, 1L, n-1L), list(rep.int(2, n), rep.int(-1, n-1L), -1L), symmetric=TRUE)
    else
      bandSparse(n, n, 0:1, list(c(1, rep.int(2, n-2L), 1), rep.int(-1, n-1L)), symmetric=TRUE)
  } else {
    w <- as.numeric(w)
    if (any(w <= 0) || length(w) != n-1L) stop("w must be a positive vector of length n-1")
    if (circular) {
      bandSparse(n, n, c(0L, 1L, n-1L), list(c(w, w[n-1L]) + c(w[n-1L], w), rep.int(-w, n-1L), -w[n-1L]), symmetric=TRUE)
    } else {
      d0 <- c(w, w[n-1L])
      d0[2:(n-1L)] <- d0[2:(n-1L)] + w[1:(n-2L)]
      bandSparse(n, n, 0:1, list(d0, rep.int(-w, n-1L)), symmetric=TRUE)
    }
  }
}

# incidence matrix for first-order random walk
D_RW1 <- function(n, w=NULL, circular=FALSE) {
  if (n < 2L) stop("RW1 model needs a sequence of at least 2 periods")
  if (is.null(w)) {
    if (circular)
      bandSparse(n, n, c(0L, 1L, 1L-n), list(rep.int(-1, n), rep.int(1, n), 1))
    else
      bandSparse(n-1L, n, 0:1, list(rep.int(-1, n-1L), rep.int(1, n-1L)))
  } else {
    w <- as.numeric(w)
    if (any(w <= 0) || length(w) != n-1L) stop("w must be a positive vector of length n-1")
    sqrt.w <- sqrt(w)
    if (circular)
      bandSparse(n, n, c(0L, 1L, 1L-n), list(-c(sqrt.w, sqrt.w[n-1L]), c(sqrt.w, sqrt.w[n-1L]), sqrt.w[n-1L]))
    else
      bandSparse(n - 1L, n, 0:1, list(-sqrt.w, sqrt.w))
  }
}

R_RW1 <- function(n, w=NULL, circular=FALSE) matrix(1, n, 1L)

Q_RW2 <- function(n) {
  if (n < 4L) stop("RW2 model needs a sequence of at least 4 periods")
  bandSparse(n, n, 0:2, list(c(1, 5, rep.int(6, n-4L), 5, 1), c(-2, rep.int(-4, n-3L), -2), rep.int(1, n-2L)), symmetric=TRUE)
}

D_RW2 <- function(n) {
  if (n < 4L) stop("RW2 model needs a sequence of at least 4 periods")
  # equal to D_RW1(n-1) %*% D_RW1(n)
  bandSparse(n - 2L, n, 0:2, list(rep.int(1, n - 2L), rep.int(-2, n - 2L), rep.int(1, n - 1L)))
}

R_RW2 <- function(n) {
  cbind(rep.int(1, n), seq_len(n))
}

Q_season <- function(n, period) {
  # TODO: check whether this works correctly if the length of the sequence is not a multiple of the period
  if (n < 2L * (period - 1L)) stop(paste("Seasonal component with period", period, "needs a sequence of at least", 2L * (period - 1L), "periods"))
  # create list of bands of sparse precision matrix for Seasonal(nP) component
  band_list <- list()
  for (b in seq_len(period)) {
    band_list <- c(band_list, list(c(seq_len(period - b + 1L), rep.int(period - (b - 1), n - 2L * period + b - 1L), (period - b + 1L):1L)))
  }
  bandSparse(n, n, 0:(period - 1L), band_list, symmetric=TRUE)
}

D_season <- function(n, period) {
  # TODO: check whether this works correctly if the length of the sequence is not a multiple of the period
  if (n < 2L * (period - 1L)) stop(paste("Seasonal component with period", period, "needs a sequence of at least", 2L * (period - 1L), "periods"))
  # create list of bands of sparse precision matrix for Seasonal(nP) component
  band_list <- list()
  for (b in seq_len(period)) {
    band_list <- c(band_list, list(rep.int(1, n - (period - 1L))))
  }
  bandSparse(n - (period - 1L), n, 0:(period - 1L), band_list)
}

R_season <- function(n, period) {
  CM <- matrix(0, n, period - 1L)
  for (b in seq_len(period - 1L)) {
    CM[seq.int(b, n, period), b] <- 1
    CM[seq.int(period, n, period), b] <- -1
  }
  CM
}

get_neighbours_list <- function(n=nrow(poly.df), poly.df, snap=sqrt(.Machine$double.eps), queen=TRUE) {
  # TODO: what if n != length(nb) (out-of-sample levels) ?
  if (n < 2L) stop("Spatial component needs at least 2 areas")
  if (!requireNamespace("spdep", quietly=TRUE)) stop("Package spdep required to construct a spatial precision matrix. Please install it.")
  # derive neighbourhood structure
  spdep::poly2nb(poly.df, snap=snap, queen=queen)
}
  
# CAR spatial precision matrix
Q_spatial <- function(n=nrow(poly.df), poly.df, snap=sqrt(.Machine$double.eps), queen=TRUE) {
  # transform into spatial precision matrix
  # cannot use spdep::nb2mat, as it returns a dense matrix
  nb2Q(get_neighbours_list(n, poly.df, snap, queen))
}

# CAR spatial incidence matrix
D_spatial <- function(n=nrow(poly.df), poly.df, snap=sqrt(.Machine$double.eps), queen=TRUE) {
  nb2D(get_neighbours_list(n, poly.df, snap, queen))
}

R_spatial <- function(n=nrow(poly.df), poly.df, snap=sqrt(.Machine$double.eps), queen=TRUE, derive.constraints=FALSE) {
  if (derive.constraints)
    derive_constraints(Q_spatial(n, poly.df, snap, queen))
  else
    R_RW1(n)
}

# convert neighbourhood structure to a sparse precision matrix
nb2Q <- function(nb) {
  d <- spdep::card(nb)  # numbers of neighbours, diagonal of the precision matrix
  if (sum(d) %% 2L != 0L) warn("not a symmetric neighbourhood structure")
  i <- j <- c(seq_along(nb), rep.int(NA_integer_, sum(d) %/% 2L))  # row and column indices
  r <- length(nb)
  for (g in seq_along(nb)) {
    if (d[g] > 0L) {
      nbg <- nb[[g]]
      nbg <- nbg[nbg > g]  # upper triangular part
      l <- length(nbg)
      if (l > 0L) {
        int <- (r+1L):(r+l)
        i[int] <- g
        j[int] <- nbg
        r <- r + l
      }
    }
  }
  sparseMatrix(i=i, j=j, x=c(d, rep.int(-1, length(i) - length(d))), dims=c(length(nb), length(nb)), symmetric=TRUE)
}

# convert neighbourhood structure to sparse oriented incidence matrix
nb2D <- function(nb) {
  d <- spdep::card(nb)  # numbers of neighbours, diagonal of the precision matrix
  if (sum(d) %% 2L != 0L) warn("not a symmetric neighbourhood structure")
  nr <- sum(d) %/% 2; nc <- length(d)  # dimension of incidence matrix
  i <- rep_each(seq_len(nr), 2L)
  j <- rep.int(NA_integer_, length(i))
  r <- 1L
  for (g in seq_along(nb)) {
    if (d[g] > 0L) {
      nbg <- nb[[g]]
      nbg <- nbg[nbg > g]  # otherwise each edge is counted twice
      for (eg in nbg) {
        j[r] <- g
        j[r + 1L] <- eg
        r <- r + 2L
      }
    }
  }
  sparseMatrix(i=i, j=j, x=rep.int(c(-1, 1), nr), dims=c(nr, nc))
}

#' Derive constraints for, i.e. null space of, a precision matrix
#' 
#' @noRd
#' @param Q l x l precision matrix.
#' @param tol tolerance in regarding a small eigenvalue to be zero.
#' @return l x r Matrix with r the number of singular values.
derive_constraints <- function(Q, tol=.Machine$double.eps^0.5) {
  test <- eigen(Q)
  zeros <- which(test$values < tol)
  if (length(zeros))
    drop0(Matrix(test$vectors[, zeros, drop=FALSE]), tol=tol)
  else
    NULL
}

# test whether a model component is (surely) a proper GMRF
is_proper_GMRF <- function(mc) {
  if (is.null(mc$factor)) return(TRUE)
  Qtypes <- sapply(attr(terms(mc$factor), "variables")[-1L], function(x) if (is.symbol(x)) "iid" else as.character(x[[1L]]))
  (mc$type == "reg") || (mc$Leroux_type != "none") || all(Qtypes %in% c("iid", "AR1"))
}

# fast GMRF prior not possible for spatial as typically n_edges > n_vertices
# fast GMRF prior only if any custom factor has (incidence) matrix D specified
# fast GMRF prior currently only used for IGMRFs
allow_fastGMRFprior <- function(mc) {
  if (is.null(mc$factor)) return(FALSE)
  factors <- attr(terms(mc$factor), "variables")
  Qtypes <- sapply(factors[-1L], function(x) if (is.symbol(x)) "iid" else as.character(x[[1L]]))
  if (mc$Leroux_type != "none" || all(Qtypes %in% c("iid", "AR1"))) return(FALSE)
  if (any(Qtypes %in% "spatial")) return(FALSE)
  if (any(sapply(factors[-1L], function(x) !is.symbol(x) && isTRUE(x$circular)))) return(FALSE)
  for (f in which(Qtypes == "custom"))
    if (is.null(factors[[f + 1L]]$D)) return(FALSE)
  TRUE
}

#' Maximize log-likelihood defined inside a sampler function
#'
#' @examples
#' \donttest{
#' n <- 1000
#' dat <- data.frame(
#'   x = rnorm(n),
#'   f = factor(sample(1:50, n, replace=TRUE))
#' )
#' df <- generate_data(
#'   ~ reg(~x, name="beta", Q0=1) + gen(~x, factor=~f, name="v"),
#'   sigma.fixed=TRUE, data=dat
#' )
#' dat$y <- df$y
#' sampler <- create_sampler(y ~ x + gen(~x, factor=~f, name="v"), data=dat)
#' opt <- maximize_llh(sampler)
#' str(opt)
#' plot(df$par$v, opt$par$v); abline(0, 1, col="red")
#' }
#'
#' @export
#' @param sampler sampler function.
#' @param method optimization method, passed to \code{\link[stats]{optim}}.
#' @param control control parameters, passed to \code{\link[stats]{optim}}.
#' @param ... other parameters passed to \code{\link[stats]{optim}}.
#' @return A list of parameter values that, provided the optimization was successful, maximize the likelihood.
maximize_llh <- function(sampler, method="BFGS", control=list(fnscale=-1), ...) {
  e <- sampler
  if (!is.null(e$Vmod)) stop("optimization not yet implemented for models with sampling variance model specified in 'formula.V'")
  ind_sigma <- e$vec_list[["sigma_"]]
  f <- function(x) {
    if (!is.null(ind_sigma) && x[ind_sigma] < 0)
      -Inf
    else
      e$llh_opt(x)
  }
  res <- optim(par=e$list2vec(sampler$start()), f, method=method, control=control, ...)
  res$par <- e$vec2list(res$par)
  res
}
