
# memory-efficient version of model.matrix
# does not first compute a full model frame
# simplifications:
# - any missing value raises an error (cf. na.action=fail in model.frame)
# - only treatment constrasts supported (i.p. no polynomial contrasts for ordered factors)
# very efficient in sparse case
# optionally uses consistent naming, e.g., var1$cat1:var2$cat2
# forbid use of category separator and : in variable names to guarantee unique column names
# option to aggregate by a factor variable (term by term, i.e. in a memory-efficient way)

# TODO
# - composite matrix as list of tabMatrices
# - contrasts: "contr.sparse" for removing level with most observations

#' Compute possibly sparse model matrix, optionally removing redundant columns
#'
#' @noRd
#' @param formula model formula; only the rhs of it is currently used.
#' @param data unit-level data frame containing all variables used in \code{formula}.
#'  These variables should not contain missing values. An error is raised in case they do.
#' @param contrasts.arg specification of contrasts for factor variables. Passed to \code{\link{model_matrix}}.
#' @param remove.redundant if \code{TRUE} redundant columns in the design matrix are removed.
#' @param redundancy.tol tolerance for detecting linear dependencies among the columns of the design matrix.
#'  If \code{NULL}, a default value is used.
#' @param sparse if \code{TRUE} a sparse is returned. This can be efficient for large datasets
#'  and a model containing categorical variables with many categories. If \code{sparse=NULL},
#'  the default, a simple heuristic determines whether a sparse or dense model matrix is returned.
#' @return Predictor matrix X, either an ordinary matrix or a sparse \code{dgCMatrix}.
model_Matrix <- function(formula, data=environment(formula), contrasts.arg=NULL, remove.redundant=TRUE, redundancy.tol=NULL, sparse=NULL) {
  X <- model_matrix(formula, data, contrasts.arg, sparse=sparse)
  if (remove.redundant) {
    # first remove trivial redundancy: zero columns
    zerocols <- which(zero_col(X))  # all-0 columns
    if (length(zerocols)) {
      X <- X[, -zerocols, drop=FALSE]
    }
    # remove remaining redundancy
    X <- remove_redundancy(X, tol=redundancy.tol)
  }
  X
}

#' Compute possibly sparse model matrix
#'
#' @export
#' @param formula model formula.
#' @param data data frame containing all variables used in \code{formula}.
#'  These variables should not contain missing values. An error is raised in case any of them does.
#' @param contrasts.arg specification of contrasts for factor variables. Currently supported are
#'  "contr.none" (no contrasts applied), "contr.treatment" (first level removed) and "contr.SAS" (last level removed).
#'  Alternatively, a named list specifying a single level per factor variable can be passed.
#' @param drop.unused.levels whether empty levels of individual factor variables should be removed.
#' @param sparse if \code{TRUE} a sparse matrix of class \code{dgCMatrix} is returned. This can be efficient
#'  for large datasets and a model containing categorical variables with many categories. If \code{sparse=NULL}, the default,
#'  whether a sparse or dense model matrix is returned is based on a simple heuristic.
#' @param drop0 whether to drop any remaining explicit zeros in resulting sparse matrix.
#' @param catsep separator for concatenating factor variable names and level names.
#'  By default it is the empty string, reproducing the labels of \code{model.matrix}.
#' @param by a vector by which to aggregate the result.
#' @param tabM if \code{TRUE} return a list of tabMatrix objects.
#' @return Design matrix X, either an ordinary matrix or a sparse \code{dgCMatrix}.
model_matrix <- function(formula, data=environment(formula), contrasts.arg=NULL,
                         drop.unused.levels=FALSE, sparse=NULL, drop0=TRUE, catsep="",
                         by=NULL, tabM=FALSE) {
  if (missing(data)) {
    trms <- terms(formula)
    if (!NROW(attr(trms, "factors"))) stop("cannot determine the sample size")
    n <- NROW(eval_in(rownames(attr(trms, "factors"))[1L], data))
  } else {
    trms <- terms(formula, data=data)
    n <- nrow(data)
  }
  if (!is.null(by)) {
    if (length(by) != n) stop("incompatible vector 'by'")
    Maggr <- aggrMatrix(by)
  }
  has_intercept <- attr(trms, "intercept") == 1L
  tmat <- terms_matrix(trms)
  if (!length(tmat)) {
    if (has_intercept) {
      if (isTRUE(sparse))
        out <- new("dgCMatrix", i=0:(n - 1L), p=c(0L, n), x=rep.int(1, n), Dim=c(n, 1L), Dimnames=list(NULL, "(Intercept)"))
      else
        out <- matrix(1, n, 1L, dimnames=list(NULL, "(Intercept)"))
      if (is.null(by))
        return(out)
      else
        return(crossprod(Maggr, out))
    } else {
      # empty design matrix
      if (isTRUE(sparse))
        return(zeroMatrix(if (is.null(by)) n else ncol(by), 0L))
      else
        return(matrix(0, nrow=if (is.null(by)) n else ncol(by), ncol=0L))
    }
  }
  tnames <- colnames(tmat)
  vnames <- rownames(tmat)
  # 1. analyse
  # which variables are quantitative
  qvar <- !catvars(trms, data)
  qterm <- apply(tmat, 2L, function(x) any(qvar[x == 1L]))
  qvar <- vnames[which(qvar)]
  q <- if (has_intercept) 1L else 0L  # nr of columns
  qd <- q  # equivalent nr of dense columns, for estimation of sparseness
  if (!is.list(contrasts.arg)) {
    contrasts.arg <- match.arg(contrasts.arg, c("contr.treatment", "contr.SAS", "contr.none"))
    if (contrasts.arg == "contr.none") contrasts.arg <- NULL
  }
  contr.list <- NULL
  first <- TRUE
  # loop over terms to determine nr of columns, and store levels to be removed
  for (k in seq_along(tnames)) {
    nc <- 1L  # number of cells (or columns in design matrix)
    # loop over quantitative variables, including possibly matrices
    for (v in intersect(vnames[tmat[, k] > 0L], qvar)) {
      qv <- eval_in(v, data)
      if (NROW(qv) != n) stop("variable lengths differ")
      if (anyNA(qv)) stop("missing(s) in variable ", v)
      if (is.matrix(qv)) nc <- nc * ncol(qv)
    }
    qd <- qd + nc
    rmlevs <- NULL
    # loop over the factor vars and compute number of levels accounting for removed cols
    for (f in setdiff(vnames[tmat[, k] > 0L], qvar)) {
      fac <- eval_in(f, data)
      if (!is.factor(fac) || drop.unused.levels) fac <- factor(fac)
      if (length(fac) != n) stop("variable lengths differ")
      if (anyNA(fac)) stop("missing(s) in variable ", f)
      levs <- attr(fac, "levels")
      ncat <- length(levs)
      if (ncat <= 1L) stop("factor variable with fewer than two levels: ", f)
      # NB if tmat[f, k] == 2L a main effect is missing and we should not remove a level!
      # also, if there is no intercept, no level should be removed for the first main factor effect
      # see the documentation of terms.formula
      if (!has_intercept && first && attr(trms, "order")[k] == 1L) {
        tmat[f, k] <- 2L  # next time we can just check for value 2, meaning not to remove a level
        first <- FALSE
      }
      # store levels to remove to be passed to fac2tabM
      if (!is.null(contrasts.arg) && tmat[f, k] == 1L) {
        if (is.list(contrasts.arg)) {
          rml <- contrasts.arg[[f]]
          if (is.null(rml))
            rml <- levs[1L]
          else
            if (!rml %in% levs) stop("invalid contrasts.arg: cannot remove level '", rml, "' from factor ", f)
        } else if (contrasts.arg == "contr.treatment") {
          rml <- levs[1L]  # remove first category
        } else if (contrasts.arg == "contr.SAS") {
          rml <- levs[ncat]  # remove last category
        }
        ncat <- ncat - 1L
        rmlevs <- c(rmlevs, setNames(rml, f))
      }
      nc <- nc * ncat
    }
    if (!is.null(rmlevs)) {
      contr.list <- c(contr.list, setNames(list(rmlevs), tnames[k]))
    }
    q <- q + nc
    # TODO also count nonzeros, i.e. the zeros in quantitative variables!
  }
  if (catsep == ":") stop("':' is not allowed as category separator in column labels")
  if (is.null(sparse)) {
    sparse <- qd < 0.5 * q  # determine whether sparse
  }
  # 2. construct
  if (!is.null(by)) {
    n <- ncol(Maggr)
  }
  if (sparse) {
    # allocate memory for i, x slots
    i <- integer(n * (ncol(tmat) + has_intercept))
    x <- numeric(length(i))
    p <- rep.int(0L, q + 1L)
  } else {
    # allocate memory for dense model matrix
    out <- matrix(NA_real_, n, q)
  }
  lab <- character(q)
  # loop over terms and fill in i,x slots of sparse model matrix
  if (has_intercept) {
    lab[1L] <- "(Intercept)"
    xk <- if (is.null(by)) 1 else colSums(Maggr)
    if (sparse) {
      i[seq_len(n)] <- 0:(n - 1L)
      x[seq_len(n)] <- xk
      p[2L] <- n
      at <- n + 1L
    } else {
      out[, 1L] <- xk
    }
    col <- 2L
  } else {
    if (sparse) at <- 1L
    col <- 1L
  }
  for (k in seq_along(tnames)) {
    countvars <- intersect(vnames[tmat[, k] > 0L], qvar)
    if (length(countvars)) {
      xk <- eval_in(countvars[1L], data)
      if (is.matrix(xk))
        labk <- paste(countvars[1L], if (is.null(colnames(xk))) seq_len(ncol(xk)) else colnames(xk), sep=catsep)
      else
        labk <- countvars[1L]
      for (v in countvars[-1L]) {
        temp <- eval_in(v, data)
        if (is.matrix(temp)) {
          xk <- t(KhatriRao(t(xk), t(temp)))
          labk <- outer(labk, paste(v, colnames(xk), sep=catsep), paste, sep=":")
        } else {
          xk <- xk * temp
          labk <- paste(labk, v, sep=":")
        }
      }
    } else {
      xk <- NULL
      labk <- NULL
    }
    facvars <- setdiff(vnames[tmat[, k] > 0L], qvar)
    if (length(facvars)) {
      if (length(countvars) && !is.matrix(xk)) {
        # TODO allow matrix; --> generalize x-slot in tabMatrix to matrix (and even dgCMatrix)
        fk <- fac2tabM(facvars, data, x=xk, contrasts=contr.list[[tnames[k]]], catsep=catsep)
      } else {
        fk <- fac2tabM(facvars, data, contrasts=contr.list[[tnames[k]]], catsep=catsep)
      }
      if (is.matrix(xk)) {
        lab[col:(col + ncol(fk)*ncol(xk) - 1L)] <- outer(colnames(fk), labk, paste, sep=":")
        fk <- t(KhatriRao(t(xk), t(fk)))  # col-index of fk runs fastest
      } else {
        lab[col:(col + ncol(fk) - 1L)] <- paste0(labk, if (!is.null(labk)) ":" else "", colnames(fk))
      }
      if (!is.null(by)) {
        fk <- crossprod(Maggr, fk)
      }
      if (sparse) {
        if (class(fk)[1L] != "dgCMatrix") {
          if (class(fk)[1L] == "tabMatrix")
            fk <- Ctab2dgC(fk)
          else
            fk <- as(fk, "CsparseMatrix")
        }
        if (l <- length(fk@i)) {
          i[at:(at + l - 1L)] <- fk@i
          x[at:(at + l - 1L)] <- fk@x
          at <- at + l
        }
        p[(col + 1L):(col + ncol(fk))] <- p[col] + fk@p[-1L]
      } else {
        if (class(fk)[1L] == "tabMatrix")
          out[, col:(col + ncol(fk) - 1L)] <- Ctab2mat(fk)
        else
          out[, col:(col + ncol(fk) - 1L)] <- as(fk, "matrix")
      }
      col <- col + ncol(fk)
    } else {
      if (is.matrix(xk)) {
        lab[col:(col + ncol(xk) - 1L)] <- labk
        if (!is.null(by)) {
          xk <- crossprod(Maggr, xk)
        }
        if (sparse) {
          i[at:(at + length(xk) - 1L)] <- rep.int(0:(n - 1L), ncol(xk))
          x[at:(at + length(xk) - 1L)] <- xk
          p[(col + 1L):(col + ncol(xk))] <- seq.int(at + n - 1L, by=n, length.out=ncol(xk))
          at <- at + length(xk)
        } else {
          out[, col:(col + ncol(xk) - 1L)] <- xk
        }
        col <- col + ncol(xk)
      } else {
        lab[col] <- labk
        if (!is.null(by)) {
          xk <- crossprod_mv(Maggr, xk)
        }
        if (sparse) {
          i[at:(at + length(xk) - 1L)] <- 0:(n - 1L)
          x[at:(at + length(xk) - 1L)] <- xk
          at <- at + length(xk)
          p[col + 1L] <- at - 1L
        } else {
          out[, col] <- xk
        }
        col <- col + 1L
      }
    }
  }
  if (sparse) {
    i <- i[seq_len(at - 1L)]
    x <- x[seq_len(at - 1L)]
    out <- new("dgCMatrix", i=i, p=p, x=x, Dim=c(n, q), Dimnames=list(NULL, lab))
    if (drop0) out <- drop0(out, is.Csparse=TRUE)
  } else {
    dimnames(out) <- list(NULL, lab)
  }
  out
}


# return a named logical vector indicating for each variable in
#   the model terms object whether it is categorical
catvars <- function(trms, data) {
  vnames <- rownames(terms_matrix(trms))
  vapply(vnames, function(x) {
      temp <- eval_in(x, data)
      is.factor(temp) || is.character(temp)
    }, FALSE
  )
}

# extract the independent variable in terms matrix from a terms object
terms_matrix <- function(trms) {
  tmat <- attr(trms, "factor")
  if (attr(trms, "response")) {
    w <- which(rowSums(tmat) == 0)
    if (length(w)) tmat <- tmat[-w, , drop=FALSE]
  }
  tmat
}
