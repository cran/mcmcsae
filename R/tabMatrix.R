# tabMatrix class: generalization of indMatrix from Matrix package
setClass("tabMatrix",
  slots = c(
    perm = "integer",  # as in indMatrix, but 0-based
    reduced = "logical",  # whether columns are removed (contrasts); a -1 in perm represents such columns
    # if reduced there can be all-0 rows
    num = "logical",  # whether the matrix is numeric or 0/1
    # TODO generalize x to matrix/dgCMatrix
    x = "numeric",  # the numeric values in case num=TRUE
    xlab = "character"  # name asociated with numeric x component (scalar, later to be generalized to colnames)
  ),
  contains = c("sparseMatrix", "generalMatrix"),
  validity = function(object) {
    n <- object@Dim[1L]
    d <- object@Dim[2L]
    perm <- object@perm
    num <- object@num
    reduced <- object@reduced
    xlab <- object@xlab
    if (length(perm) != n) return(paste("length of 'perm' slot must be", n))
    if (num && length(object@x) != n) return(paste("length of 'x' slot must be", n))
    if (length(xlab) > 1L) return("length of 'xlab' must not be greater than 1")
    if (n > 0L && (any(perm > d - 1L) || (!reduced && any(perm < 0L)))) return("'perm' slot is not a valid index")
    TRUE
  }
)

# general constructor
tabMatrix <- function(Dim, Dimnames=NULL, reduced=FALSE, perm=integer(), num=FALSE, x=numeric(), xlab=NULL) {
  out <- Ctab(Dim=Dim, reduced=reduced, perm=perm, num=num, x=x)
  if (!is.null(Dimnames)) attr(out, "Dimnames") <- Dimnames
  if (!is.null(xlab)) attr(out, "xlab") <- xlab
  out
}

# constructor from factor
setAs("factor", "tabMatrix",
  function(from) {
    levs <- attr(from, "levels")
    tabMatrix(
      Dim = c(length(from), length(levs)),
      Dimnames = list(NULL, levs),
      reduced = FALSE,
      perm = as.integer(from) - 1L,
      num = FALSE,
      x = numeric(0L),
      xlab = character(0L)
    )
  }
)

# colSums method for tabMatrix (colSums promoted to S4 generic by Matrix)
setMethod("colSums", "tabMatrix", function(x, na.rm=FALSE, dims=1, ...)
  if (x@num)
    fast_aggrC(x@x, x@perm + 1L, x@Dim[2L])
  else
    tabulate(x@perm + 1L, nbins=x@Dim[2L])
)

setMethod("rowSums", "tabMatrix", function(x, na.rm=FALSE, dims=1, ...) {
  if (x@num)
    out <- x@x
  else
    out <- rep.int(1, x@Dim[1L])
  if (x@reduced) out[x@perm == -1L] <- 0
  out
})

setMethod("isDiagonal", "tabMatrix", function(object) {
  d <- object@Dim
  if (d[1L] != d[2L]) return(FALSE)
  if (object@reduced) {
    sub <- which(object@perm != -1L)
    identical(object@perm[sub], (0:(d[1L] - 1L))[sub])
  } else {
    identical(object@perm, 0:(d[1L] - 1L))
  }
})

setMethod("diag", "tabMatrix", function(x=1, nrow, ncol) {
  d <- x@Dim
  if (d[1L] <= d[2L]) {
    if (x@num)
      out <- ifelse(x@perm == 0:(x@d[1L] - 1L), x@x, 0)
    else
      out <- ifelse(x@perm == 0:(x@d[1L] - 1L), 1, 0)
  } else {
    if (x@num)
      out <- ifelse(x@perm[seq_len(d[2L])] == 0:(d[2L] - 1L), x@x[seq_len(d[2L])], 0)
    else
      out <- ifelse(x@perm[seq_len(d[2L])] == 0:(d[2L] - 1L), 1, 0)
  }
  if (x@reduced) out[x@perm == -1L] <- 0
  out
})

setMethod("t", "tabMatrix", function(x) {
  if (x@reduced) {
    sub <- which(x@perm != -1L)
    if (x@num)
      sparseMatrix(i=x@perm[sub], j=(0:(x@Dim[1L] - 1L))[sub], x=x@x[sub], dims=rev(x@Dim), dimnames=rev(x@Dimnames), index1=FALSE, check=FALSE)
    else
      sparseMatrix(i=x@perm[sub], j=(0:(x@Dim[1L] - 1L))[sub], x=rep.int(1, length(sub)), dims=rev(x@Dim), dimnames=rev(x@Dimnames), index1=FALSE, check=FALSE)
  } else {
    if (x@num)
      sparseMatrix(i=x@perm, j=0:(x@Dim[1L] - 1L), x=x@x, dims=rev(x@Dim), dimnames=rev(x@Dimnames), index1=FALSE, check=FALSE)
    else
      sparseMatrix(i=x@perm, j=0:(x@Dim[1L] - 1L), x=rep.int(1, x@Dim[1L]), dims=rev(x@Dim), dimnames=rev(x@Dimnames), index1=FALSE, check=FALSE)
  }
})

setMethod("coerce", c("tabMatrix", "ddiMatrix"), function(from, to="ddiMatrix", strict=TRUE) {
  if (!isDiagonal(from)) stop("not a diagonal tabMatrix")
  if (from@reduced) {
    x <- if (from@num) from@x else rep.int(1, from@Dim[1L])
    x[from@perm == -1L] <- 0
    out <- Cdiag(x)
  } else {
    if (from@num)
      out <- Cdiag(from@x)
    else
      out <- CdiagU(from@Dim[1L])
  }
  if (!(is.null(from@Dimnames[[1L]]) && is.null(from@Dimnames[[2L]]))) attr(out, "Dimnames") <- from@Dimnames
  out
})

setMethod("coerce", c("tabMatrix", "matrix"), function(from, to="matrix", strict=TRUE) {
  out <- Ctab2mat(from)
  if (!(is.null(from@Dimnames[[1L]]) && is.null(from@Dimnames[[2L]]))) attr(out, "dimnames") <- from@Dimnames
  out
})

setMethod("coerce", c("tabMatrix", "CsparseMatrix"), function(from, to="CsparseMatrix", strict=TRUE) {
  out <- Ctab2dgC(from)
  if (!(is.null(from@Dimnames[[1L]]) && is.null(from@Dimnames[[2L]]))) attr(out, "Dimnames") <- from@Dimnames
  out
})

setMethod("coerce", c("ddiMatrix", "tabMatrix"), function(from, to="tabMatrix", strict=TRUE)
  tabMatrix(Dim=from@Dim, Dimnames=from@Dimnames,
    reduced=FALSE, perm=0:(from@Dim[1L] - 1L), num = from@diag == "N",
    x = if (from@diag == "N") from@x else numeric()
  )
)

setMethod("coerce", c("matrix", "tabMatrix"), function(from, to="tabMatrix", strict=TRUE) {
  perm <- apply(from, 1L, function(x) which(x != 0))
  if (is.integer(perm)) {
    x <- from[cbind(seq_len(nrow(from)), perm)]
    num <- any(x != 1)
    if (!num) x <- numeric()
    reduced <- FALSE
  } else {
    if (!is.list(perm)) stop("not a tabMatrix")
    l <- sapply(perm, length)
    if (any(!(l %in% 0:1))) stop("not a tabMatrix")
    l <- l == 1L
    perm[!l] <- 0L
    perm <- as.integer(perm)
    reduced = TRUE
    xnz <- from[cbind(seq_len(nrow(from)), perm)[l, , drop=FALSE]]
    num <- any(xnz != 1)
    if (num) {
      x <- rep.int(0, length(perm))
      x[l] <- xnz
    } else {
      x <- numeric()
    }
  }
  tabMatrix(Dim=dim(from), Dimnames=dimnames(from),
    reduced=reduced, perm=perm - 1L, num=num, x=x
  )
})

setMethod("coerce", c("dgCMatrix", "tabMatrix"), function(from, to="tabMatrix", strict=TRUE) {
  perm <- rep.int(-1L, nrow(from))
  x <- vector("double", nrow(from))
  nnz <- from@p[-1L] - from@p[-length(from@p)]
  ind <- 1L
  for (co in seq_along(nnz)) {
    if (nnz[co]) {
      ix_ind <- ind:(ind + nnz[co] - 1L)
      perm[from@i[ix_ind][from@x[ix_ind] != 0] + 1L] <- co - 1L
      x[from@i[ix_ind][from@x[ix_ind] != 0] + 1L] <- from@x[ix_ind][from@x[ix_ind] != 0]
      ind <- ind + nnz[co]
    }
  }
  if (all(x == 1)) {
    x <- numeric(0L)
    num <- FALSE
  } else {
    num <- TRUE
  }
  tabMatrix(Dim=from@Dim, Dimnames=from@Dimnames,
    reduced=any(perm == -1L), perm=perm, num=num, x=x
  )
})

tab_is_zero <- function(x) (x@reduced && all(x@perm == -1L)) || (x@num && all(x@x == 0))

dgC_is_tabMatrix <- function(M) anyDuplicated(M@i[M@x != 0]) == 0

# see whether tabMatrix is actually a permutation matrix
tab_isPermutation <- function(M) {
  M@Dim[1L] == M@Dim[2L] && !M@num && identical(sort(M@perm), 0:(M@Dim[1L] - 1L))
}

#' S4 method for row and column subsetting a 'tabMatrix'
#'
#' @keywords internal
#' @param x a tabMatrix object.
#' @param i integer vector indicating the rows to select. Not used for column subsetting.
#' @param j integer vector indicating the columns to select. Not used for row subsetting.
#' @param ... not used.
#' @param drop whether to return a vector in case of a single selected row or column.
#' @return the selected rows/columns as a tabMatrix or a vector.
#' @name tabMatrix-indexing
NULL

# row selection
#' @rdname tabMatrix-indexing
setMethod("[", c(x="tabMatrix", i="index", j="missing", drop="logical"), function (x, i, j, ..., drop=TRUE) {
  if (is.character(i)) stop("indexing by name currently not supported for tabMatrix")
  i <- as.integer(i)
  if (any(i < 0L)) stop("negative integer indexing currently not supported for tabMatrix")
  if (any(i == 0L | i > x@Dim[1L])) stop("out of range")
  if (drop && length(i) == 1L) {
    out <- rep.int(0, x@Dim[2L])
    if (x@perm[i] >= 0L)
      out[x@perm[i] + 1L] <- if (x@num) x@x[i] else 1
    out
  } else {
    tabMatrix(Dim=c(length(i), x@Dim[2L]), Dimnames=list(x@Dimnames[[1L]][i], x@Dimnames[[2L]]),
      reduced=x@reduced, perm=x@perm[i], num=x@num, x=if (x@num) x@x[i] else numeric(0L)
    )
  }
})

# column selection
#' @rdname tabMatrix-indexing
setMethod("[", c(x="tabMatrix", i="missing", j="index", drop="logical"), function (x, i, j, ..., drop=TRUE) {
  if (is.character(j)) stop("indexing by name currently not supported for tabMatrix")
  j <- as.integer(j)
  #if (any(j < 0L)) stop("negative integer indexing currently not supported for tabMatrix")
  if (all(j < 0L)) j <- seq_len(x@Dim[2L])[j]
  #if (any(j == 0L | j > x@Dim[2L])) stop("out of range")
  if (any(j <= 0L | j > x@Dim[2L])) stop("out of range")
  if (any(duplicated(j))) stop("duplicated column indices currently not supported for tabMatrix")
  if (drop && length(j) == 1L) {
    out <- rep.int(0, x@Dim[1L])
    ind <- which(x@perm == j - 1L)
    if (length(ind)) out[ind] <- if (x@num) x@x[ind] else 1
    out
  } else {
    reduced <- x@reduced
    perm <- match(x@perm, j - 1L)
    if (anyNA(perm)) {
      perm[is.na(perm)] <- 0L
      reduced <- TRUE
    }
    tabMatrix(Dim=c(x@Dim[1L], length(j)), Dimnames=list(x@Dimnames[[1L]], x@Dimnames[[2L]][j]),
      reduced=reduced, perm=perm - 1L, num=x@num, x=if (x@num) x@x else numeric(0L)
    )
  }
})

setMethod("nnzero", "tabMatrix", function(x, na.counted=NA) {
  if (x@reduced) {
    if (x@num)
      sum(x@x != 0 & x@perm != -1L)
    else
      sum(x@perm != -1L)
  } else {
    if (x@num)
      sum(x@x != 0)
    else
      x@Dim[1L]
  }
})

#' S4 method for generic 'anyNA' and signature 'tabMatrix'
#'
#' @keywords internal
#' @param x a tabMatrix object.
#' @param recursive not used.
#' @return whether the tabMatrix object contains missings or not.
setMethod("anyNA", "tabMatrix", function(x, recursive=FALSE)
  if (x@num) anyNA(x@x) else FALSE
)

setMethod("show", "tabMatrix", function(object) str(object))

# expand a vector, using the indicator part of tabM
expand_mv <- function(tabM, v) Ctab_numeric_prod(tabM, v, TRUE)

remove_levels <- function(f, l=1L) {
  lvs <- attr(f, "levels")
  if (length(l) != 1L || l < 1L || l > length(lvs)) stop("incorrect input")
  lvsnew <- lvs[-l]
  structure(match(lvs, lvsnew)[f], levels=lvsnew, class="factor")
}

# factor (interaction) to tabMatrix
# fvars looked up in data, which can either be a data.frame or an environment
# contrasts: either "contr.treatment", "contr.SAS" or an integer vector of the same length as
#   fvars specifying for each factor variable which level (by name) is considered the baseline.
#   If left unspecified no levels are removed.
# drop.unused.levels: drop unused levels of each factor
# drop: drop empty cells in the interaction (which may result in more columns to be dropped)
fac2tabM <- function(fvars, data, enclos=.GlobalEnv, x=numeric(), xlab=character(), drop.unused.levels=FALSE, drop=FALSE, contrasts=NULL, varsep=":", catsep="$", lex.order=FALSE) {
  if (missing(fvars) || !length(fvars)) stop("unexpected 'fvars' argument")
  for (f in seq_along(fvars)) {
    fac <- eval_in(fvars[f], data, enclos)
    if (!is.factor(fac) || drop.unused.levels) fac <- factor(fac)
    if (!is.null(contrasts)) {
      # TODO instead of changing levels here, compute the indices of the final codes vector
      #      which should be set to -1; this would be a lot faster for long factors
      if (length(contrasts) == 1L && contrasts == "contr.treatment") {  # R default: first level is baseline
        fac <- remove_levels(fac, 1L)
      } else {
        if (length(contrasts) == 1L && contrasts == "contr.SAS") {  # 'SAS' convention: last level is baseline
          fac <- remove_levels(fac, length(attr(fac, "levels")))
        } else {  # user-specified base
          if (fvars[f] %in% names(contrasts)) {
            m <- match(contrasts[fvars[f]], attr(fac, "levels"))
            if (length(m) != 1L || is.na(m)) stop("invalid contrasts")
            fac <- remove_levels(fac, m)
          }
        }
      }
    }
    if (f == 1L) {
      out <- fac
      attr(out, "levels") <- paste(fvars[1L], attr(out, "levels"), sep=catsep)
    } else {
      out <- interaction(out, fac, drop=drop, sep=paste0(varsep, fvars[f], catsep), lex.order=lex.order)
    }
  }
  if (length(xlab))
    labs <- paste(attr(out, "levels"), xlab, sep=":")
  else
    labs <- attr(out, "levels")
  out <- as.integer(out) - 1L
  if (!is.null(contrasts)) out[is.na(out)] <- -1L
  tabMatrix(Dim=c(length(out), length(labs)), Dimnames=list(NULL, labs),
    reduced=!is.null(contrasts), perm=out, num=(length(x) > 0L), x=as.double(x),
    xlab=xlab
  )
}

# formula to list of tabMatrices
# ... passed to fac2tabM
tables2tabM <- function(formula, data, ...) {
  n <- nrow(data)
  trms <- terms(formula, data=data)
  tmat <- terms_matrix(trms)
  if (!length(tmat)) {
    if (intercept_only(formula))
      return(list(new("tabMatrix", perm=rep.int(0L, n), reduced=FALSE, num=TRUE, x=rep.int(1, n), Dim=c(n, 1L))))
    else
      stop("empty formula")
  }
  tnames <- colnames(tmat)
  vnames <- rownames(tmat)
  qvar <- !catvars(trms, data)  # quantitative variables
  qvar <- vnames[which(qvar)]
  out <- setNames(vector(mode="list", length(tnames)), tnames)
  for (k in seq_along(tnames)) {
    countvars <- intersect(vnames[tmat[, k] > 0L], qvar)
    if (length(countvars)) {
      xk <- eval_in(countvars[1L], data)
      if (is.matrix(xk)) stop("matrix variables not yet allowed")
      labk <- countvars[1L]
      for (v in countvars[-1L]) {
        temp <- eval_in(v, data)
        xk <- xk * temp
        labk <- paste(labk, v, sep=":")
      }
    } else {
      xk <- NULL
      labk <- NULL
    }
    facvars <- setdiff(vnames[tmat[, k] > 0L], qvar)
    if (length(facvars)) {
      if (length(countvars)) {
        fk <- fac2tabM(facvars, data, x=xk, xlab=labk, contrasts=NULL, ...)
      } else {
        fk <- fac2tabM(facvars, data, contrasts=NULL, ...)
      }
    } else {
      if (length(countvars)) {
        fk <- new("tabMatrix", perm=rep.int(0L, n), reduced=FALSE, num=TRUE, x=xk, xlab=labk,
          Dim=c(n, 1L), Dimnames=list(NULL, labk))
      } else {
        # this should not happen, as the intercept only case is handled above
        fk <- new("tabMatrix", perm=rep.int(0L, n), reduced=FALSE, num=TRUE, x=rep.int(1, n), Dim=c(n, 1L))
      }
    }
    out[[k]] <- fk
  }
  out
}
