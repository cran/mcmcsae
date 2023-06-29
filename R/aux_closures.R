
# Return a closure for summation of symmetric matrices of a fixed (zero-)structure.
# supported matrices are matrix, dsCMatrix, ddiMatrix
# at most one of M0, M2 can be null
# M0 is completely fixed
# M1 and M2 can change in value, but not in zero-structure
# function returned computes M0 + w1*M1 + w2*M2
# sparse: should the matrix sum template and output be sparse? for now only works when !is.null(M0)
make_mat_sum <- function(M0=NULL, M1, M2=NULL, sparse=NULL) {
  # TODO multiply M1 and M2 by a positive random number, to exclude the possibility of coincidental cancellations in the sum
  if (!is.null(M0)) {
    if (all(class(M0)[1L] != c("ddiMatrix", "dsCMatrix", "matrix"))) stop("unsupported matrix type")
    if (class(M0)[1L] == "dsCMatrix") M0 <- drop0(M0, is.Csparse=TRUE)
  }
  if (is.null(M2)) {
    if (is.null(M0)) stop("at least one of M0 and M2 must be non-NULL")
    template <- economizeMatrix(M0 + M1, sparse=sparse, symmetric=TRUE)
  } else {
    if (is.null(M0))
      template <- economizeMatrix(M1 + M2, symmetric=TRUE)
    else
      template <- economizeMatrix(M0 + M1 + M2, sparse=sparse, symmetric=TRUE)
  }
  if (is.null(M0)) {
    classes <- paste0(substr(class(M1)[1L], 1L, 3L), substr(class(M2)[1L], 1L, 3L))
    switch(classes,
      matmat = {
        rm(template)
        update <- function(M1, M2, w1=1, w2=1) w1 * M1 + w2 * M2
      },
      matddi = {
        rm(template)
        update <- function(M1, M2, w1=1, w2=1) add_diagC(w1 * M1, w2 * ddi_diag(M2))
      },
      ddimat = {
        rm(template)
        update <- function(M1, M2, w1=1, w2=1) add_diagC(w2 * M2, w1 * ddi_diag(M1))
      },
      matdsC = {
        IM <- cbind(M2@i + 1L, rep.int(seq_len(M2@Dim[2L]), diff(M2@p)))  # one-based index
        symm <- which(IM[, 1L] != IM[, 2L])  # rows that should be added in order to symmetrize the result
        rm(template)
        update <- function(M1, M2, w1=1, w2=1) {
          out <- w1 * M1
          out[IM] <- out[IM] + w2 * M2@x
          out[IM[symm, 2:1, drop=FALSE]] <- out[IM[symm, 2:1, drop=FALSE]] + w2 * M2@x[symm]  # for symmetry
          out
        }
      },
      dsCmat = {
        IM <- cbind(M1@i + 1L, rep.int(seq_len(M1@Dim[2L]), diff(M1@p)))
        symm <- which(IM[, 1L] != IM[, 2L])
        rm(template)
        update <- function(M1, M2, w1=1, w2=1) {
          out <- w2 * M2
          out[IM] <- out[IM] + w1 * M1@x
          out[IM[symm, 2:1, drop=FALSE]] <- out[IM[symm, 2:1, drop=FALSE]] + w1 * M1@x[symm]
          out
        }
      },
      ddiddi = {
        update <- function(M1, M2, w1=1, w2=1) {
          out <- template
          attr(out, "x") <- w1 * ddi_diag(M1) + w2 * ddi_diag(M2)
          out
        }
      },
      ddidsC=, dsCddi=, dsCdsC = {
        UnitDiag1 <- UnitDiag2 <- FALSE
        if (isUnitDiag(M1)) {
          M1 <- expandUnitDiag(M1)
          UnitDiag1 <- TRUE
        }
        if (isUnitDiag(M2)) {
          M2 <- expandUnitDiag(M2)
          UnitDiag2 <- TRUE
        }
        attr(M1, "x") <- rep.int(1, length(M1@x))
        attr(M2, "x") <- rep.int(2, length(M2@x))
        template <- economizeMatrix(M1 + M2, symmetric=TRUE, sparse=TRUE)
        ind1 <- which(template@x %in% c(1,3)) - 1L  # zero-based index as used in C++ code
        ind2 <- which(template@x %in% c(2,3)) - 1L
        nx <- length(template@x)
        update <- function(M1, M2, w1=1, w2=1) {
          out <- template
          attr(out, "x") <- sparse_sum_x(nx, ind1, ind2, M1@x, M2@x, UnitDiag1, UnitDiag2, w1, w2)
          out
        }
      },
      {
        warn("possibly inefficient (sparse) matrix addition")
        update <- function(M1, M2, w1=1, w2=1) w1 * M1 + w2 * M2
      }
    )
  } else {  # including fixed contribution M0
    classes <- paste0(substring(class(template)[1L], 1L, 3L), substring(class(M1)[1L], 1L, 3L), substring(class(M2)[1L], 1L, 3L))
    switch(classes,
      # TODO always represent diagonal matrices by numeric vectors instead of matrix or ddiMatrix
      matmatNUL=, ddimatNUL = {
        template <- as.matrix(M0)
        update <- function(M1, M2, w1=1, w2=1) template + w1 * M1
      },
      matddiNUL = {
        template <- as.matrix(M0)
        update <- function(M1, M2, w1=1, w2=1) add_diagC(template, w1 * ddi_diag(M1))
      },
      matdsCNUL = {
        IM <- cbind(M1@i + 1L, rep.int(seq_len(M1@Dim[2L]), diff(M1@p)))
        symm <- which(IM[, 1L] != IM[, 2L])
        template <- as.matrix(M0)
        update <- function(M1, M2, w1=1, w2=1) {
          out <- template
          out[IM] <- out[IM] + w1 * M1@x
          out[IM[symm, 2:1, drop=FALSE]] <- out[IM[symm, 2:1, drop=FALSE]] + w1 * M1@x[symm]
          out
        }
      },
      dsCmatNUL = {
        template <- as.matrix(M0)
        update <- function(M1, M2, w1=1, w2=1) template + w1 * M1
      },
      dsCdsCNUL=, dsCddiNUL = {
        UnitDiag1 <- FALSE
        if (isUnitDiag(M1)) {
          M1 <- expandUnitDiag(M1)
          UnitDiag1 <- TRUE
        }
        if (isUnitDiag(M0)) M0 <- expandUnitDiag(M0)
        M0x <- M0@x
        attr(M0, "x") <- rep.int(1, length(M0@x))
        attr(M1, "x") <- rep.int(2, length(M1@x))
        template <- economizeMatrix(M0 + M1, symmetric=TRUE, sparse=TRUE)
        ind1 <- which(template@x %in% c(2,3))
        ind2 <- integer(0L); M2x <- numeric(0L)
        attr(template, "x")[ind1] <- template@x[ind1] - 2
        attr(template, "x")[template@x == 1] <- M0x
        ind1 <- ind1 - 1L  # zero-based index as used in C++ code
        # now template equals M0, but possibly with redundant zeros to accommodate M1's contribution
        nx <- length(template@x)
        rm(M0x)
        update <- function(M1, M2, w1=1, w2=1) {
          out <- template
          attr(out, "x") <- out@x + sparse_sum_x(nx, ind1, ind2, M1@x, M2x, UnitDiag1, FALSE, w1, w2)
          out
        }
      },
      ddiddiNUL=, ddidsCNUL = {  # the 2nd case dsC is actually diagonal
        template <- as(M0, "diagonalMatrix")
        UnitDiag1 <- isUnitDiag(M1)
        if (isUnitDiag(template)) template <- expandUnitDiag(template)
        update <- function(M1, M2, w1=1, w2=1) {
          out <- template
          attr(out, "x") <- out@x + w1 * if (UnitDiag1) 1 else M1@x
          out
        }
      },
      matmatmat = {
        template <- as.matrix(M0)
        update <- function(M1, M2, w1=1, w2=1) template + w1 * M1 + w2 * M2
      },
      matmatddi = {
        template <- as.matrix(M0)
        update <- function(M1, M2, w1=1, w2=1) add_diagC(template + w1 * M1, w2 * ddi_diag(M2))
      },
      matddimat = {
        template <- as.matrix(M0)
        update <- function(M1, M2, w1=1, w2=1) add_diagC(template + w2 * M2, w1 * ddi_diag(M1))
      },
      matddiddi = {
        template <- as.matrix(M0)
        update <- function(M1, M2, w1=1, w2=1) add_diagC(template, w1 * ddi_diag(M1) + w2 * ddi_diag(M2))
      },
      matmatdsC = {
        IM <- cbind(M2@i + 1L, rep.int(seq_len(M2@Dim[2L]), diff(M2@p)))
        symm <- which(IM[, 1L] != IM[, 2L])
        template <- as.matrix(M0)
        update <- function(M1, M2, w1=1, w2=1) {
          out <- template + w1 * M1
          out[IM] <- out[IM] + w2 * M2@x
          out[IM[symm, 2:1, drop=FALSE]] <- out[IM[symm, 2:1, drop=FALSE]] + w2 * M2@x[symm]
          out
        }
      },
      matdsCmat = {
        IM <- cbind(M1@i + 1L, rep.int(seq_len(M1@Dim[2L]), diff(M1@p)))
        symm <- which(IM[, 1L] != IM[, 2L])
        template <- as.matrix(M0)
        update <- function(M1, M2, w1=1, w2=1) {
          out <- template + w2 * M2
          out[IM] <- out[IM] + w1 * M1@x
          out[IM[symm, 2:1, drop=FALSE]] <- out[IM[symm, 2:1, drop=FALSE]] + w1 * M1@x[symm]
          out
        }
      },
      matdsCddi = {
        IM <- cbind(M1@i + 1L, rep.int(seq_len(M1@Dim[2L]), diff(M1@p)))
        symm <- which(IM[, 1L] != IM[, 2L])
        template <- as.matrix(M0)
        update <- function(M1, M2, w1=1, w2=1) {
          out <- template
          out[IM] <- out[IM] + w1 * M1@x
          out[IM[symm, 2:1, drop=FALSE]] <- out[IM[symm, 2:1, drop=FALSE]] + w1 * M1@x[symm]
          add_diagC(out, w2 * ddi_diag(M2))
        }
      },
      matddidsC = {
        IM <- cbind(M2@i + 1L, rep.int(seq_len(M2@Dim[2L]), diff(M2@p)))
        symm <- which(IM[, 1L] != IM[, 2L])
        template <- as.matrix(M0)
        update <- function(M1, M2, w1=1, w2=1) {
          out <- add_diagC(template, w1 * ddi_diag(M1))
          out[IM] <- out[IM] + w2 * M2@x
          out[IM[symm, 2:1, drop=FALSE]] <- out[IM[symm, 2:1, drop=FALSE]] + w2 * M2@x[symm]
          out
        }
      },
      matdsCdsC = {
        IM1 <- cbind(M1@i + 1L, rep.int(seq_len(M1@Dim[2L]), diff(M1@p)))
        symm1 <- which(IM1[, 1L] != IM1[, 2L])
        IM2 <- cbind(M2@i + 1L, rep.int(seq_len(M2@Dim[2L]), diff(M2@p)))
        symm2 <- which(IM2[, 1L] != IM2[, 2L])
        template <- as.matrix(M0)
        update <- function(M1, M2, w1=1, w2=1) {
          out <- template
          out[IM1] <- out[IM1] + w1 * M1@x
          out[IM1[symm1, 2:1, drop=FALSE]] <- out[IM1[symm1, 2:1, drop=FALSE]] + w1 * M1@x[symm1]
          out[IM2] <- out[IM2] + w2 * M2@x
          out[IM2[symm2, 2:1, drop=FALSE]] <- out[IM2[symm2, 2:1, drop=FALSE]] + w2 * M2@x[symm2]
          out
        }
      },
      dsCddiddi=, dsCdsCddi=, dsCddidsC=, dsCdsCdsC = {
        UnitDiag1 <- UnitDiag2 <- FALSE
        if (isUnitDiag(M1)) {
          M1 <- expandUnitDiag(M1)
          UnitDiag1 <- TRUE
        }
        if (isUnitDiag(M2)) {
          M2 <- expandUnitDiag(M2)
          UnitDiag2 <- TRUE
        }
        if (isUnitDiag(M0)) M0 <- expandUnitDiag(M0)
        M0x <- M0@x
        attr(M0, "x") <- rep.int(1, length(M0@x))
        attr(M1, "x") <- rep.int(2, length(M1@x))
        attr(M2, "x") <- rep.int(4, length(M2@x))
        template <- economizeMatrix(M0 + M1 + M2, symmetric=TRUE, sparse=TRUE)
        ind1 <- which(template@x %in% c(2,3,6,7))
        ind2 <- which(template@x %in% c(4,5,6,7))
        attr(template, "x")[ind1] <- template@x[ind1] - 2
        attr(template, "x")[ind2] <- template@x[ind2] - 4
        attr(template, "x")[template@x == 1] <- M0x
        ind1 <- ind1 - 1L  # zero-based index as used in C++ code
        ind2 <- ind2 - 1L
        # now template equals M0, but possibly with redundant zeros to accommodate M1's contribution
        nx <- length(template@x)
        rm(M0x)
        update <- function(M1, M2, w1=1, w2=1) {
          out <- template
          attr(out, "x") <- out@x + sparse_sum_x(nx, ind1, ind2, M1@x, M2@x, UnitDiag1, UnitDiag2, w1, w2)
          out
        }
      },
      ddiddiddi = {
        template <- as(M0, "diagonalMatrix")
        if (isUnitDiag(template)) template <- expandUnitDiag(template)
        update <- function(M1, M2, w1=1, w2=1) {
          out <- template
          attr(out, "x") <- out@x + w1 * ddi_diag(M1) + w2 * ddi_diag(M2)
          out
        }
      },
      {
        warn("possibly inefficient (sparse) matrix addition")
        template <- M0
        if (is.null(M2))
          update <- function(M1, M2, w1=1, w2=1) template + w1 * M1
        else
          update <- function(M1, M2, w1=1, w2=1) template + w1 * M1 + w2 * M2
      }
    )
  }
  rm(M0, M1, M2, classes)
  update
}


make_det <- function(M, perm=FALSE) {
  if (nrow(M) <= 1000L) {
    # for moderate dimensions eigenvalues are the fastest way for simple determinant updates
    ev0 <- eigen(M, only.values=TRUE)$values
    rm(M, perm)
    function(w1, w2) sum(log(w1 * ev0 + w2))
  } else {
    detchol <- build_chol(M, perm)
    d <- nrow(M)
    rm(perm)
    function(w1, w2) d * log(w1) + 2 * c(determinant(detchol$update(M, w2/w1)$cholM)$modulus)
  }
}
