
# only system A, Lt and L solves used so far
.solve.systems <- setNames(1:9, c("A", "LDLt", "LD", "DLt", "L", "Lt", "D", "P", "Pt"))

# define R function to set defaults
# ordering=1L (AMD only) seems to give same permutations as Matrix::Cholesky
# ordering=0L is the cholmod default and also tries nested dissection
Cholesky_dsC <- function(A, perm=TRUE, LDL=FALSE, super=NA, Imult=0, ordering=.opts$chol.ordering)
  cCHM_dsC_Cholesky(A, perm, LDL, super, Imult, ordering)


# in case of permutation/pivoting, the decomposition is PtLLtP
build_chol <- function(M, perm=NULL, LDL=FALSE, super=NA, Imult=0, ordering=.opts$chol.ordering) {

  if (LDL) stop("LDL decompositions not (yet) supported")
  type <- class(M)[1L]
  if (all(type != c("numeric", "matrix", "ddiMatrix", "dsCMatrix"))) stop("'", type, "' is not a supported matrix class")
  size <- if (type == "numeric") length(M) else dim(M)[1L]

  if (type == "dsCMatrix") {
    if (is.null(perm)) {
      # TODO find better decision rule for use of permutation in cholesky factorization
      perm <- size > 100L
    }
  } else {
    perm <- FALSE
  }

  # define update, solve and Ltimes methods with arguments
  # update
  #   parent: new matrix to be decomposed; in case of sparse matrix the same zero pattern is assumed
  #   mult: mult*I is added to parent before decomposition
  # solve
  #   rhs: right hand side, of type numeric, matrix or dgCMatrix
  #   system: the type of system to be solved, see Matrix package
  #   systemP: if TRUE modifies system Lt to LtP and L to PtL; has no effect if PERM=FALSE or system="A" or type="ddiMatrix"
  # Ltimes left-multiplies a vector/matrix by L (in case of permutations P'L which need not be lower triangular!) or its transpose
  switch(type,
    numeric=, ddiMatrix = {  # diagonal represented as vector case (including trivial scalar case)
      if (type == "ddiMatrix" && M@diag == "U") {  # trivial case, identity matrix, no updates
        if (Imult != 0) stop("unit ddiMatrix assumed fixed")
        cholM <- M
        update <- function(parent, mult=0) if (mult != 0) stop("unit ddiMatrix assumed fixed")
        Ltimes <- function(x, transpose=TRUE) x
        solve <- function(rhs, system="A", systemP=FALSE) rhs
      } else {
        if (type == "ddiMatrix") {
          cholM <- sqrt(M@x + Imult)  # a vector is sufficient as internal representation
          update <- function(parent, mult=0) cholM <<- sqrt(parent@x + mult)
        } else {
          cholM <- sqrt(M + Imult)
          update <- function(parent, mult=0) cholM <<- sqrt(parent + mult)
        }
        Ltimes <- function(x, transpose=TRUE) cholM * x
        solve <- function(rhs, system="A", systemP=FALSE) {
          switch(class(rhs)[1L],
            numeric=, matrix =
              switch(system,
                A = rhs / cholM^2,
                Lt=, L = rhs / cholM
              ),
            dgCMatrix =
              switch(system,
                A = {
                  attr(rhs, "x") <- rhs@x / cholM[rhs@i + 1L]^2
                  rhs
                },
                Lt=, L = {
                  attr(rhs, "x") <- rhs@x / cholM[rhs@i + 1L]
                  rhs
                }
              ),
            stop("'", class(rhs)[1L], "' not supported as right hand side in solve")
          )
        }
      }
    },
    dsCMatrix = {
      cholM <- Cholesky_dsC(M, perm=perm, LDL=LDL, super=super, Imult=Imult, ordering=ordering)
      perm <- is.unsorted(cholM@perm, strictly = TRUE)
      if (perm) {  # permutation matrix unaffected by updates
        P <- cholM@perm + 1L
        iP <- invPerm(cholM@perm, zero.p=TRUE)
      }
      if (.opts$chol.inplace) {
        # requires that zero-pattern does not change!
        update <- function(parent, mult=0) cCHM_update_inplace(cholM, parent, mult)
      } else {
        # for forked parallel processing cholM should probably be stored in the state list p
        update <- function(parent, mult=0) cholM <<- .updateCHMfactor(cholM, parent, mult)
      }
      if (perm) {
        # TODO check whether the permutation is handled correctly in the transpose=FALSE case
        Ltimes <- function(x, transpose=TRUE) {
          L <- as(cholM, "Matrix")
          if (transpose) {
            if (is.vector(x)) crossprod_mv(L, x[P]) else crossprod_mv(L, x[P, , drop=FALSE])
          } else {
            if (is.vector(x)) (L %m*v% x)[iP] else (L %m*v% x)[iP, , drop=FALSE]
          }
        }
        solve <- function(rhs, system="A", systemP=FALSE) {
          switch(class(rhs)[1L],
            numeric =
              if (systemP)
                switch(system,
                  Lt = cCHMf_solve(cholM, rhs, 6L)[iP],
                  L = cCHMf_solve(cholM, rhs[P], 5L)
                )
              else
                cCHMf_solve(cholM, rhs, .solve.systems[system]),
            matrix =
              if (systemP)
                switch(system,
                  Lt = cCHMf_solve_matrix(cholM, rhs, 6L)[iP, , drop=FALSE],
                  L = cCHMf_solve_matrix(cholM, rhs[P, , drop=FALSE], 5L)
                )
              else
                cCHMf_solve_matrix(cholM, rhs, .solve.systems[system]),
            dgCMatrix =
              if (systemP)
                switch(system,
                  Lt = cCHMf_spsolve(cholM, rhs, 6L)[iP, , drop=FALSE],
                  L = cCHMf_spsolve(cholM, rhs[P, , drop=FALSE], 5L)
                )
              else
                cCHMf_spsolve(cholM, rhs, .solve.systems[system]),
            stop("'", class(rhs)[1L], "' not supported as right hand side in solve")
          )
        }
      } else {
        Ltimes <- function(x, transpose=TRUE) {
          L <- as(cholM, "Matrix")
          if (is.vector(x)) {
            if (transpose) crossprod_mv(L, x) else L %m*v% x
          } else {
            if (transpose) crossprod_mm(L, x) else L %m*m% x
          }
        }
        solve <- function(rhs, system="A", systemP=FALSE)
          switch(class(rhs)[1L],
            numeric = cCHMf_solve(cholM, rhs, .solve.systems[system]),
            matrix = cCHMf_solve_matrix(cholM, rhs, .solve.systems[system]),
            dgCMatrix = cCHMf_spsolve(cholM, rhs, .solve.systems[system]),
            stop("'", class(rhs)[1L], "' not supported as right hand side in solve")
          )
      }
    },
    matrix = {  # LDL, super ignored
      # NB Ccholesky currently returns upper triangular matrix, i.e. Lt
      if (Imult == 0)
        cholM <- Ccholesky(M)
      else
        cholM <- Ccholesky(add_diagC(M, rep.int(Imult, size)))
      update <- function(parent, mult=0) {
        if (mult == 0)
          cholM <<- Ccholesky(parent)
        else
          cholM <<- Ccholesky(add_diagC(parent, rep.int(mult, size)))
      }
      Ltimes <- function(x, transpose=TRUE) {
        if (is.vector(x)) {
          if (transpose) cholM %m*v% x else crossprod_mv(cholM, x)
        } else {
          if (transpose) cholM %m*m% x else crossprod_mm(cholM, x)
        }
      }
      solve <- function(rhs, system="A", systemP=FALSE) {
        switch(class(rhs)[1L],
          numeric =
            switch(system,
              A = Cbacksolve(cholM, Cforwardsolve(cholM, rhs)),
              Lt = Cbacksolve(cholM, rhs),
              L = Cforwardsolve(cholM, rhs),
              stop("unsupported solve system")
            ),
          matrix =
            switch(system,
              A = CbacksolveM(cholM, CforwardsolveM(cholM, rhs)),
              Lt = CbacksolveM(cholM, rhs),
              L = CforwardsolveM(cholM, rhs),
              stop("unsupported solve system")
            ),
          # disallow dgCMatrix rhs with dense M; in this case it is better to force rhs to be dense as well, as the result is typically dense
          #dgCMatrix = {  # this works, but may be slow; TODO instead of dgC, return the result as a dense matrix
          #  nc <- dim(rhs)[2L]
          #  out <- rhs
          #  for (i in seq_len(nc)) out[, i] <- solve(rhs[, i], system)
          #  out
          #}
          stop("'", class(rhs)[1L], "' not supported as right hand side in solve")
        )
      }
    }
  )  # END switch(type, ...)

  rm(M, type, LDL, super, Imult)
  environment()
}
