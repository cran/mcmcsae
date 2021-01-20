
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
  if (!(type %in% c("numeric", "matrix", "ddiMatrix", "dsCMatrix"))) stop("build_chol: unsupported matrix class '", type, "'")
  size <- if (type == "numeric") length(M) else dim(M)[1L]

  if (is.null(perm)) {
    # TODO find better decision rule for use of permutation in cholesky factorization
    perm <- size > 100L
    # perm can still change below depending on type
  }

  # define update, solve and crossprodL methods with arguments
  # update
  #   parent: new matrix to be decomposed; in case of sparse matrix the same zero pattern is assumed
  #   mult: mult*I is added to parent before decomposition
  # solve
  #   rhs: right hand side, of type numeric, matrix or dgCMatrix
  #   system: the type of system to be solved, see Matrix package
  #   systemP: if TRUE modifies system Lt to LtP and L to PtL; has no effect if PERM=FALSE or system="A" or type="ddiMatrix"
  # crossprodL takes the crossprod of L (in case of permutations P'L which need not be lower triangular!) with rhs
  #   rhs: for now assumed to be a numeric vector
  switch(type,
    numeric=, ddiMatrix = {  # diagonal represented as vector case (including trivial scalar case)
      perm <- FALSE
      if (type == "ddiMatrix" && M@diag == "U") {  # trivial case, identity matrix, no updates
        if (Imult != 0) stop("unit ddiMatrix assumed fixed")
        cholM <- M
        update <- function(parent, mult=0) if (mult != 0) stop("unit ddiMatrix assumed fixed")
        crossprodL <- function(rhs) rhs
        solve <- function(rhs, system="A", systemP=FALSE) rhs
      } else {
        if (type == "ddiMatrix") {
          cholM <- sqrt(M@x + Imult)  # a vector is sufficient as internal representation
          update <- function(parent, mult=0) cholM <<- sqrt(parent@x + mult)
        } else {
          cholM <- sqrt(M + Imult)
          update <- function(parent, mult=0) cholM <<- sqrt(parent + mult)
        }
        crossprodL <- function(rhs) cholM * rhs
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
      perm <- !identical(cholM@perm, 0:(length(cholM@perm) - 1L))
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
        crossprodL <- function(rhs) crossprod_mv(as(cholM, "Matrix"), rhs[P])
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
        crossprodL <- function(rhs) crossprod_mv(as(cholM, "Matrix"), rhs)
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
      # NB chol stores upper triangular matix, i.e. Lt
      if (Imult == 0) {
        cholM <- chol.default(M, pivot=perm)
      } else {
        cholM <- chol.default(add_diagC(M, rep.int(Imult, size)), pivot=perm)
      }
      if (perm) {  # NB this assumes that permutation never changes....CHECK!!!
        P <- attr(cholM, "pivot")
        iP <- invPerm(P)
      }
      update <- function(parent, mult=0) {
        if (mult == 0)
          cholM <<- chol.default(parent, pivot=perm)
        else
          cholM <<- chol.default(add_diagC(parent, rep.int(mult, size)), pivot=perm)
        if (perm) {  # permutation may change; alternatively compute P and iP in solve function as needed
          P <<- attr(cholM, "pivot")
          iP <<- invPerm(P)
        }
      }
      if (perm) {
        crossprodL <- function(rhs) cholM %m*v% rhs[P]
        solve <- function(rhs, system="A", systemP=FALSE) {
          switch(class(rhs)[1L],
            numeric =
              switch(system,
                A = Cbacksolve(cholM, Cforwardsolve(cholM, rhs[P]))[iP],
                Lt = if (systemP) Cbacksolve(cholM, rhs)[iP] else Cbacksolve(cholM, rhs),
                L = if (systemP) Cforwardsolve(cholM, rhs[P]) else Cforwardsolve(cholM, rhs),
                stop("unsupported solve system")
              ),
            matrix =
              switch(system,
                A = CbacksolveM(cholM, CforwardsolveM(cholM, rhs[P, , drop=FALSE]))[iP, , drop=FALSE],
                Lt = if (systemP) CbacksolveM(cholM, rhs)[iP, , drop=FALSE] else CbacksolveM(cholM, rhs),
                L = if (systemP) CforwardsolveM(cholM, rhs[P, , drop=FALSE]) else CforwardsolveM(cholM, rhs),
                stop("unsupported solve system")
              ),
            stop("'", class(rhs)[1L], "' not supported as right hand side in solve")
          )
        }
      } else {
        crossprodL <- function(rhs) cholM %m*v% rhs
        solve <- function(rhs, system="A", systemP=FALSE) {
          switch(class(rhs)[1L],
            numeric =
              switch(system,
                A = Cbacksolve(cholM, Cforwardsolve(cholM, rhs)),
                Lt = Cbacksolve(cholM, rhs),
                L = Cforwardsolve(cholM, rhs),
                #LDLt = backsolve(cholM, forwardsolve(cholM, rhs, upper.tri=TRUE, transpose=TRUE)),
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
    }
  )  # END switch(type, ...

  rm(M, type, LDL, super, Imult)
  environment()
}
