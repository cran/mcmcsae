
#' Set up a a function for direct sampling from a constrained multivariate normal distribution
#'
#' @keywords internal
#' @param mbs block component containing several model components.
#' @param X design matrix.
#' @param Q (structured) precision matrix.
#' @param R equality restriction matrix.
#' @param r rhs vector for equality constraints \eqn{R'x = r}, where \eqn{R'} denotes the transpose of R.
#' @param sampler sampler object as created by \code{\link{create_sampler}}.
#' @param name name of the cMVN vector parameter.
#' @param chol.control options for Cholesky decomposition, see \code{\link{chol_control}}.
#' @returns An environment with precomputed quantities and functions for sampling
#'   from a multivariate normal distribution subject to equality constraints.
create_cMVN_sampler <- function(mbs, X, Q, R=NULL, r=NULL, sampler, name="x", chol.control) {

  if (name == "") stop("empty name")
  q <- ncol(Q)

  # TODO create a function for handling R, as this code is duplicated from TMVN_sampler
  if (is.null(R)) {
    # special case, mainly useful for testing purposes
    cholQ <- build_chol(Q, control=chol.control)
  } else {
    R <- economizeMatrix(R, vec.as.diag=FALSE, check=TRUE)
    if (nrow(R) != q || ncol(R) > q) stop("incompatible constraint matrix 'R'")
    if (is.null(r)) {
      r <- numeric(ncol(R))
    } else {
      r <- as.numeric(r)
      if (length(r) != ncol(R)) stop("length of 'r' should equal the number of columns of 'R'")
    }
    # TODO
    #  check that R has full column rank
    #  parameters eps1, eps2 to modify Q.expanded for more numerical stability
    #   add eps1*I to Q block and -eps2*I to zero block
    #   this seems to help improve the accuracy with which the constraints are satisfied
    #   especially after updates in draw
    Q.expanded <- economizeMatrix(
      rbind(
        cbind(Q, R),
        cbind(t(R), Matrix(0, ncol(R), ncol(R)))
      ), sparse=TRUE, symmetric=TRUE
    )  # support dsCMatrix only
    # seems that besides LDL=TRUE we need perm=FALSE for indefinite matrix??
    chol.control$perm <- FALSE
    cholQ <- build_chol(Q.expanded, control=chol.control, LDL=TRUE)  # indefinite system, need LDLt
  }

  
  # X, QT passed from block's draw function
  draw <- function(p, Xy, X, QT, Imult=0) {

    # Xy is rhs, i.e. X' Qn ytilde (+ possibly prior reg term if b0 != 0)
    if (is.null(p[["sigma_"]])) sigma <- 1 else sigma <- p[["sigma_"]]
    u <- Xy + sigma * crossprod_mv(X, sampler$drawMVNvarQ(p))
    for (mc in mbs) {
      if (mc[["type"]] == "gen")
        u[mc$block.i] <- u[mc$block.i] + sigma * mc$drawMVNvarQ(p)
      else
        u[mc$block.i] <- u[mc$block.i] + mc$drawMVNvarQ(p)
    }

    if (is.null(R)) {
      cholQ$update(QT)
    } else {
      attr(Q.expanded, "x")[seq_along(QT@x)] <- QT@x
      cholQ$update(Q.expanded, Imult)
    }
 
    p[[name]] <- cholQ$solve(c(u, r))[1:q]
    p
  }

  environment()
}
