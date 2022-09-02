#' Solve Ax=b by preconditioned conjugate gradients
#'
#' @export
#' @keywords internal
#' @param b right hand side vector.
#' @param env environment containing at least a function \code{A_times} that computes
#'  the matrix-vector product Ax for some input vector x, and a function \code{M_solve}
#'  that computes M^-1 x for some preconditioner matrix M.
#' @param x start value for the conjugate gradient algorithm.
#' @param max.it maximum number of iterations.
#' @param e total squared error stop criterion.
#' @param verbose whether progress information is shown.
#' @param ... any parameters passed to \code{A_times} and \code{M_solve}.
#' @return The (approximated) solution to Ax=b.
#' @references
#'  M.R. Hestenes and E. Stiefel (1952).
#'    Methods of conjugate gradients for solving linear systems.
#'    Journal of Research of the National Bureau of Standards 49(6), 409-436.
# TODO
#   - matrix rhs
#   - maybe in Rcpp
CG <- function(b, env, x=NULL, max.it=NULL, e=NULL, verbose=FALSE, ...) {
  if (is.null(x)) x <- 0*b
  if (is.null(max.it)) max.it <- length(b)
  if (is.null(e)) e <- length(b)*1e-12
  r <- b - env$A_times(x, ...)
  z <- env$M_solve(r, ...)
  p <- z
  k <- 0L
  repeat {
    Ap <- env$A_times(p, ...)
    rz <- dotprodC(r, z)
    alpha <-  rz / dotprodC(p, Ap)
    x <- x + alpha * p
    # if alpha*p sufficiently small exit; TODO better criterion
    #if (max(abs(alpha * p)) < 1e-6 || k > max.it) break; this may prevent alpha becoming NA (for small e)?
    r <- r - alpha * Ap
    # if r sufficiently small exit; THIS criterion does not seem to work for non-psd preconditioner
    #if (mean(abs(r)) < sqrt(.Machine$double.eps)) break
    z <- env$M_solve(r, ...)
    if (dotprodC(z, z) < e || k > max.it) {
      # criterion used in Nishimura and Suchard
      if (k > max.it) warning("max.it iterations reached with discrepancy ", dotprodC(z, z), immediate.=TRUE)
      break
    }
    #if (verbose) cat("iteration ", k, "   |r| = ", sqrt(sum(r^2)), "   |z| = ", sqrt(sum(z^2)), "\n")
    beta <- dotprodC(r, z) / rz
    p <- z + beta * p
    k <- k + 1L
  }
  x
}


#' Set up conjugate gradient sampler
#'
# for N(Q^-1 X' Q0 y, Q^-1)
# where Q = X' Q0 X + oplus_k QA_k x QV_k is q x q and (for mc_gen) QA_k = DA_k' DA_k may be singular
# NB experimental feature for now, needs more testing
#' @export
#' @keywords internal
#' @param mbs block component containing several model components.
#' @param X design matrix.
#' @param sampler sampler object as created by \code{\link{create_sampler}}.
#' @param max.it passed to \code{\link{CG}}.
#' @param stop.criterion passed to \code{\link{CG}}.
#' @param verbose passed to \code{\link{CG}}.
#' @param preconditioner one of  "GMRF", "GMRF2", "GMRF3" and "identity".
#' @param scale scale parameter for GMRF3 preconditioner.
#' @return An environment with precomputed quantities and functions for multiplication
#'   by A and by an (inverse) preconditioning matrix.
#' @references
#'  A. Nishimura and M.A. Suchard (2018).
#'    Prior-preconditioned conjugate gradient method for accelerated Gibbs sampling in
#'    large n & large p sparse Bayesian regression.
#'    arXiv:1810.12437.
# TODO
# - input checks
# - allow updates of QA (priorA or Leroux), i.p. QA = DA' diag(w) DA where w is updated
# - allow updates of X (model component mec)
# - use X's of underlying model components
#   and for each component use instead of X0 and XA instead of X, using mixed product relations
#   for Khatri-Rao etc.; X0 will typically be dense and XA tabMatrix
# - p >> n case
setup_CG_sampler <- function(mbs, X, sampler, max.it=NULL, stop.criterion=NULL, verbose=NULL,
                             preconditioner=c("GMRF", "GMRF2", "GMRF3", "identity"),
                             scale) {
  n <- nrow(X)
  q <- ncol(X)
  if (is.null(max.it)) max.it <- q
  if (is.null(stop.criterion)) stop.criterion <- 1e-12 * q
  if (is.null(verbose)) verbose <- FALSE
  preconditioner <- match.arg(preconditioner)

  if (sampler$family$family == "gaussian") {
    if (sampler$modeled.Q) {
      Q_x <- switch(sampler$Q0.type,
        unit=, diag = function(p, x) p[["Q_"]] * x,
        symm = function(p, x) p[["QM_"]] %m*v% x
      )
      cholQ <- switch(sampler$Q0.type,
        unit=, diag = build_chol(runif(n, 0.9, 1.1)),
        symm = build_chol(crossprod_sym(Cdiag(runif(0.9, 1.1)), sampler$Q0))
      )
    } else {
      Q_x <- switch(sampler$Q0.type,
        unit = function(p, x) x,
        diag = function(p, x) sampler$Q0@x * x,
        symm = function(p, x) sampler$Q0 %m*v% x
      )
      cholQ <- switch(sampler$Q0.type,
        unit = build_chol(Diagonal(n)),
        diag = build_chol(Cdiag(sampler$Q0@x)),
        symm = build_chol(sampler$Q0)
      )
    }
  } else {
    if (sampler$family$link == "probit") {
      Q_x <- function(p, x) x
      cholQ <- build_chol(Diagonal(n))
    } else {
      Q_x <- function(p, x) p[["Q_"]] * x
      cholQ <- build_chol(runif(n, 0.9, 1.1))
    }
  }

  # set up / precompute Cholesky factors
  # for gen these factors are updated in each MCMC iteration
  cholQV <- list()
  for (mc in mbs) {
    if (mc$type == "gen") {
      if (mc$var == "unstructured")
        cholQV[[mc$name]] <- build_chol(rWishart(1L, mc$q0, diag(mc$q0))[,,1L])
      else
        cholQV[[mc$name]] <- build_chol(runif(mc$q0, 0.5, 1.5))
    } else {
      if (mc$type != "reg") stop("TBI")
      # TODO account for b0, and optimize for zero Q0 (default)
      cholQV[[mc$name]] <- build_chol(mc$Q0)
    }
  }

  # QT = oplus_k QA x QV is created in parent's draw function
  # TODO check why sigma^2 factor in 2nd term should not be there
  if (sampler$sigma.fixed)
    A_times <- function(x, X, QT, cholQV, p)
      crossprod_mv(X, Q_x(p, X %m*v% x)) + QT %m*v% x
  else
    A_times <- function(x, X, QT, cholQV, p)
      crossprod_mv(X, Q_x(p, X %m*v% x)) + (QT %m*v% x) * p[["sigma_"]]^2

  # add block index vectors to mbs components
  n_vec_list <- 0L
  for (mc in mbs) {
    mc$i.bl <- (n_vec_list + 1L):(n_vec_list + mc[["q"]])
    n_vec_list <- n_vec_list + mc[["q"]]
    
    if (mc$type == "gen") {
      # TODO handle DA in mc_gen
      if (is.null(mc$DA)) {
        mc$DA <- mc$QA
        if (class(mc$DA)[1L] != "ddiMatrix" || mc$DA@diag != "U") stop("cannot derive DA")
      }
      # TODO
      # - duplicate code, see setup_priorGMRFsampler in mc_gen
      # - cholDD needs update in case of local shrinkage priorA
      mc$cholDD <- build_chol(tcrossprod(mc$DA), perm=FALSE)  # sometimes perm=TRUE may be more efficient
    } else {
      # regression usually uses non- or weakly-informative prior -->
      # use simple regression posterior variances in preconditioner, as suggested in NS paper
      # TODO more efficient computation of diagonal values only
      mc$gamma <- 2*diag(solve(crossprod_sym(mc$X, sampler$Q0)))
    }
  }
  rm(n_vec_list)

  switch(preconditioner,
    GMRF =
      # multiply by D+ DA x V, the (pseudo-)inverse of preconditioner M
      M_solve <- function(x, X, QT, cholQV, p) {
        out <- vector("numeric", q)
        for (mc in mbs) {
          if (mc$type == "gen") {
            temp <- mc$DA %m*m% t.default(cholQV[[mc$name]]$solve(matrix(x[mc$i.bl], nrow=mc$q0)))
            out[mc$i.bl] <- as.numeric(t.default(crossprod_mm(mc$DA, mc$cholDD$solve(temp))))
          } else {
            out[mc$i.bl] <- mc$gamma * x[mc$i.bl]
          }
        }
        out
      },
    GMRF2 =
      # multiply by D+ D+' x QV^-1, cf. Nishimura and Suchard
      M_solve <- function(x, X, QT, cholQV, p) {
        out <- vector("numeric", q)
        for (mc in mbs) {
          if (mc$type == "gen") {
            temp <- mc$DA %m*m% t.default(cholQV[[mc$name]]$solve(matrix(x[mc$i.bl], nrow=mc$q0)))
            temp <- mc$cholDD$solve(temp)
            out[mc$i.bl] <- as.numeric(t.default(crossprod_mm(mc$DA, mc$cholDD$solve(temp))))
          } else {
            out[mc$i.bl] <- mc$gamma * x[mc$i.bl]
          }
        }
        out
      },
    GMRF3 =
      # multiply by D+ D+' x QV^-1, cf. Nishimura and Suchard + scaling
      M_solve <- function(x, X, QT, cholQV, p) {
        out <- vector("numeric", q)
        for (mc in mbs) {
          if (mc$type == "gen") {
            temp <- mc$DA %m*m% t.default(cholQV[[mc$name]]$solve(matrix(scale * x[mc$i.bl], nrow=mc$q0)))
            temp <- mc$cholDD$solve(temp)
            out[mc$i.bl] <- as.numeric(t.default(crossprod_mm(mc$DA, mc$cholDD$solve(temp))))
          } else {
            out[mc$i.bl] <- mc$gamma * x[mc$i.bl]
          }
        }
        out
      },
    identity =
      M_solve <- function(x, X, QT, cholQV, p) x
  )
  self <- environment()
  # X, QT passed from block's draw function
  draw <- function(p, Xy, X, QT, sampler, start=NULL) {
    # Xy is rhs, i.e. X' Qn ytilde (+ possibly prior reg term if b0 != 0)
    y1 <- Crnorm(n)
    if (sampler$modeled.Q) {
      if (sampler$Q0.type == "symm")
        cholQ$update(p[["QM_"]])
      else
        cholQ$update(p[["Q_"]])
    }
    if (is.null(p[["sigma_"]])) sigma <- 1 else sigma <- p[["sigma_"]]
    u <- Xy + sigma * crossprod_mv(X, cholQ$Ltimes(y1, transpose=FALSE))
    for (mc in mbs) {
      if (mc$type == "gen") {
        y2 <- Crnorm(mc$q0 * nrow(mc$DA))
        dim(y2) <- c(mc$q0, nrow(mc$DA))
        cholQV[[mc$name]]$update(sigma^2 * p[[mc$name_Qv]])
        u[mc$i.bl] <- u[mc$i.bl] + as.numeric(cholQV[[mc$name]]$Ltimes(y2, transpose=FALSE) %m*m% mc$DA)
      } else {
        if (mc$informative.prior) {
          y2 <- Crnorm(mc$q)
          u[mc$i.bl] <- u[mc$i.bl] + cholQV[[mc$name]]$Ltimes(y2, transpose=FALSE)
        }
      }
    }
    # solve Qtot v = u using preconditioned conjugate gradients, where Qtot = X'QX + sigma^2 QT
    out <- CG(u, self, start, max.it=max.it, e=stop.criterion, verbose=verbose, X=X, QT=QT, cholQV=cholQV, p=p)
    if (verbose) cat("discrepancies: ", summary(A_times(out, X, QT, cholQV, p) - u), "\n")
    out
  }  # END function draw

  self
}
