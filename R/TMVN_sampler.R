#' Set up a sampler object for sampling from a possibly truncated and degenerate multivariate normal distribution
#'
#' This function sets up an object for multivariate normal sampling based on a specified precision matrix
#' or cholesky decomposition thereof. Linear equality and inequality restrictions are supported.
#' For sampling under inequality restrictions three algorithms are available. The default in that case is
#' an exact Hamiltonian Monte Carlo algorithm (Pakman and Paninski, 2014). Alternatively, a Gibbs sampling
#' algorithm can be used (Rodriguez-Yam et al., 2004). The third option is a data augmentation method
#' to sample from a smooth approximation to the truncated multivariate normal distribution (Souris et al., 2018).
#'
#' The componentwise Gibbs sampler uses univariate truncated normal samplers as described
#' in Botev and L'Ecuyer (2016). These samplers are implemented in R package \pkg{TruncatedNormal},
#' but here translated to C++ for an additional speed-up.
#'
#' @examples
#' \donttest{
#' S <- cbind(diag(2), c(-1, 1), c(1.1, -1))  # inequality matrix
#' # S'x >= 0 represents the wedge x1 <= x2 <= 1.1 x1
#' # example taken from Pakman and Paninski (2014)
#' # 1. exact Hamiltonian Monte Carlo (Pakman and Paninski, 2014)
#' sampler <- create_TMVN_sampler(Q=diag(2), mu=c(4, 4), S=S, method="HMC")
#' sim <- MCMCsim(sampler)
#' summary(sim)
#' plot(as.matrix(sim$x), pch=".")
#' # 2. Gibbs sampling approach (Rodriguez-Yam et al., 2004)
#' sampler <- create_TMVN_sampler(Q=diag(2), mu=c(4, 4), S=S, method="Gibbs")
#' sim <- MCMCsim(sampler)
#' summary(sim)
#' plot(as.matrix(sim$x), pch=".")
#' # 3. soft TMVN approximation (Souris et al., 2018)
#' sampler <- create_TMVN_sampler(Q=diag(2), mu=c(4, 4), S=S, method="softTMVN")
#' sim <- MCMCsim(sampler)
#' summary(sim)
#' plot(as.matrix(sim$x), pch=".")
#' }
#'
#' @export
#' @author Harm Jan Boonstra and Grzegorz Baltissen, who helped implement a first version
#'  of the TMVNV Gibbs sampler.
#' @param Q precision matrix (unconstrained) (n x n); used in case \code{cholQ} is not supplied.
#' @param perm whether permutation/pivoting should be used to build a Cholesky object.
#' @param mu mean of the (unconstrained) multivariate normal distribution.
#' @param Xy alternative to specifying mu; in this case \code{mu} is computed as \code{Q^{-1} Xy}.
#' @param update.Q whether \code{Q} is updated for each draw.
#' @param update.mu whether \code{mu} is updated for each draw. By default equal to \code{update.Q}.
#' @param name name of the TMVN vector parameter.
#' @param coef.names optional labels for the components of the vector parameter.
#' @param R equality restriction matrix.
#' @param r rhs vector for equality constraints \code{R^T x = r}.
#' @param S inequality restriction matrix.
#' @param s rhs vector for inequality constraints \code{S^T x >= s}.
#' @param lower alternative to \code{s} for two-sided inequality restrictions \code{lower <= S^T x <= upper}.
#' @param upper alternative to \code{s} for two-sided inequality restrictions \code{lower <= S^T x <= upper}.
#' @param method sampling method. The options are "direct" for direct sampling from the
#'  unconstrained or equality constrained multivariate normal (MVN). For inequality constrained
#'  MVN sampling three methods are supported: "HMC" for (exact) Hamiltonian Monte Carlo, "Gibbs" for a
#'  component-wise Gibbs sampling approach, and "softTMVN" for a data augmentation method that samples
#'  from a smooth approximation to the truncated MVN.
#' @param reduce whether to a priori restrict the simulation to the subspace defined by the
#'  equality constraints.
#' @param T.HMC the duration of a Hamiltonian Monte Carlo simulated particle trajectory. If a vector of
#'  length 2 is supplied it is interpreted as an interval from which the duration is drawn uniformly,
#'  independently in each HMC iteration.
#' @param print.info whether information about violations of inequalities and bounces off inequality walls
#'  is printed to the screen. This sometimes provides useful diagnostic information.
#' @param sharpness for method 'softTMVN', the sharpness of the soft inequalities; the larger the better
#'  the approximation of exact inequalities. It must be either a scalar value or a vector of length
#'  equal to the number of inequality restrictions.
#' @param useV for method 'softTMVN' whether to base computations on variance instead of precision
#'  matrices.
#' @return An environment for sampling from a possibly degenerate and truncated multivariate normal
#'  distribution.
#' @references
#'  Z.I. Botev and P. L'Ecuyer (2016).
#'    Simulation from the Normal Distribution Truncated to an Interval in the Tail.
#'    in VALUETOOLS.
#'
#'  Y. Cong, B. Chen and M. Zhou (2017).
#'    Fast simulation of hyperplane-truncated multivariate normal distributions.
#'    Bayesian Analysis 12(4), 1017-1037.
#'
#'  A. Pakman and L. Paninski (2014).
#'    Exact Hamiltonian Monte Carlo for truncated multivariate gaussians.
#'    Journal of Computational and Graphical Statistics 23(2), 518-542.
#'
#'  G. Rodriguez-Yam, R.A. Davis and L.L. Scharf (2004).
#'    Efficient Gibbs sampling of truncated multivariate normal with application to constrained linear regression.
#'    Unpublished manuscript.
#'
#'  H. Rue and L. Held (2005).
#'    Gaussian Markov Random Fields.
#'    Chapman & Hall/CRC.
#'
#'  A. Souris, A. Bhattacharya and P. Debdeep (2018).
#'    The Soft Multivariate Truncated Normal Distribution.
#'    arXiv:1807.09155.
#
# equalities: R'x = r
# inequalities: S'x >= s
#
create_TMVN_sampler <- function(Q, perm=NULL,
                                mu=NULL, Xy=NULL, update.Q=FALSE, update.mu=update.Q,
                                name="x", coef.names=NULL,
                                R=NULL, r=NULL, S=NULL, s=NULL, lower=NULL, upper=NULL,
                                method=NULL, reduce=(method=="Gibbs" && !is.null(R)),
                                T.HMC=pi/2, print.info=FALSE, sharpness=100, useV=FALSE) {

  if (name == "") stop("empty name")
  store_default <- function(prior.sampler=FALSE) name  # for direct use of create_TMVN_sampler

  if (!update.Q) {
    Q <- economizeMatrix(Q, symmetric=TRUE, drop.zeros=TRUE)
  }
  n <- nrow(Q)

  if (!is.null(coef.names)) {
    if (length(coef.names) != n) stop("incompatible length of 'coef.names'")
    coef.names <- list(coef.names)
    names(coef.names) <- name
  }

  if (is.null(R)) {
    rm(R, r)
    reduce <- FALSE
    eq <- FALSE
  } else {
    if (nrow(R) != n || ncol(R) > n) stop("incompatible constraint matrix 'R'")
    if (is.null(r)) {
      r <- rep.int(0, ncol(R))
    } else {
      r <- as.vector(r)
      if (length(r) != ncol(R)) stop("length of 'r' should equal the number of columns of 'R'")
    }
    # TODO check that R has full column rank
    eq <- TRUE
  }

  if (is.null(S)) {
    rm(S, s)
    ineq <- FALSE
  } else {
    if (!is.null(lower)) {
      # use lower and upper bounds, and translate to Sx >= s
      if (is.null(upper)) stop("lower and upper bound should be specified together")
      if (!is.null(s)) warning("argument 's' ignored in combination with 'lower' and 'upper'", immediate.=TRUE)
      if (length(lower) != length(upper) || any(lower >= upper)) stop("'lower' and 'upper' incompatible")
      if (length(lower) != ncol(S)) stop("'S' incompatible with 'lower' and 'upper'")
      S <- cbind(S, -S)
      s <- c(lower, -upper)
    }
    if (nrow(S) != n) stop("incompatible constraint matrix 'S'")
    ncS <- ncol(S)
    if (is.null(s)) {
      s <- rep.int(0, ncS)
    } else {
      s <- as.vector(s)
      if (length(s) != ncS) stop("length of 's' should equal the number of columns of 'S'")
    }
    ineq <- TRUE
  }
  rm(lower, upper)

  if (ineq) {
    method <- match.arg(method, c("HMC", "Gibbs", "softTMVN"))  # method direct not for inequalities
    if (method %in% c("HMC", "Gibbs")) {
      base_which <- base::which  # faster in case Matrix::which is loaded
      eps <- sqrt(.Machine$double.eps)
      negeps <- -eps
    }
  } else {
    method <- match.arg(method, c("direct", "HMC", "Gibbs"))
  }
  if (eq && method == "Gibbs" && !reduce) {
    warning("'reduce' has been set to TRUE for Gibbs method with equalities", immediate.=TRUE)
    reduce <- TRUE
  }

  if (reduce || (ineq && method == "Gibbs")) name.tr <- paste0(name, "_tr_")

  use.cholV <- FALSE
  if (reduce) {
    # transform to subspace defined by R

    # in this case assume Q is provided (for now)
    if (is.null(Q)) stop("precision matrix 'Q' should be provided when 'reduce=TRUE'")
    if (update.Q || update.mu) stop("'update.Q=TRUE' and 'update.mu=TRUE' not supported in combination with 'reduce=TRUE'")

    mu <- if (is.null(mu)) rep.int(0, n) else as.numeric(mu)
    if (length(mu) != n) stop("incompatible dimensions 'Q' and 'mu'")

    QRofR <- qr(economizeMatrix(R, sparse=FALSE))  # check the rank
    # transformation z = Q'x where Q is the orthogonal Q matrix of QRofR
    z1 <- solve(t(qr.R(QRofR, complete=FALSE)), r)  # first, fixed part of z
    QofR <- economizeMatrix(qr.Q(QRofR, complete=TRUE))
    n1 <- ncol(R)  # nr of equality constraints, n1 >= 1
    n2 <- n - n1
    if (n2 < 1L) stop("degenerate case: as many equality restrictions as variables")
    Q1 <- QofR[, seq_len(n1), drop=FALSE]
    seq.free <- seq_len(n2)
    Q2 <- QofR[, n1 + seq.free, drop=FALSE]
    rm(QofR)
    mu_z1 <- crossprod_mv(Q1, mu)
    mu_z2 <- crossprod_mv(Q2, mu)
    # NB cholQ will now refer to z2 subspace
    cholQ <- build_chol(crossprod_sym(Q2, Q))  # chol factor L
    mu_z2_given_z1 <- mu_z2 - cholQ$solve(crossprod_mv(Q2, (Q %m*v% (Q1 %m*v% (z1 - mu_z1)))))
    if (method == "softTMVN") {
      Q <- crossprod_sym(Q2, Q)
    } else {
      rm(Q)
    }
    # fixed part of x: backtransform mu_z2_given_z1
    x0 <- Q1 %m*v% z1 + Q2 %m*v% mu_z2_given_z1
    mu <- rep.int(0, n2)  # this now refers to the mean of the transformed variable z2 - mu_z2_given_z1

    if (ineq) {
      # transform s and S to the frame defined by z2 - mu_z2_given_z1 (with vanishing mean)
      s <- s - crossprod_mv(S, x0)
      S <- economizeMatrix(crossprod(Q2, S), drop.zeros=TRUE)
      if (method == "HMC") {
        # define VS as Q^-1 S = V S
        VS <- economizeMatrix(cholQ$solve(S), allow.tabMatrix=FALSE)
      }
    }

    rm(Q1, z1, mu_z1, mu_z2, mu_z2_given_z1)
    # also remove R, r, n, n1, n2 ? and seq.free if method != "Gibbs"

    # now we have
    #   x0 as offset
    #   chol_Q to draw MVN variates for z2 around its mean
    #   Q2 to backtransform

  } else {  # !reduce

    if (method == "direct" && !update.Q && n <= 1000L) {  # special case where we use cholV instead of cholQ
      tryCatch(
        suppressWarnings({
          V <- economizeMatrix(solve(Q), symmetric=TRUE)
          cholV <- economizeMatrix(chol(V))
        }),
        error = function(err) {
          # TODO using determinant is not a very good way to check for non-pd
          d <- determinant(Q, logarithm=FALSE)
          if (d$sign * d$modulus < .Machine$double.eps) {
            stop("Non-positive-definite matrix in model component `", name, "`. Perhaps increasing the coefficients' prior precision helps.")
          } else {
            stop(err)
          }
        }
      )
      # remove any cached Cholesky factorizations
      if (class(Q)[1L] == "dsCMatrix") attr(Q, "factors") <- list()
      if (class(V)[1L] == "dsCMatrix") attr(V, "factors") <- list()
      use.cholV <- TRUE
    } else {
      cholQ <- build_chol(Q, perm)
    }
    rm(perm)
    if (method != "softTMVN") rm(Q)
    if (is.null(mu)) {
      if (is.null(Xy)) {
        mu <- rep.int(0, n)
      } else {
        mu <- if (use.cholV) V %m*v% Xy else cholQ$solve(Xy)
      }
    } else {
      mu <- as.numeric(mu)
    }
    if (length(mu) != n) stop("incompatible dimensions 'Q' and 'mu'")

    if (ineq) {
      if (method == "HMC" && print.info) {
        # set up a vector of wall bounce counts
        bounces <- setNames(rep.int(0L, ncS), if (is.null(colnames(S))) seq_len(ncS) else colnames(S))
        display.n <- seq_len(min(10L, ncS))  # number of counts to display
      }
      # currently tabMatrix disallowed because used as solve rhs below
      S <- economizeMatrix(S, sparse=if (class(cholQ$cholM)[1L] == "matrix") FALSE else NULL, allow.tabMatrix=FALSE)
      if (method == "HMC" && !update.Q) {
        # define VS as Q^-1 S = V S
        VS <- economizeMatrix(cholQ$solve(S), sparse=if (is.matrix(cholQ$cholM)) FALSE else NULL, allow.tabMatrix=FALSE)  # NB alternatively transform after projection
      }
    }
    if (method == "Gibbs") {
      n2 <- n
      seq.free <- seq_len(n)
    }

    if (eq) {
      # update.Q, dense cholQ --> faster to have dense R (and S) as well
      R <- economizeMatrix(R, sparse=if (update.Q && class(cholQ$cholM)[1L] == "matrix") FALSE else NULL, allow.tabMatrix=FALSE)
      if (update.Q) {
        bigR <- prod(dim(R)) > 100000L && ncol(R) > 5L  # if bigR a different (typically faster for large R) projection method is used
        #cholRVR <- build_chol(crossprod_sym2(R, cholQ$solve(R)))  # cholesky template object for (hopefully) faster projection
        cholRVR <- build_chol(crossprod_sym2(cholQ$solve(R, system="L", systemP=TRUE)))
      } else {
        if (method != "softTMVN") {
          VR <- if (use.cholV) economizeMatrix(V %*% R) else cholQ$solve(R)
          VR.RVRinv <- economizeMatrix(VR %*% solve(crossprod_sym2(R, VR)))
          rm(VR)
        }
      }
      if (method == "HMC" && !update.Q) {
        # project mean mu on eq constraint surface
        mu <- mu + VR.RVRinv %m*v% (r - crossprod_mv(R, mu))
        if (ineq) {
          # project inequality matrix on equality constraint surface, use Q as a metric
          VS <- economizeMatrix(VS - VR.RVRinv %*% crossprod(R, VS), sparse=if (is.matrix(cholQ$cholM)) FALSE else NULL, allow.tabMatrix=FALSE)
        }
      }
    }  # END if (eq)
  }  # END !reduce
  if (use.cholV && !update.mu) rm(V)
  rm(Xy)

  if (ineq && method == "Gibbs") {
    # transform to unit covariance matrix frame
    # inequalities U'v >= u
    # drop very small numbers in U as they give problems in the Gibbs sampler
    # use transpose for faster dgCMatrix col access
    #Ut <- economizeMatrix(t(zapsmall(cholQ$solve(S, system="L", systemP=TRUE), digits=10L)), drop.zeros=TRUE)
    Ut <- economizeMatrix(t(drop0(cholQ$solve(S, system="L", systemP=TRUE), tol=1e-10)))
    u <- s - crossprod_mv(S, mu)
    sparse <- class(Ut)[1L] == "dgCMatrix"
  }

  if (method == "HMC") {
    # Hamiltonian H = 1/2 x'Mx - g'x + 1/2 p'M^{-1}p
    # precision matrix M, covariance matrix M^{-1}, mu=M^{-1}g must obey equality restrictions
    # Hamilton's equations: dx/dt = M^{-1}p, p = M dx/dt
    #                       dp/dt = - Mx + g
    if (length(T.HMC) == 2L) {
      drawT <- TRUE
      runif_T <- eval(substitute(function() runif(1L, l, u), list(l=T.HMC[1L], u=T.HMC[2L])))
    } else {
      if (length(T.HMC) != 1L || T.HMC <= 0) stop("invalid 'T.HMC'")
      drawT <- FALSE
    }
    if (ineq) {
      if (update.Q) {
        simplified <- FALSE
      } else {
        refl.fac <- 2 / colSums(S * VS)  # 2 / vector of normal' Q normal for all inequalities
        s.adj <- s - crossprod_mv(S, mu)  # if positive, mu violates inequalities
        abs.s.adj <- abs(s.adj)
        simplified <- !eq && identical(VS, S)  # simplified=TRUE is the case described in Pakman and Paninski: identity Q and no equality constraints
        if (simplified) rm(VS)
      }
      twopi <- 2 * pi
    }
  }

  zero.mu <- !update.mu && all(mu == 0)
  if (zero.mu) rm(mu)

  if (method == "softTMVN") {
    if (.opts$PG.approx) {
      mPG <- as.integer(.opts$PG.approx.m)
      if (!length(mPG) %in% c(1L, n)) stop("invalid value for option 'PG.approx.m'")
      rPolyaGamma <- function(b, c) CrPGapprox(ncS, b, c, mPG)
    } else {
      if (!requireNamespace("BayesLogit", quietly=TRUE)) stop("please install package 'BayesLogit' and try again")
      rpg <- BayesLogit::rpg
      rPolyaGamma <- function(b, c) rpg(ncS, b, c)
    }
    if (useV) {
      # 'dual' version of MVN sampling
      V <- chol2inv(chol(Q))
      cholQ <- build_chol(Q)
      rm(Q)
      # define Diagonal matrix template to hold omega.tilde for use in sparse matrix templated sum
      D.omega.tilde.inv <- Diagonal(x=runif(ncS, 0.9, 1.1))
      matsum <- make_mat_sum(M0=economizeMatrix(crossprod_sym(S, V), symmetric=TRUE, drop.zeros=TRUE), M1=D.omega.tilde.inv)
      ch <- build_chol(matsum(D.omega.tilde.inv))
      VS <- economizeMatrix(V %*% S, drop.zeros=TRUE)
      if (eq && !reduce) {
        VR <- economizeMatrix(V %*% R, drop.zeros=TRUE)
        RVR <- economizeMatrix(crossprod(R, VR), drop.zeros=TRUE, symmetric=TRUE)
        SVR <- economizeMatrix(crossprod(S, VR), drop.zeros=TRUE)
        matsum_RVxR <- make_mat_sum(M0=RVR, M1=-crossprod_sym2(SVR, ch$solve(SVR)))
        cholRVxR <- build_chol(matsum_RVxR(M1=crossprod_sym2(SVR, ch$solve(SVR)), w1=-1))
      }
    } else {
      # more convenient to use transpose of S for use in crossprod_sym
      St <- t(S)
      rm(S)
      if (!zero.mu) Qmu <- Q %m*v% mu
      matsum <- make_mat_sum(M0=Q, M1=crossprod_sym(St, runif(nrow(St), 0.9, 1.1)))
      ch <- build_chol(matsum(crossprod_sym(St, runif(nrow(St), 0.9, 1.1))))
      if (eq && !reduce) {
        cholRVR <- build_chol(crossprod_sym2(R, ch$solve(R)))
      }
    }
  }  # END if (method == "softTMVN")

  if (method == "direct") {
    draw <- function(p) {}
    if (update.Q)
      draw <- add(draw, quote(cholQ$update(Q, Imult)))
    if (zero.mu)
      if (use.cholV)
        draw <- add(draw, bquote(coef <- drawMVN_cholV(.(n), cholV, scale)))
      else
        draw <- add(draw, quote(coef <- drawMVN_cholQ(cholQ, sd=scale)))
    else {
      if (update.mu)
        if (use.cholV)
          draw <- add(draw, bquote(coef <- V %m*v% Xy + drawMVN_cholV(.(n), cholV, scale)))
        else
          draw <- add(draw, quote(coef <- drawMVN_cholQ(cholQ, Xy, sd=scale)))
      else
        if (use.cholV)
          draw <- add(draw, bquote(coef <- mu + drawMVN_cholV(.(n), cholV, scale)))
        else
          draw <- add(draw, quote(coef <- mu + drawMVN_cholQ(cholQ, sd=scale)))
    }
    if (eq && !reduce) {
      if (update.Q) {
        if (bigR) {
          draw <- add(draw, quote(cholV.R <- cholQ$solve(R, system="L", systemP=TRUE)))
          draw <- add(draw, quote(cholRVR$update(crossprod_sym2(cholV.R))))
          draw <- add(draw, bquote(p[[.(name)]] <- coef + cholQ$solve(R %m*v% cholRVR$solve(r - crossprod_mv(R, coef)))))
          # alternative (this seems not faster, unless R is already dense, as cholV.R is usually less sparse than R):
          # draw <- add(draw, bquote(p[[.(name)]] <- coef + cholQ$solve(cholV.R %m*v% cholRVR$solve(r - crossprod_mv(R, coef)), system="Lt", systemP=TRUE)))
        } else {
          draw <- add(draw, quote(VR <- cholQ$solve(R)))
          draw <- add(draw, quote(cholRVR$update(crossprod_sym2(R, VR))))
          draw <- add(draw, bquote(p[[.(name)]] <- coef + VR %m*v% cholRVR$solve(r - crossprod_mv(R, coef))))
        }
      } else {
        draw <- add(draw, bquote(p[[.(name)]] <- coef + VR.RVRinv %m*v% (r - crossprod_mv(R, coef))))
      }
    } else {
      if (reduce) {  # backtransform
        draw <- add(draw, bquote(p[[.(name)]] <- x0 + Q2 %m*v% coef))
      } else {
        draw <- add(draw, bquote(p[[.(name)]] <- coef))
      }
    }
    draw <- add(draw, quote(p))
    # set function signature (avoiding check NOTE)
    if (update.Q) {
      # total precision matrix to use in chol is Q + Imult*I
      formals(draw) <- c(alist(p=), alist(scale=1), alist(Q=), alist(Imult=0), alist(Xy=))
    } else {
      if (update.mu)
        formals(draw) <- c(alist(p=), alist(scale=1), alist(Xy=))
      else
        formals(draw) <- c(alist(p=), alist(scale=1))
    }
  }  # END if (method == "direct")

  if (method == "softTMVN") {
    draw <- function(p) {
      # 1. draw Polya-Gamma mixture precisions
      if (reduce)  # equalities have already been taken care of
        x <- p[[name.tr]]
      else
        x <- p[[name]]
      discr <- if (useV) crossprod_mv(S, x) - s else St %m*v% x - s
      omega <- rPolyaGamma(1, sharpness * discr)
      # 2. draw vector of coefficients
      omega.tilde <- sharpness * sharpness * omega
      if (useV) {
        attr(D.omega.tilde.inv, "x") <- 1 / omega.tilde
        XX_V <- matsum(D.omega.tilde.inv)
        ch$update(XX_V)
        alpha <- sharpness * s + 0.5/omega
        if (zero.mu)
          y1 <- drawMVN_cholQ(cholQ)
        else
          y1 <- mu + drawMVN_cholQ(cholQ)
        y2 <- Crnorm(length(s)) * sqrt(1/omega)
        if (reduce) {
          p[[name.tr]] <- y1 + VS %m*v% ch$solve((alpha - sharpness * crossprod_mv(S, y1) - y2) * (1 / sharpness))
          p[[name]] <- x0 + Q2 %m*v% p[[name.tr]]
        } else {
          coef <- y1 + VS %m*v% ch$solve((alpha - sharpness * crossprod_mv(S, y1) - y2) * (1 / sharpness))
          if (eq) {
            cholRVxR$update(matsum_RVxR(M1=crossprod_sym2(SVR, ch$solve(SVR)), w1=-1))
            temp <- V %m*v% (R %m*v% cholRVxR$solve(r - crossprod_mv(R, coef)))
            p[[name]] <- coef + temp - V %m*v% (S %m*v% ch$solve(crossprod_mv(S, temp)))
          } else
            p[[name]] <- coef
        }
      } else {
        XX <- crossprod_sym(St, omega.tilde)
        XX_Q <- matsum(XX)
        ch$update(XX_Q)
        alpha <- 0.5 * sharpness + s * omega.tilde
        if (zero.mu)
          Xy <- crossprod_mv(St, alpha)
        else
          Xy <- crossprod_mv(St, alpha) + Qmu
        if (reduce) {
          p[[name.tr]] <- drawMVN_cholQ(ch, Xy)
          p[[name]] <- x0 + Q2 %m*v% p[[name.tr]]
        } else {
          coef <- drawMVN_cholQ(ch, Xy)
          if (eq) {
            VR <- ch$solve(R)
            cholRVR$update(crossprod_sym2(R, VR))
            p[[name]] <- coef + VR %m*v% cholRVR$solve(r - crossprod_mv(R, coef))
          } else
            p[[name]] <- coef
        }
      }
      p
    }
  }  # END if (method == "softTMVN")

  if (method == "Gibbs") {
    if (ineq) {
      draw <- function(p) {
        # draw from truncated univariate normal full conditionals
        v <- p[[name.tr]]
        ustar <- u - Ut %m*v% v
        if (sparse) {
          v <- Crtmvn_Gibbs(v, Ut, ustar, eps)
        } else {
          for (i in seq.free) {
            Ui <- Ut[, i, drop=TRUE]
            ustar <- ustar + Ui * v[i]
            Ui_plus <- base_which(Ui > eps)
            a <- if (length(Ui_plus)) max(ustar[Ui_plus] / Ui[Ui_plus]) else -Inf
            Ui_minus <- base_which(Ui < negeps)
            b <- if (length(Ui_minus)) min(ustar[Ui_minus] / Ui[Ui_minus]) else Inf
            if (a < b)
              v[i] <- Crtuvn(a, b)
            else {
              # this seems a numerically stable way to deal with numerical inaccuracy:
              if (a < v[i])
                v[i] <- a
              else if (b > v[i])
                v[i] <- b
            }
            ustar <- ustar - Ui * v[i]
          }
        }
        p[[name.tr]] <- v
      }
    } else {
      draw <- function(p) {
        v <- Crnorm(n2)
      }
    }
    if (eq) {
      # TODO precompute Q2 %*% Lt^-1
      draw <- add(draw, bquote(p[[.(name)]] <- x0 + Q2 %m*v% cholQ$solve(v, system="Lt", systemP=TRUE)))
    } else {
      if (zero.mu)
        draw <- add(draw, bquote(p[[.(name)]] <- cholQ$solve(v, system="Lt", systemP=TRUE)))
      else
        draw <- add(draw, bquote(p[[.(name)]] <- mu + cholQ$solve(v, system="Lt", systemP=TRUE)))
    }
    draw <- add(draw, quote(p))
  }  # END if (method == "Gibbs") {

  if (method == "HMC") {
    draw <- function(p) {
      if (update.Q) {
        VS <- cholQ$solve(S)
        if (eq) {
          if (bigR) {
            cholV.R <- cholQ$solve(R, system="L", systemP=TRUE)
            cholRVR$update(crossprod_sym2(cholV.R))
            VS <- VS - cholQ$solve(R %*% cholRVR$solve(crossprod(R, VS)))
          } else {
            VR <- cholQ$solve(R)
            cholRVR$update(crossprod_sym2(R, VR))
            VS <- VS - VR %*% cholRVR$solve(crossprod(R, VS))
          }
        }
        refl.fac <- 2 / colSums(S * VS)  # 2 / vector of normal' Q normal for all inequalities
      }
      if (update.mu) {
        mu <- cholQ$solve(Xy)
        if (eq)
          if (update.Q)
            if (bigR)
              mu <- mu + cholQ$solve(R %m*v% cholRVR$solve(r - crossprod_mv(R, mu)))
            else
              mu <- mu + VR %m*v% cholRVR$solve(r - crossprod_mv(R, mu))
            else
              mu <- mu + VR.RVRinv %m*v% (r - crossprod_mv(R, mu))
            s.adj <- s - crossprod_mv(S, mu)  # if positive, mu violates inequalities
            abs.s.adj <- abs(s.adj)
      }

      # 1. draw p from N(0, M^{-1}); instead draw v=dx/dt from N(0, M)
      #    this is the initial velocity (starting from final position x of previous draw)
      v <- drawMVN_cholQ(cholQ, sd=scale)  # draw velocity at start time t=0
      if (reduce) {  # equalities have already been taken care of
        x <- p[[name.tr]]
      } else {
        x <- p[[name]]
        if (eq) {
          if (update.Q) {
            if (bigR)
              v <- v - cholQ$solve(R %m*v% cholRVR$solve(crossprod_mv(R, v)))
            else
              v <- v - VR %m*v% cholRVR$solve(crossprod_mv(R, v))
          } else {
            v <- v - VR.RVRinv %m*v% crossprod_mv(R, v)
          }
        }
      }
      if (drawT) T.HMC <- runif_T()
      # solution: x(t) = mu + v0*sin(t) + (x0 - mu)*cos(t)
      # inequalities: S'mu + S'v sin(t) + S'(x0 - mu) cos(t) >= s
      #           --> S'mu + u cos(t + phi) >= s
      if (ineq) {
        if (print.info) {
          viol <- crossprod_mv(S, x) < s
          if (any(viol))
            cat("\nviolated constraints:", names(bounces)[viol])  # print violated constraints
        }
        repeat {
          Sv <- crossprod_mv(S, v)
          if (zero.mu)
            Sx <- crossprod_mv(S, x)
          else
            Sx <- crossprod_mv(S, x - mu)
          u <- sqrt(Sv^2 + Sx^2)

          # which inequality walls can be reached by the exact solution?
          i_hit <- base_which(u > abs.s.adj)
          # remove wall passings into the feasible region
          # at t=0: crossprod_mv(S, x) < s
          n_hit <- length(i_hit)
          if (n_hit == 0L) break  # no walls can be hit -> immediately return end position

          #phi <- atan2(-Sv[i_hit], Sx[i_hit])  # -pi < phi < pi
          ui <- 1 / u[i_hit]
          phi <- -sign(Sv[i_hit]) * acos(Sx[i_hit] * ui)  # faster than atan2
          # compute collision times, for which u cos(t + phi) = s - S' mu
          # 0 <= acos(x) <= pi
          # solutions:
          #  - phi + acos(rhs/u) + 2 pi n
          #  - phi - acos(rhs/u) + 2 pi n
          t_hit <- acos(s.adj[i_hit] * ui)

          # -pi < t_hit - phi < 2 pi
          # -2 pi < -t_hit - phi < pi
          t_hit <- c( t_hit - phi, - t_hit - phi )
          ind <- base_which(t_hit < negeps)
          t_hit[ind] <- t_hit[ind] + twopi
          t_hit[t_hit < eps] <- 0

          # ignore inequality boundaries at which the particle currently is (including previously hit walls) while moving to the interior
          h0 <- which.min(t_hit)
          th <- t_hit[h0]
          h <- if (h0 > n_hit) h0 - n_hit else h0  # inequality w.r.t i_hit where first collision occurs
          while (th == 0 && Sv[i_hit[h]] > 0) {
            # the collision has already occurred, and we should ignore th=0 otherwise the particle would stick to the wall
            t_hit[h0] <- twopi
            h0 <- which.min(t_hit)
            th <- t_hit[h0]
            h <- if (h0 > n_hit) h0 - n_hit else h0  # inequality w.r.t i_hit where first collision occurs
          }
          if (th >= T.HMC) break  # simulation time ends before next wall hit

          # next wall hit happens within simulation period
          h <- i_hit[h]  # index w.r.t. all inequalities
          if (print.info) bounces[h] <- bounces[h] + 1L
          if (th > 0) {
            costh <- cos(th)
            sinth <- sin(th)
            v0 <- v
            if (zero.mu) {
              v <- v0 * costh - x * sinth  # velocity at time th
              x <- v0 * sinth + x * costh  # x at time th, new starting position
            } else {
              delta <- x - mu
              v <- v0 * costh - delta * sinth  # velocity at time th
              x <- mu + v0 * sinth + delta * costh  # x at time th, new starting position
            }
            T.HMC <- T.HMC - th  # update remaining simulation time
          }

          Sh <- get_col(S, h)  # TODO use sparse vector if possible
          vproj <- dotprodC(Sh, v)
          if (vproj < 0) {  # a hit
            alpha <- vproj * refl.fac[h]
            normal <- if (simplified) Sh else get_col(VS, h)
            v <- v - alpha * normal
          }  # else do nothing: passage to the feasible side of the inequality wall

        }  # END repeat
        if (print.info) {
          most_bounces <- sort(bounces, decreasing=TRUE)[display.n]
          cat("\nmost bounces:", paste0(names(most_bounces), " (", most_bounces, ") "))
        }
      }  # END if (ineq)

      if (reduce) {
        if (zero.mu)
          p[[name.tr]] <- v * sin(T.HMC) + x * cos(T.HMC)  # final state (new draw for x)
        else
          p[[name.tr]] <- mu + v * sin(T.HMC) + (x - mu) * cos(T.HMC)  # final state (new draw for x)
        p[[name]] <- x0 + Q2 %m*v% p[[name.tr]]
      } else {
        if (zero.mu)
          p[[name]] <- v * sin(T.HMC) + x * cos(T.HMC)  # final state (new draw for x)
        else
          p[[name]] <- mu + v * sin(T.HMC) + (x - mu) * cos(T.HMC)  # final state (new draw for x)
      }

      p
    }  # END function draw
    # set function signature (avoiding check NOTE)
    if (update.Q) {
      formals(draw) <- c(alist(p=), alist(scale=1), alist(cholQ=), alist(Xy=))
    } else {
      if (update.mu)
        formals(draw) <- c(alist(p=), alist(scale=1), alist(Xy=))
      else
        formals(draw) <- c(alist(p=), alist(scale=1))
    }
  }  # END if (method == "HMC")

  # start function
  if (ineq && method == "Gibbs") {
    start <- function(p=list(), scale=1) {
      if (is.null(p[[name.tr]])) {
        if (is.null(p[[name]])) {
          if (!requireNamespace("lintools", quietly=TRUE)) stop("please install package lintools and try again")
          temp <- -as(Ut, "dgTMatrix")
          # sparse version
          A <- data.frame(
            row = temp@i,  # constraint index
            col = temp@j,  # variable index
            coef = temp@x
          )
          x0 <- rnorm(n2, sd=scale)
          for (i in 1:10) {  # at most 10 trials
            epsilon <- scale * rexp(ncS)
            res <- lintools::sparse_project(x=x0, A, b=-(u + epsilon), neq=0L, base=0L, sorted=FALSE)
            scale <- 0.3 * scale
            if (anyNA(res$x)) {
              x0 <- rnorm(n2, sd=scale)
            } else {
              if (all(Ut %m*v% res$x >= u)) break
              x0 <- res$x + rnorm(n2, sd=scale)
            }
          }
          p[[name.tr]] <- res$x
        } else {
          # a start value at the original scale is provided; transform to projected form
          # TODO check that the start value satisfies all constraints
          if (eq) {
            z2.start <- crossprod_mv(Q2, p[[name]] - x0)
          } else {
            if (zero.mu)
              z2.start <- p[[name]]
            else
              z2.start <- p[[name]] - mu
          }
          p[[name.tr]] <- cholQ$crossprodL(z2.start)
        }
      }
      p
    }
  } else {
    # TODO here too use lintools to find a better starting value
    if (reduce) {
      start <- function(p=list(), scale=1) {
        if (!is.null(p[[name.tr]])) {
          return(p)
        } else {
          if (!is.null(p[[name]])) {
            p[[name.tr]] <- crossprod_mv(Q2, p[[name]] - x0)
            return(p)
          }
        }
      }
    } else {
      start <- function(p=list(), scale=1) {
        if (!is.null(p[[name]])) {
          p[[name]] <- as.numeric(p[[name]])
          return(p)
        }
      }
    }
    # TODO check that user-provided start value has right format and(?) obeys constraints
    if (use.cholV) {
      if (zero.mu) {
        start <- add(start, bquote(coef <- drawMVN_cholV(.(n), cholV, scale)))
      } else {
        start <- add(start, bquote(coef <- mu + drawMVN_cholV(.(n), cholV, scale)))
      }
    } else {
      if (zero.mu || update.Q) {
        start <- add(start, quote(coef <- drawMVN_cholQ(cholQ, sd=scale)))
      } else {
        start <- add(start, quote(coef <- mu + drawMVN_cholQ(cholQ, sd=scale)))
      }
    }
    if (eq && !reduce) {
      if (update.Q) {
        # TODO: do we really need a start funtion in this case, or call draw
        if (bigR) {
          start <- add(start, quote(cholV.R <- cholQ$solve(R, system="L", systemP=TRUE)))
          start <- add(start, quote(cholRVR$update(crossprod_sym2(cholV.R))))
          start <- add(start, bquote(p[[.(name)]] <- coef + cholQ$solve(R %m*v% cholRVR$solve(r - crossprod_mv(R, coef)))))
        } else {
          start <- add(start, quote(VR <- cholQ$solve(R)))
          start <- add(start, quote(cholRVR$update(crossprod_sym2(R, VR))))
          start <- add(start, bquote(p[[.(name)]] <- coef + VR %m*v% cholRVR$solve(r - crossprod_mv(R, coef))))
        }
      } else {
        if (method == "softTMVN") {
          if (useV) {
            start <- add(start, quote(VxR <- VR - VS %*% ch$solve(SVR)))  # TODO for sparse matrices use matsum template
            start <- add(start, quote(cholRVxR$update(crossprod(R, VxR))))
            start <- add(start, bquote(p[[.(name)]] <- coef + VxR %m*v% cholRVxR$solve(r - crossprod_mv(R, coef))))
          } else
            start <- add(start, bquote(p[[.(name)]] <- coef + ch$solve(R) %m*v% cholRVR$solve(r - crossprod_mv(R, coef))))
        } else {
          # NB it doesn't matter if mu is already projected
          start <- add(start, bquote(p[[.(name)]] <- coef + VR.RVRinv %m*v% (r - crossprod_mv(R, coef))))
        }
      }
    } else {
      if (reduce) {  # backtransform
        start <- add(start, bquote(p[[.(name.tr)]] <- coef))
        start <- add(start, bquote(p[[.(name)]] <- x0 + Q2 %m*v% coef))
      } else {
        start <- add(start, bquote(p[[.(name)]] <- coef))
      }
    }
    #check_constraints(p[[name]])
    start <- add(start, quote(p))
  }

  #check_constraints <- function(x) {
  #  if (eq) {
  #    if (!isTRUE(all.equal(crossprod_mv(R, x), r))) stop("equality restriction(s) violated")
  #  }
  #  if (ineq) {
  #    if (!all(crossprod_mv(S, x) >= s - eps)) stop("inequality restriction(S) violated")
  #  }
  #}

  # TODO remove more unused quantities
  rm(n, use.cholV)
  if (method != "HMC") rm(T.HMC)
  if (method != "HMC" || !ineq) rm(print.info)

  environment()
}
