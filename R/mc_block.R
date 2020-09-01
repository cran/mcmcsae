
# mcs: sublist of mod components to be sampled in a block
# e: environment of create_sampler
create_mc_block <- function(mcs, e=parent.frame()) {
  type <- "block"
  name <- "coef_"
  debug <- any(sapply(mcs, `[[`, "debug"))
  vec_list <- list()

  if (.opts$auto.order.block) {
    # order the components such that sparse matrices come first; may help find a better Cholesky permutation
    o <- unname(which(vapply(mcs, function(m) isDiagonal(m$X), TRUE)))  # start with diagonal matrices
    if (length(o))
      o <- c(o, seq_along(mcs)[-o][order(vapply(mcs[-o], function(m) sparsity(m$X), 1), decreasing=TRUE)])
    else
      o <- order(vapply(mcs, function(m) sparsity(m$X), 1), decreasing=TRUE)
    mcs <- mcs[o]
    rm(o)
  }

  # start with empty matrices
  X <- matrix(0, e$n, 0L)  # block design matrix; used to produce (weighted) X'X matrix
  QT <- matrix(0, 0L, 0L)  # (template for) joint prior precision matrix
  ind <- 0L
  for (mc in mcs) {
    if (mc$type == "gen" && mc$gl)
      X <- cbind(X, cbind(mc$X, Matrix(0, e$n, mc$glp$q)))
    else
      X <- cbind(X, mc$X)
    if (mc$type == "reg") {
      QT <- bdiag(QT, mc$Q0)
    } else {
      QT <- bdiag(QT, mc$Q)
      rm("Q", envir=mc)  # individual Q matrices no longer needed (we still have kron_prod closures)
    }
    vec_list[[mc$name]] <- (ind + 1L):ncol(X)
    ind <- ncol(X)
  }
  X <- economizeMatrix(X)
  QT <- economizeMatrix(QT, sparse=TRUE, symmetric=TRUE)

  if (e$modeled.Q) {
    XX <- crossprod_sym(X, crossprod_sym(Diagonal(x=runif(e$n, 0.9, 1.1)), e$Q0))
  } else {
    XX <- economizeMatrix(crossprod_sym(X, e$Q0), symmetric=TRUE, drop.zeros=TRUE)
  }

  # the x-slot of QT is simply the concatenation of the x-slots of the model components, so updating is easy
  # derive constraint matrix, if any
  if (any(sapply(mcs, function(mc) !is.null(mc$R)))) {
    R <- Matrix(0, 0L, 0L)
    for (mc in mcs) {
      if (is.null(mc$R)) {
        if (mc$type == "gen" && mc$gl)
          R <- rbind(R, Matrix(0, mc$q + mc$glp$q, ncol(R)))
        else
          R <- rbind(R, Matrix(0, mc$q, ncol(R)))
      } else {
        if (mc$type == "gen" && mc$gl)
          R <- bdiag(R, mc$glp$R)
        else
          R <- bdiag(R, mc$R)
      }
    }
    if (nrow(R) != ncol(X)) stop("incompatible dimensions of design and constraint matrices")
    # TODO remove individual R matrices as they are not needed in the single block sampler
    #      add support for additional constraints defined over the whole coefficient vector
    # In the case of constraints the XX + Q matrix often becomes singular
    # sampling from such a IGMRF can be done by first adding a multiple of tcrossprod(R) to it (Rue and Held, 2005)
    # most convenient is to add it to XX;
    # But this is not always required. Better add as little as is necessary to get pd (takes up lots of memory for lengthy random walks ...)
    # Another option is to add a multiple of I to XX + Q and correct with MH
  } else {
    R <- NULL
  }

  if (any(sapply(mcs, function(mc) !is.null(mc$S)))) {
    S <- Matrix(0, 0L, 0L)
    for (mc in mcs) {
      if (is.null(mc$S))
        S <- rbind(S, Matrix(0, mc$q, ncol(S)))
      else
        S <- bdiag(S, mc$S)
    }
    if (nrow(S) != ncol(X)) stop("incompatible dimensions of design and constraint matrices")
    # TODO remove individual S matrices as they are not needed in the single block sampler
    #      add support for additional constraints defined over the whole coefficient vector
  } else {
    S <- NULL
  }

  sparse_template(environment(), modeled.Q=e$modeled.Q, add.outer.R=e$control$add.outer.R, prior.only=e$prior.only)

  if (e$single.block) {
    if (e$compute.weights) {
      # form the q_all x m matrix corresponding to the linear predictor as represented componentwise in linpred
      # TODO do we really need to store both X and t(X) in this case?
      if (all(sapply(e$linpred, is.character)) && all(unlist(e$linpred) == "X")) {
        linpred <- economizeMatrix(t(X), strip.names=FALSE)
      } else {
        linpred <- economizeMatrix(t(do.call(cbind, e$linpred[names(vec_list)])), strip.names=FALSE)
        # TODO as this matrix is passed to Cholesky's solve method it must be matrix or dgCMatrix
      }
    }
  }

  get_Qvector <- function(p, tau) {
    Qvector <- NULL
    for (mc in mcs)
      if (mc$type == "reg")
        Qvector <- c(Qvector, tau * mc$Q0@x)
      else
        Qvector <- c(Qvector, p[[mc$name_Q]])
    Qvector
  }

  if (e$prior.only) return(environment())

  # BEGIN draw function
  draw <- function(p) {}
  if (debug) draw <- add(draw, quote(browser()))
  if (e$single.block) {
    draw <- add(draw, quote(p$e_ <- e$y_eff()))
  } else {
    for (m in seq_along(mcs)) {
      # residuals could also be computed using the block mb$X,
      #   but this way it is usually faster due to more efficient matrix types
      # TODO composite matrix type
      if (e$e.is.res)
        draw <- add(draw, bquote(p$e_ <- p[["e_"]] + mcs[[.(m)]]$linpred(p)))
      else
        draw <- add(draw, bquote(p$e_ <- p[["e_"]] - mcs[[.(m)]]$linpred(p)))
    }
  }
  draw <- add(draw, bquote(tau <- .(if (e$sigma.fixed) 1 else quote(1 / p[["sigma_"]]^2))))
  # update the block-diagonal joint precision matrix
  draw <- add(draw, quote(attr(QT, "x") <- get_Qvector(p, tau)))
  if (e$modeled.Q) {
    if (e$Q0.type == "symm")
      draw <- add(draw, quote(XX <- crossprod_sym(X, p[["QM_"]])))
    else
      draw <- add(draw, quote(XX <- crossprod_sym(X, p[["Q_"]])))
  }
  draw <- add(draw, quote(Qlist <- update(XX, QT, 1, 1/tau)))
  if (e$single.block && !e$modeled.Q && e$family$link != "probit") {
    Xy <- crossprod_mv(X, e$Q0 %m*v% e$y_eff())  # precompute
  } else {
    draw <- add(draw, quote(Xy <- crossprod_mv(X, e$Q_e(p))))
  }
  draw <- add(draw, bquote(coef <- MVNsampler$draw(p, .(if (e$sigma.fixed) 1 else bquote(p[["sigma_"]])), Q=Qlist$Q, Imult=Qlist$Imult, Xy=Xy)[[.(name)]]))
  # split coef and assign to the separate coefficient batches
  draw <- add(draw, quote(
    for (mc in mcs) {
      if (mc$type == "gen" && mc$gl) {
        u <- coef[vec_list[[mc$name]]]
        p[[mc$name]] <- u[mc$i.v]
        p[[mc$name_gl]] <- u[mc$i.alpha]
      } else {
        p[[mc$name]] <- coef[vec_list[[mc$name]]]
      }
      if (e$e.is.res)
        p$e_ <- p[["e_"]] - mc$linpred(p)
      else
        p$e_ <- p[["e_"]] + mc$linpred(p)
    }
  ))
  if (e$compute.weights) {
    # TODO %m*m% function for matrix-matrix multiplication; solve-sparse method that returns dense
    draw <- add(draw, quote(p$weights_ <- X %*% as.matrix(MVNsampler$cholQ$solve(linpred))))
    if (e$modeled.Q) {
      if (e$Q0.type == "symm")
        draw <- add(draw, quote(p$weights_ <- p[["QM_"]] %*% p[["weights_"]]))
      else
        draw <- add(draw, quote(p$weights_ <- p[["Q_"]] * p[["weights_"]]))
    } else {
      if (e$Q0.type != "unit") {
        draw <- add(draw, quote(p$weights_ <- e$Q0 %*% p[["weights_"]]))
      }
    }
  }
  draw <- add(draw, quote(p))
  # END draw function

  start <- function(p) {}
  start <- add(start, bquote(coef <- MVNsampler$start(p, .(e$scale_y))[[.(name)]]))
  start <- add(start, quote(
    for (mc in mcs) {
      if (mc$type == "gen" && mc$gl) {
        u <- coef[vec_list[[mc$name]]]
        if (is.null(p[[mc$name]])) p[[mc$name]] <- u[mc$i.v]
        if (is.null(p[[mc$name_gl]])) p[[mc$name_gl]] <- u[mc$i.alpha]
      } else {
        if (is.null(p[[mc$name]])) p[[mc$name]] <- coef[vec_list[[mc$name]]]
      }
    }
  ))
  # TODO check formats of user-provided start values
  start <- add(start, quote(p))

  rm(ind, mc)
  environment()
}
