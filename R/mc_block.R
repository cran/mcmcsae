
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
  X <- matrix(0, e$n, 0L)  # block design matrix
  ind <- 0L
  for (mc in mcs) {
    if (mc$type == "gen" && mc$gl)
      X <- cbind(X, mc$X, zeroMatrix(e$n, mc$glp$q))
    else
      X <- cbind(X, mc$X)
    vec_list[[mc$name]] <- (ind + 1L):ncol(X)
    ind <- ncol(X)
  }
  rm(ind)
  X <- economizeMatrix(X)

  # template for (updating) blocked precision matrix
  # each mc$Q for gen is either ddi or dsC; the same holds true for mc$Q0 for reg and mec
  get_Q <- function(mc) if (mc$type == "gen") mc$Q else mc$Q0
  if (all(sapply(mcs, function(mc) class(get_Q(mc))[[1L]] == "ddiMatrix"))) {
    QT <- Cdiag(unlist(lapply(mcs, function(mc) ddi_diag(get_Q(mc))), use.names=FALSE))
  } else {
    x <- NULL
    i <- NULL
    size <- 0L
    p <- 0L
    for (mc in mcs) {
      Q <- get_Q(mc)
      if (class(Q)[1L] == "ddiMatrix") {
        x <- c(x, ddi_diag(Q))
        i <- c(i, size + seq_len(nrow(Q)) - 1L)
        p <- c(p, p[length(p)] + seq_len(nrow(Q)))
      } else {
        x <- c(x, Q@x)
        i <- c(i, size + Q@i)
        p <- c(p, p[length(p)] + Q@p[-1L])
      }
      size <- size + nrow(Q)
    }
    QT <- new("dsCMatrix", i=i, p=p, x=x, uplo="U", Dim=c(size, size))
    rm(size, x, i, p)
  }
  rm(get_Q)
  # individual Q matrices no longer needed (we still have kron_prod closures)
  for (mc in mcs) if (mc$type == "gen") rm("Q", envir=mc)
  get_Qvector <- function(p, tau) {
    Qvector <- NULL
    for (mc in mcs)
      if (mc$type == "gen")
        Qvector <- c(Qvector, p[[mc$name_Q]])
      else
        Qvector <- c(Qvector, tau * mc$Q0@x)
    Qvector
  }

  if (is.null(e$control$CG)) {
    if (e$modeled.Q)
      XX <- crossprod_sym(X, crossprod_sym(Diagonal(x=runif(e$n, 0.9, 1.1)), e$Q0))
    else
      XX <- economizeMatrix(crossprod_sym(X, e$Q0), symmetric=TRUE, drop.zeros=TRUE)

    # derive constraint matrix, if any
    if (any(sapply(mcs, function(mc) !is.null(mc$R)))) {
      R <- zeroMatrix(0L, 0L)
      for (mc in mcs) {
        if (is.null(mc$R)) {
          if (mc$type == "gen" && mc$gl)
            R <- rbind(R, zeroMatrix(mc$q + mc$glp$q, ncol(R)))
          else
            R <- rbind(R, zeroMatrix(mc$q, ncol(R)))
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
      S <- zeroMatrix(0L, 0L)
      for (mc in mcs) {
        if (is.null(mc$S))
          S <- rbind(S, zeroMatrix(mc$q, ncol(S)))
        else
          S <- bdiag(S, mc$S)
      }
      if (nrow(S) != ncol(X)) stop("incompatible dimensions of design and constraint matrices")
      # TODO remove individual S matrices as they are not needed in the single block sampler
      #      add support for additional constraints defined over the whole coefficient vector
    } else {
      S <- NULL
    }

    sparse_template(environment(), update.XX=e$modeled.Q || any(sapply(mcs, `[[`, "type") == "mec"),
                    add.outer.R=e$control$add.outer.R, prior.only=e$prior.only)

  } else {
    # TODO check that the only constraints are GMRF constraints
    CG_sampler <- do.call(setup_CG_sampler, c(list(mbs=mcs, X=X, sampler=e), e$control$CG))
  }

  if (e$compute.weights) {
    # form the q_all x m matrix corresponding to the linear predictor as represented componentwise in linpred
    # TODO do we really need to store both X and t(X) in this case?
    if (is.null(e$linpred)) {
      linpred <- economizeMatrix(t(X), strip.names=FALSE)
    } else {
      linpred <- economizeMatrix(t(do.call(cbind, lapply(e$linpred[names(vec_list)], function(x) environment(x)$Xnew))), allow.tabMatrix=FALSE)
    }
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
  if (!is.null(S)) {  # need to reconstruct coef_ as input to TMVN sampler
    # TODO store coef_ component and only replace the subcomponents with PX
    #      and check whether this works in case of gen component with gl=TRUE
    draw <- add(draw, quote(
      for (mc in mcs) p[["coef_"]][vec_list[[mc$name]]] <- p[[mc$name]]
    ))
  }
  # update mec component columns of X
  # TODO more efficient update of only those elements that can change (for dgC or matrix X)
  for (mc in mcs)
    if (mc$type == "mec")
      draw <- add(draw, bquote(X[, vec_list[[.(mc$name)]]] <- p[[.(mc$name_X)]]))
  if (e$single.block && !e$modeled.Q && !any(sapply(mcs, `[[`, "type") == "mec") && e$family$link != "probit")
    Xy <- crossprod_mv(X, e$Q0 %m*v% e$y_eff())
  else
    draw <- add(draw, quote(Xy <- crossprod_mv(X, e$Q_e(p))))
  if (is.null(e$control$CG)) {
    if (e$modeled.Q) {
      if (e$Q0.type == "symm")
        draw <- add(draw, quote(XX <- crossprod_sym(X, p[["QM_"]])))
      else {
        cps_template <- NULL
        if (inherits(X, "dgCMatrix")) {
          tryCatch(
            cps_template <- sparse_crossprod_sym_template(X),
            error = function(e) {
              # template too large
              NULL
            }
          )
        }
        if (is.null(cps_template)) {
          draw <- add(draw, quote(XX <- crossprod_sym(X, p[["Q_"]])))
        } else {
          draw <- add(draw, quote(XX <- cps_template(p[["Q_"]])))
        }
      }
    } else if (any(sapply(mcs, `[[`, "type") == "mec")) {
      draw <- add(draw, quote(XX <- crossprod_sym(X, e$Q0)))
    }
    draw <- add(draw, quote(Qlist <- update(XX, QT, 1, 1/tau)))
    draw <- add(draw, bquote(coef <- MVNsampler$draw(p, .(if (e$sigma.fixed) 1 else quote(p[["sigma_"]])), Q=Qlist$Q, Imult=Qlist$Imult, Xy=Xy)[[.(name)]]))
  } else {
    draw <- add(draw, quote(CGstart <- vector("numeric", ncol(X))))
    draw <- add(draw, quote(for (mc in mcs) CGstart[vec_list[[mc$name]]] <- p[[mc$name]]))
    draw <- add(draw, quote(coef <- CG_sampler$draw(p, Xy, X, QT, e, start=CGstart)))
  }
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
    # TODO solve-sparse method that returns dense
    draw <- add(draw, quote(p$weights_ <- X %m*m% as.matrix(MVNsampler$cholQ$solve(linpred))))
    if (e$modeled.Q) {
      if (e$Q0.type == "symm")
        draw <- add(draw, quote(p$weights_ <- p[["QM_"]] %m*m% p[["weights_"]]))
      else
        draw <- add(draw, quote(p$weights_ <- p[["Q_"]] * p[["weights_"]]))
    } else {
      if (e$Q0.type != "unit") {
        draw <- add(draw, quote(p$weights_ <- e$Q0 %m*m% p[["weights_"]]))
      }
    }
  }
  draw <- add(draw, quote(p))
  # END draw function

  start <- function(p) {}
  if (is.null(e$control$CG)) {
    start <- add(start, bquote(coef <- MVNsampler$start(p, e$scale_sigma)[[.(name)]]))
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
  } else {
    start <- add(start, quote(
      for (mc in mcs) {
        if (mc$type == "gen" && mc$fastGMRFprior)
          p[[mc$name]] <- mc$rprior(p)[[mc$name]]
        else
          p[[mc$name]] <- Crnorm(mc$q, sd=e$scale_sigma)
      }
    ))
  }
  # TODO check formats of user-provided start values
  start <- add(start, quote(p))

  rm(mc)
  environment()
}
