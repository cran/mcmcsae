# mc: model component
# update.XX: whether XX is updated in each draw
# add.outer.R: whether outer product of constraints should be added (for positive definiteness in blocked components)
# this function is called for its side effect of assigning a MVN sampler object
#   to environment mc, as well as an update method and possibly mat_sum and kron_prod methods
sparse_template <- function(mc, update.XX=FALSE, add.outer.R=FALSE) {

  if (mc[["type"]] == "block") {
    Q <- mc$QT
    keep.kp <- FALSE
  } else {
    # build a Kronecker product template
    if (mc$var == "unstructured") {
      Qv <- rWishart(1L, df=mc$prior$df, Sigma=diag(mc$q0))[,,1L]
    } else if (mc$var == "scalar" && !is.null(mc$Q0)) {
      Qv <- scale_mat(mc$Q0, runif(if (mc$usePX && mc$PX$vector) mc$q0 else 1L, 0.25, 0.75))
    } else if (mc$var == "diagonal" || (mc$usePX && mc$PX$vector)) {
      Qv <- runif(mc$q0, 0.25, 0.75)
    } else {  # scalar var, scalar or no PX
      Qv <- runif(1L, 0.25, 0.75)
    }
    if (mc$gl) {
      QA <- mc$glp$QA.ext
    } else if (mc$Leroux_update) {
      QA <- mc$mat_sum_Leroux(mc$QA, mc$idL, w2=runif(1L, 0.25, 0.75))
    } else if (is.null(mc$priorA)) {
      QA <- mc$QA
    } else {
      QA <- crossprod_sym(mc$DA, runif(mc$lD, 0.75, 1.25))
    }

    kron_prod <- build_kron(QA, Qv, q2=mc$q0, M1.fixed = is.null(mc$priorA) && !mc$Leroux_update)
    Q <- kron_prod(QA, Qv)
    if (mc$in_block) assign("Q", Q, envir=mc)  # only for single-time use in create_mc_block
    keep.kp <- TRUE
  }

  if (mc[["type"]] == "block" || !mc$in_block) {
    if (mc[["type"]] == "gen" && mc$gl) {
      XX <- mc$glp$XX.ext
      R <- mc$glp$R
    } else {
      XX <- mc$XX
      R <- mc$R
    }

    if (mc$e$control$expanded.cMVN.sampler) {
      if (update.XX) {
        mat_sum <- make_mat_sum(M1=XX, M2=Q)
        XX_Q <- mat_sum(XX, Q)
      } else {
        mat_sum <- make_mat_sum(M0=XX, M1=Q)
        XX_Q <- mat_sum(Q)
      }
      mc$MVNsampler <- create_cMVN_sampler(
        mbs=mc[["mcs"]], X=mc[["X"]], Q=XX_Q,
        R=R, r=mc$r,
        sampler=mc$e,
        name=mc$name, chol.control=mc$e$control$chol.control
      )
    } else {
      mc$MVNsampler <- NULL
      R0 <- R
      if (is.null(add.outer.R)) {
        # unresolved singularities usually only occur when at least two components are improper
        add.outer.R <- mc[["type"]] == "block" && sum(!sapply(mc$mcs, is_proper_GMRF)) >= 2L
        if (add.outer.R && !is.null(R)) {
          # it should not always be necessary to use the full matrix R
          # leave out the largest of the individual model components' contribution:
          cols <- sapply(mc$mcs, function(x) if (is.null(x$R)) 0L else ncol(x$R))
          cols <- cols[cols != 0L]
          if (length(cols) > 1L) {
            inds <- c(0L, cumsum(cols))
            maxR <- which.max(cols)
            R0 <- R[, -((inds[maxR] + 1L):inds[maxR + 1L]), drop=FALSE]
          }
        }
      }
      if (!add.outer.R || is.null(R)) {
        if (update.XX) {
          mat_sum <- make_mat_sum(M1=XX, M2=Q)
          XX_Q <- mat_sum(XX, Q)
        } else {
          mat_sum <- make_mat_sum(M0=XX, M1=Q)
          XX_Q <- mat_sum(Q)
        }
        tryCatch(
          suppressWarnings(
            mc$MVNsampler <- create_TMVN_sampler(
              Q=XX_Q, update.Q=TRUE, name=mc$name,
              R=R, r=mc$r, S=mc$S, s=mc$s,
              chol.control=mc$e$control$chol.control
            )
          ),
          error = function(e) {
            #if (grepl("Cholesky", e$message)) {  # TODO are there other possible messages for singular matrix conditions?
            # the problem is most likely a non-positive-definite XX + Q because of unresolved
            #   singularities in XX (due to levels without observations, or blocking and noninformative regression prior)
            # try to correct this by adding a multiple of tcrossprod(R) to XX (Rue and Held, 2005)
            #   the effect of which is cancelled by imposing these restrictions (R) during sampling
            # unfortunately, this typically leads to less sparse XX, XX+Q, and therefore slower Cholesky, or even out-of-memory
            # TODO add tcrossprod(R0) where R0 is the minimal selection of columns s.t. XX + R0 R0' is positive definite
            #   i.e. restrict R0 to a selection of constraints that involve all non-observed levels
            #   and perhaps use weights to further reduce the condition number
            if (is.null(R)) {
              if (mc[["type"]] == "block") {
                noninf <- sapply(mc$mcs, function(x) isFALSE(x$informative.prior))
                if (any(noninf))
                  warn("error in setting up a multivariate normal sampler, possibly due to ",
                       "a singular precision matrix; you may try to replace the non-informative prior ",
                       "of model component(s) '", paste(names(noninf)[noninf], collapse="', '"),
                       "' by a (weakly) informative prior")
              }
              stop(e)
            } else {
              if (class(XX_Q)[1L] == "dsCMatrix" && XX_Q@Dim[1L] >= 1000L) {
                message("adding outer product of constraint matrix to the ",
                        "precision matrix to improve its condition number; note that this may be ",
                        "inefficient if the resulting matrix is much less sparse")
              }
            }
          }
        )
      }
      if (is.null(mc$MVNsampler)) {  # 2nd attempt or add.outer.R
        if (is.null(R)) stop("singular precision matrix and no constraints")
        if (update.XX) {
          mat_sum <- make_mat_sum(M0=tcrossprod(R0), M1=XX, M2=Q)
          XX_Q <- mat_sum(XX, Q)
        } else {
          mat_sum <- make_mat_sum(M0=economizeMatrix(XX + tcrossprod(R0), symmetric=TRUE), M1=Q)
          XX_Q <- mat_sum(Q)
        }
        tryCatch(
          suppressWarnings(
            mc$MVNsampler <- create_TMVN_sampler(
              Q=XX_Q, update.Q=TRUE, name=mc$name,
              R=R, r=mc$r, S=mc$S, s=mc$s,
              chol.control=mc$e$control$chol.control
            )
          ),
          error = function(e) {
            if (mc[["type"]] == "block") {
              noninf <- sapply(mc$mcs, function(x) isFALSE(x$informative.prior) ||
                                 (x[["type"]] == "gen" && x$gl && isFALSE(x$glp$informative.prior))
              )
              if (any(noninf))
                warn("error in setting up a multivariate normal sampler, possibly due to ",
                     "a singular precision matrix; you may try to replace the non-informative prior ",
                     "of model component(s) '", paste(names(noninf)[noninf], collapse="', '"),
                     "' by a (weakly) informative prior")
            }
            stop(e)
          }
        )
      }
    }

    keep.ms <- TRUE
    if (mc[["type"]] == "gen" && mc$unit_Q) {  # unit QA, scalar Qv (with unit Q0), and no constraints
      if (isUnitDiag(XX)) {
        # TODO cleaner and more efficient handling of unit XX case (unit_Q + unit XX --> scalar update)
        XX.expanded <- expandUnitDiag(XX)
        mc$update <- function(XX, QA, Qv, w) list(Q=XX.expanded, Imult=Qv*w)
        keep.kp <- keep.ms <- FALSE
      } else {
        if (class(XX)[1L] == class(if (update.XX) mat_sum(XX, Q) else mat_sum(Q))[1L]) {
          keep.kp <- keep.ms <- FALSE
          mc$update <- function(XX, QA, Qv, w) list(Q=XX, Imult=Qv*w)
        } else {  # keep mat_sum, kron_prod because for some reason XX has type incompatible with ch
          if (update.XX)
            mc$update <- function(XX, QA, Qv, w) list(Q=mat_sum(XX, kron_prod(QA, Qv * w)), Imult=0)
          else
            mc$update <- function(XX, QA, Qv, w) list(Q=mat_sum(kron_prod(QA, Qv * w)), Imult=0)
        }
      }
    } else {
      if (mc[["type"]] == "block" || mc$q0 == 1L) {  # scalar Qv, or block
        if (mc[["type"]] != "block") keep.kp <- FALSE
        if (update.XX)
          mc$update <- function(XX, QA, Qv, w) list(Q=mat_sum(XX, QA, w2=Qv*w), Imult=0)
        else
          mc$update <- function(XX, QA, Qv, w) list(Q=mat_sum(QA, w1=Qv*w), Imult=0)
      } else {  # q0 > 1, no block
        if (update.XX)
          mc$update <- function(XX, QA, Qv, w) list(Q=mat_sum(XX, kron_prod(QA, Qv), w2=w), Imult=0)
        else
          mc$update <- function(XX, QA, Qv, w) list(Q=mat_sum(kron_prod(QA, Qv), w1=w), Imult=0)
      }
    }
    if (keep.ms) mc$mat_sum <- mat_sum
  }

  if (mc[["type"]] != "block" && mc$in_block) keep.kp <- TRUE
  if (keep.kp) mc$kron_prod <- kron_prod

}
