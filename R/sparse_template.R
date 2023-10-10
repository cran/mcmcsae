# mc: model component
# update.XX: whether XX is updated in each draw
# add.outer.R: whether outer product of constraints should be added (for positive definiteness in blocked components)
sparse_template <- function(mc, update.XX=FALSE, add.outer.R=FALSE, prior.only=FALSE) {

  if (mc$type == "block") {
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

  if (mc$type == "block" || !mc$in_block) {
    if (mc$type == "gen" && mc$gl) {
      XX <- mc$glp$XX.ext
      R <- mc$glp$R
    } else {
      XX <- mc$XX
      R <- mc$R
    }

    mc$MVNsampler <- NULL
    if (!add.outer.R || is.null(R)) {
      if (update.XX)
        mat_sum <- make_mat_sum(M1=XX, M2=Q)
      else
        mat_sum <- make_mat_sum(M0=XX, M1=Q)
      XX_Q <- if (update.XX) mat_sum(XX, Q) else mat_sum(Q)
      tryCatch(
        suppressWarnings(
          mc$MVNsampler <- create_TMVN_sampler(Q=XX_Q, update.Q=TRUE, name=mc$name, R=R, S=mc$S, chol.control=mc$e$control$chol.control)
        ),
        error = function(e) {
          # TODO really need this? examples? besides, determinant check is not a very reliable test for singularity
          if (determinant(XX, logarithm=FALSE)$modulus < .Machine$double.eps) {
            # the problem is most likely a non-positive-definite XX + Q because of unresolved
            #   singularities in XX (due to levels without observations, or blocking)
            # try to correct this by adding a multiple of tcrossprod(R) to XX (Rue and Held, 2005)
            #   the effect of which is cancelled by imposing these restrictions (R) during sampling
            # unfortunately, this typically leads to less sparse XX, XX+Q, and therefore slower Cholesky, or even out-of-memory
            # TODO add tcrossprod(R0) where R0 is the minimal selection of columns s.t. XX + R0 R0' is positive definite
            #   i.e. restrict R0 to a selection of constraints that involve all non-observed levels
            if (is.null(R)) stop("singular precision matrix and no constraints")
          } else {
            e
          }
        }
      )
    }
    if (is.null(mc$MVNsampler)) {  # 2nd attempt or add.outer.R
      if (is.null(R)) stop("singular precision matrix and no constraints")
      if (update.XX)
        mat_sum <- make_mat_sum(M0=tcrossprod(R), M1=XX, M2=Q)
      else
        mat_sum <- make_mat_sum(M0=economizeMatrix(XX + tcrossprod(R), symmetric=TRUE), M1=Q)
      XX_Q <- if (update.XX) mat_sum(XX, Q) else mat_sum(Q)
      mc$MVNsampler <- create_TMVN_sampler(Q=XX_Q, update.Q=TRUE, name=mc$name, R=R, S=mc$S, chol.control=mc$e$control$chol.control)
    }

    keep.ms <- TRUE
    if (mc$type == "gen" && mc$unit_Q) {  # unit QA, scalar Qv (with unit Q0), and no constraints
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
      if (mc$type == "block" || mc$q0 == 1L) {  # scalar Qv, or block
        if (mc$type != "block") keep.kp <- FALSE
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

  if (mc$type != "block" && mc$in_block) keep.kp <- TRUE
  if (keep.kp) mc$kron_prod <- kron_prod

}
