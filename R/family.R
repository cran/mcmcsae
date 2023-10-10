#' Functions for specifying a sampling distribution and link function
#'
#' These functions are intended for use in the \code{family} argument of \code{\link{create_sampler}}.
#' In future versions these functions may gain additional arguments, but currently the corresponding
#' functions \code{gaussian} and \code{binomial} can be used as well.
#'
#' @param link the name of a link function. Currently the only allowed link functions are:
#'  \code{"identity"} for (log-)Gaussian sampling distributions, \code{"logit"} (default) and \code{"probit"}
#'  for binomial distributions and \code{"log"} for negative binomial sampling distributions.
#' @param K number of categories for multinomial model; this must be specified for prior predictive sampling.
#' @param shape.vec optional formula specification of unequal shape parameter for gamma family
#' @param shape.prior prior for gamma shape parameter. Supported prior distributions:
#'  \code{\link{pr_fixed}} with a default value of 1, \code{\link{pr_exp}} and
#'  \code{\link{pr_gamma}}. The current default is that of a fixed shape
#'  equal to 1, i.e. \code{pr_fixed(value=1)}.
#' @param shape.MH.type the type of Metropolis-Hastings algorithm employed
#'  in case the shape parameter is to be inferred. The two choices currently
#'  supported are "RW" for a random walk proposal on the log-shape scale
#'  and "gamma" for an approximating gamma proposal, found using an iterative
#'  algorithm. In the latter case, a Metropolis-Hastings accept-reject step is
#'  currently omitted, so the sampling algorithm is an approximate one,
#'  though one that is usually quite accurate and efficient.
#' @return A family object.
#' @name mcmcsae-family
#' @references
#'  J.W. Miller (2019).
#'    Fast and Accurate Approximation of the Full Conditional for Gamma Shape Parameters.
#'    Journal of Computational and Graphical Statistics 28(2), 476-480.
NULL

#' @export
#' @rdname mcmcsae-family
f_gaussian <- function(link="identity") {
  link <- match.arg(link)
  list(family="gaussian", link=link, linkinv=identity)
}

#' @export
#' @rdname mcmcsae-family
f_binomial <- function(link=c("logit", "probit")) {
  link <- match.arg(link)
  if (link == "logit")
    linkinv <- make.link(link)$linkinv
  else
    linkinv <- pnorm
  list(family="binomial", link=link, linkinv=linkinv)
}

#' @export
#' @rdname mcmcsae-family
f_negbinomial <- function(link="logit") {
  link <- match.arg(link)
  list(family="negbinomial", link=link, linkinv=make.link(link)$linkinv)
}

#' @export
#' @rdname mcmcsae-family
f_poisson <- function(link="log") {
  link <- match.arg(link)
  list(family="poisson", link=link, linkinv=make.link(link)$linkinv)
}

#' @export
#' @rdname mcmcsae-family
f_multinomial <- function(link="logit", K=NULL) {
  link <- match.arg(link)
  if (!is.null(K)) {
    K <- as.integer(K)
    if (length(K) != 1L) stop("number of categories 'K' must be a scalar integer")
    if (K < 2L) stop("number of categories 'K' must be at least 2")
  }
  list(family="multinomial", link=link, linkinv=make.link(link)$linkinv, K=K)
}

#' @export
#' @rdname mcmcsae-family
# TODO use control function for shape.vec, shape.prior and shape.MH.type
f_gamma <- function(link="log", shape.vec = ~ 1, shape.prior = pr_fixed(1),
                    shape.MH.type = c("RW", "gamma")) {
  family <- "gamma"  # later change to name <- "gamma"
  link <- match.arg(link)
  linkinv=make.link(link)$linkinv
  if (!inherits(shape.vec, "formula")) stop("'shape.vec' must be a formula")
  if (shape.prior$type == "fixed") {
    if (shape.prior$value <= 0) stop("gamma shape parameter must be positive")
  } else if (shape.prior$type == "exp") {  # special case of gamma
    shape.prior <- pr_gamma(shape=1, rate=1/shape.prior$scale)
  } else if (shape.prior$type != "gamma") {
    stop("unsupported prior for gamma shape parameter")
  }
  alpha.fixed <- shape.prior$type == "fixed"
  shape.prior$init(1L)  # scalar parameter
  alpha.scalar <- intercept_only(shape.vec)
  shape.MH.type <- match.arg(shape.MH.type)
  alpha0 <- NULL
  get_shape <- function() stop("please call method 'init' first")
  init <- function(data) get_shape <<- make_get_shape(data)
  make_get_shape <- function(data) {
    if (alpha.scalar) {
      alpha0 <<- 1
    } else {
      alpha0 <<- model_matrix(update.formula(shape.vec, ~ . - 1), data)
      if (ncol(alpha0) != 1L) stop("'shape.vec' must contain a single numeric variable")
      alpha0 <<- as.numeric(alpha0)
    }
    if (alpha.fixed) {
      alpha0 <<- alpha0 * shape.prior$value
      function(p) alpha0
    } else {
      if (alpha.scalar)
        function(p) p[["gamma_shape_"]]
      else
        function(p) alpha0 * p[["gamma_shape_"]]
    }
  }
  if (!alpha.fixed) {
    make_draw_shape <- function(y) {
      # set up sampler for full conditional posterior for alpha, given linear predictor
      n <- length(y)
      #if (shape.prior$type == "fixed") return(function(p) shape.prior$value)
      # MH within Gibbs
      switch(shape.MH.type,
        RW = {
          # random walk update of log(alpha), with normal stdev tau
          tau <- 0.2  # start value of RW stdev
          adapt <- function(ar) {
            if (ar[["gamma_shape_"]] < .2)
              tau <<- tau * runif(1L, 0.6, 0.9)
            else if (ar[["gamma_shape_"]] > .7)
              tau <<- tau * runif(1L, 1.1, 1.5)
          }
          f <- function(p) {
            alpha <- p[["gamma_shape_"]]
            alpha.star <- alpha * exp(rnorm(1L, sd=tau))
          }
          if (alpha.scalar) {
            sumlogy <- sum(log(y))
            f <- add(f, quote(
              log.ar <- shape.prior[["shape"]] * log(alpha.star/alpha) - shape.prior[["rate"]] * (alpha.star - alpha) +
                n * (lgamma(alpha) - lgamma(alpha.star) + alpha.star * log(alpha.star) - alpha * log(alpha)) +
                (alpha.star - alpha) * (sumlogy - sum(p[["e_"]]) - sum(y * exp(-p[["e_"]])))
              ))
          } else {
            f <- add(f, quote(alpha.vec <- alpha0 * alpha))
            f <- add(f, quote(alpha.star.vec <- alpha0 * alpha.star))
            f <- add(f, quote(
              log.ar <- shape.prior[["shape"]] * log(alpha.star/alpha) - shape.prior[["rate"]] * (alpha.star - alpha) +
                sum(lgamma(alpha.vec) - lgamma(alpha.star.vec)) +
                sum(alpha.star.vec * log(alpha.star.vec * y)) -
                sum(alpha.vec * log(alpha.vec * y)) +
                sum((alpha.vec - alpha.star.vec) * (y * exp(-p[["e_"]]) + p[["e_"]]))
            ))
          }
          add(f, quote(if (log(runif(1L)) < log.ar) alpha.star else alpha))
        },
        gamma = {
          # Miller's gamma approximation of the shape's full conditional
          #   iteration starts with approximate gamma density derived using Stirling's formula
          if (alpha.scalar) {
            A0 <- shape.prior[["shape"]] + 0.5 * n
            B00 <- shape.prior[["rate"]] - sum(log(y)) - n
            function(p) {
              A <- A0
              B0 <- B00 + sum(y * exp(-p[["e_"]])) + sum(p[["e_"]])
              B <- B0
              #A <- shape.prior$shape + n
              #B0 <- shape.prior$rate + sum(y * exp(-p[["e_"]])) - sumlogy + sum(p[["e_"]]) - n
              #B <- B0 + n
              for (i in 1:10) {
                a <- A/B
                A <- shape.prior[["shape"]] + n*a*(a * trigamma(a) - 1)
                B <- B0 + (A - shape.prior[["shape"]])/a + n*(digamma(a) - log(a))
                if (abs(a/(A/B) - 1) < 1e-8) break
              }
              rgamma(1L, A, B)
            }
            # To add an MH correction step:
            # (does not seem necessary, as the approximation is often excellent)
            #alpha.star <- rgamma(1L, A, B)
            #alpha <- p[["gamma_shape_"]]
            #log.ar <- n * (lgamma(alpha) - lgamma(alpha.star) + alpha.star * log(alpha.star) - alpha * log(alpha)) +
            #  (alpha.star - alpha) * (sumlogy - sum(p[["e_"]]) - sum(y * exp(-p[["e_"]]))) +
            #  (A - shape.prior$shape) * log(alpha/alpha.star) - (B - shape.prior$rate) * (alpha - alpha.star)
            #if (log(runif(1L)) < log.ar) alpha.star else alpha
          } else {
            ala0 <- sum(alpha0 * log(alpha0))
            A0 <- shape.prior[["shape"]] + 0.5 * n
            B00 <- shape.prior[["rate"]] - sum(alpha0 * (log(y) + 1))
            function(p) {
              A <- A0
              B0 <- B00 + sum(alpha0 * (y * exp(-p[["e_"]]) + p[["e_"]]))
              B <- B0
              for (i in 1:10) {
                a <- A/B
                a.vec <- a * alpha0
                A <- shape.prior[["shape"]] - sum(a.vec) + sum(a.vec * a.vec * trigamma(a.vec))
                B <- B0 + (A - shape.prior[["shape"]])/a - log(a) * sum(alpha0) +
                  sum(alpha0 * digamma(a.vec)) - ala0
                if (abs(a/(A/B) - 1) < 1e-8) break
              }
              rgamma(1L, A, B)
            }
            # possibly add MH correction step
          }
        }
      )
    }
  }
  make_llh <- function(y) {
    n <- length(y)
    if (alpha.fixed) {
      alpha <- get_shape()
      if (alpha.scalar) {
        llh_0 <- (alpha - 1) * sum(log(y))
        llh_0 <- llh_0 + n * (alpha * log(alpha) - lgamma(alpha))
        function(p) llh_0 - alpha * sum(p[["e_"]] + y * exp(-p[["e_"]]))
      } else {
        llh_0 <- sum((alpha - 1) * log(y))
        llh_0 <- llh_0 + sum(alpha * log(alpha) - lgamma(alpha))
        function(p) llh_0 - sum(alpha * (p[["e_"]] + y * exp(-p[["e_"]])))
      }
    } else {
      llh_0 <- -sum(log(y))
      if (alpha.scalar) {
        function(p) {
          alpha <- get_shape(p)
          (1 - alpha) * llh_0 + n * (alpha * log(alpha) - lgamma(alpha)) - alpha * sum(p[["e_"]] + y * exp(-p[["e_"]]))
        }
      } else
        function(p) {
          alpha <- get_shape(p)
          llh_0 + sum(alpha * log(alpha * y) - lgamma(alpha)) - sum(alpha * (p[["e_"]] + y * exp(-p[["e_"]])))
        }
    }
  }
  make_llh_i <- function(y) {
    n <- length(y)
    if (alpha.fixed) {
      alpha <- get_shape()
      pllh_0 <- alpha * log(alpha) - lgamma(alpha) + (alpha - 1) * log(y)
      if (alpha.scalar)
        function(draws, i, e_i) {
          nr <- dim(e_i)[1L]
          rep_each(pllh_0[i], nr) - alpha * (e_i + rep_each(y[i], nr) * exp(-e_i))
        }
      else
        function(draws, i, e_i) {
          nr <- dim(e_i)[1L]
          rep_each(pllh_0[i], nr) - rep_each(alpha[i], nr) * (e_i + rep_each(y[i], nr) * exp(-e_i))
        }
    } else {
      if (alpha.scalar)
        function(draws, i, e_i) {
          nr <- dim(e_i)[1L]
          alpha <- as.numeric(as.matrix.dc(draws[["gamma_shape_"]], colnames=FALSE))
          alpha * log(alpha) - lgamma(alpha) + (alpha - 1) * rep_each(log(y[i]), nr) +
            - alpha * (e_i + rep_each(y[i], nr) * exp(-e_i))
        }
      else
        function(draws, i, e_i) {
          nr <- dim(e_i)[1L]
          alpha_i <- outer(as.numeric(as.matrix.dc(draws[["gamma_shape_"]], colnames=FALSE)), alpha0[i])
          alpha_i * log(alpha_i) - lgamma(alpha_i) + (alpha_i - 1) * rep_each(log(y[i]), nr) +
            - alpha_i * (e_i + rep_each(y[i], nr) * exp(-e_i))
        }
    }
  }
  make_rpredictive <- function(newdata) {
    if (is.null(newdata)) {
      # in-sample prediction/replication, linear predictor,
      # or custom X case, see prediction.R
      rpredictive <- function(p, lp) {
        alpha <- get_shape(p)
        rgamma(length(lp), shape=alpha, rate=alpha * exp(-lp))
      }
    } else {
      newn <- nrow(newdata)
      get_newshape <- make_get_shape(newdata)
      rpredictive <- function(p, lp) {
        alpha <- get_newshape(p)
        rgamma(newn, shape=alpha, rate=alpha * exp(-lp))
      }
    }
  }
  environment()
}
