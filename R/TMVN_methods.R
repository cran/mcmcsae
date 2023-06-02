#' Functions for specifying the method and corresponding options for sampling
#' from a possibly truncated and degenerate multivariate normal distribution
#'
#' These functions are intended for use in the \code{method} argument of \code{\link{create_TMVN_sampler}}.
#' In future versions these functions may gain additional arguments.
#'
#' @param diagnostic whether information about violations of inequalities and bounces
#'  off inequality walls (for 'HMC' and 'HMCZigZag' methods) or gradient events
#'  (for 'HMCZigZag') is printed to the screen.
#' @param slice if \code{TRUE}, a Gibbs within slice sampler is used.
#' @param Tsim the duration of a Hamiltonian Monte Carlo simulated particle trajectory.
#'  This can be specified as either a single positive numeric value for a fixed
#'  simulation time, or as a function that is applied in each MCMC iteration to
#'  generates a simulation time.
#' @param max.events maximum number of events (reflections off inequality walls
#'  and for method 'HMCZigZag' also gradient events). Default is unlimited.
#'  Specifying a finite number may speed up the sampling but may also result
#'  in a biased sampling algorithm.
#' @param rate vector of Laplace rate parameters for method 'HMCZigZag'. It must be
#'  a positive numeric vector of length one or the number of variables.
#' @param prec.eq positive numeric vector of length 1 or the number of equality restrictions,
#'  to control the precision with which the equality restrictions are imposed; the larger
#'  \code{prec.eq} the more precisely they will be imposed.
#' @param adapt experimental feature: if \code{TRUE} the rate parameter will be adapted
#'  in an attempt to make the sampling algorithm more efficient.
#' @param sharpness for method 'softTMVN', the sharpness of the soft inequalities; the larger the better
#'  the approximation of exact inequalities. It must a positive numeric vector of length
#'  one or the number of inequality restrictions.
#' @param useV for method 'softTMVN' whether to base computations on variance instead of precision
#'  matrices.
#' @param debug if \code{TRUE} a breakpoint is set at the beginning of the TMVN sampling
#'  function. Mainly intended for developers.
#' @return A method object, for internal use only.
#' @name mcmcsae-TMVN-method
NULL


#' @export
#' @rdname mcmcsae-TMVN-method
m_direct <- function() {
  list(method = "direct")
}

#' @export
#' @param diagnostic whether information about violations of inequalities is printed to the screen
#' @rdname mcmcsae-TMVN-method
m_Gibbs <- function(slice=FALSE, diagnostic=FALSE, debug=FALSE) {
  list(method = "Gibbs", slice=slice, diagnostic=diagnostic, debug=debug)
}

#' @export
# @param diagnostic whether information about violations of inequalities and bounces off inequality walls
#  is printed to the screen. This sometimes provides useful diagnostic information.
#' @rdname mcmcsae-TMVN-method
m_HMC <- function(Tsim=pi/2, max.events=.Machine$integer.max, diagnostic=FALSE, debug=FALSE) {
  if (is.function(Tsim)) {
    Tsim_value <- Tsim()
  } else {
    Tsim_value <- as.numeric(Tsim)[1L]
    Tsim <- function() Tsim_value
  }
  if (!is.numeric(Tsim_value) || length(Tsim_value) != 1L || Tsim_value <= 0) stop("'Tsim' should be a single positive numeric value, or a function to generate one")
  max.event <- as.integer(max.events)
  if (length(max.events) != 1L || max.events <= 0) stop("'max.events must be a positive integer'")
  list(method = "HMC", Tsim=Tsim, max.events=max.events, diagnostic=diagnostic, debug=debug)
}

#' @export
#' @rdname mcmcsae-TMVN-method
m_HMCZigZag <- function(Tsim=1, rate=1, prec.eq=NULL, diagnostic=FALSE,
                        max.events=.Machine$integer.max, adapt=FALSE, debug=FALSE) {
  if (is.function(Tsim)) {
    Tsim_value <- Tsim()
  } else {
    Tsim_value <- as.numeric(Tsim)[1L]
    Tsim <- function() Tsim_value
  }
  if (!is.numeric(Tsim_value) || length(Tsim_value) != 1L || Tsim_value <= 0) stop("'Tsim' should be a single positive numeric value, or a function to generate one")
  if (any(rate <= 0)) stop("Laplace distribution scale parameters must be positive")
  max.event <- as.integer(max.events)
  if (length(max.events) != 1L || max.events <= 0) stop("'max.events must be a positive integer'")
  list(method = "HMCZigZag", Tsim=Tsim, rate=rate, prec.eq=prec.eq, diagnostic=diagnostic,
       max.events=max.events, adapt=adapt, debug=debug)
}

#' @export
#' @rdname mcmcsae-TMVN-method
m_softTMVN <- function(sharpness=100, useV=FALSE, debug=FALSE) {
  list(method = "softTMVN", sharpness=sharpness, useV=useV, debug=debug)
}
