#' Functions for specifying a sampling distribution and link function
#'
#' These functions are intended for use in the \code{family} argument of \code{\link{create_sampler}}.
#' In future versions these functions may gain additional arguments, but currently the corresponding
#' functions \code{gaussian} and \code{binomial} can be used as well.
#'
#' @param link the name of a link function. Currently the only allowed link functions are:
#'  \code{"identity"} for (log-)Gaussian sampling distributions, \code{"logit"} (default) and \code{"probit"}
#'  for binomial distributions and \code{"log"} for negative binomial sampling distributions.
#' @return A family object.
#' @name mcmcsae-family
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
