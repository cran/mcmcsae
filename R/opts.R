
.opts <- new.env()

#' Set global options relating to computational details
#'
#' @export
#' @param auto.order.block whether Gibbs blocks should be ordered automatically in such a
#'  way that those with the most sparse design matrices come first. This way of ordering
#'  can make Cholesky updates more efficient.
#' @param chol.inplace whether sparse Cholesky updates should re-use the same memory location.
#' @param chol.ordering an integer passed to CHOLMOD routines determining which reordering
#'  schemes are tried to limit sparse Cholesky fill-in.
#' @param PG.approx whether Polya-Gamma draws for logistic binomial models are
#'  approximated by a hybrid gamma convolution approach. If not, \code{BayesLogit::rpg}
#'  is used, which is exact for some values of the shape parameter.
#' @param PG.approx.m if \code{PG.approx=TRUE}, the number of explicit gamma draws in the
#'  sum-of-gammas representation of the Polya-Gamma distribution. The remainder (infinite)
#'  convolution is approximated by a single moment-matching gamma draw. Special values are:
#'  \code{-2L} for a default choice depending on the value of the shape parameter and
#'  balancing performance and accuracy, \code{-1L} for a moment-matching normal approximation,
#'  and \code{0L} for a moment-matching gamma approximation.
#' @param CRT.approx.m scalar integer specifying the degree of approximation to sampling
#'  from a Chinese Restaurant Table distribution. The approximation is based on Le Cam's theorem.
#'  Larger values yield a slower but more accurate sampler.
#' @return This function sets or resets options in the option environment \code{.opts}.
#' @references
#'  D. Bates, M. Maechler, B. Bolker and S.C. Walker (2015).
#'    Fitting Linear Mixed-Effects Models Using lme4.
#'    Journal of Statistical Software 67(1), 1-48.
#'
#'  Y. Chen, T.A. Davis, W.W. Hager and S. Rajamanickam (2008).
#'    Algorithm 887: CHOLMOD, supernodal sparse Cholesky factorization and update/downdate.
#'    ACM Transactions on Mathematical Software 35(3), 1-14.
set_opts <- function(auto.order.block=TRUE, chol.inplace=TRUE, chol.ordering=0L,
                     PG.approx=TRUE, PG.approx.m=-2L, CRT.approx.m=20L) {
  .opts$auto.order.block <- auto.order.block
  .opts$chol.inplace <- chol.inplace
  .opts$chol.ordering <- as.integer(chol.ordering)
  .opts$PG.approx <- PG.approx
  .opts$PG.approx.m <- as.integer(PG.approx.m)
  .opts$CRT.approx.m <- as.integer(CRT.approx.m)
}
set_opts()


check_ny <- function(ny, n, data) {
  if (is.null(ny)) return(1L)  # Bernoulli-logistic model
  if (is.character(ny)) {
    if (length(ny) != 1L) stop("wrong input for 'ny'")
    ny <- eval_in(ny, data)
  } else {
    if (!is.numeric(ny)) stop("wrong input for 'ny'")
    if (!(length(ny) %in% c(1L, n))) stop("wrong length for 'ny'")
    if (anyNA(ny)) stop("missings in 'ny' not allowed")
    if (any(ny < 0)) stop("'ny' cannot be negative")
  }
  if (any(ny - as.integer(ny) > sqrt(.Machine$double.eps))) warning("non-integral values in 'ny' have been rounded")
  as.integer(round(ny))
  # or should we allow non-integral ny? allowed for model fitting but not for prediction by rbinom
  # note that y in create_sampler is currently not rounded, but y <= ny is checked
  # maybe add argument round and round only in case of prediction
}
