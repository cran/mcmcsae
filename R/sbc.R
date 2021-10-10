
#' Simulation based calibration
#'
#' @examples
#' \dontrun{
#' # this example takes a long time
#' n <- 10L  # data size
#' dat <- data.frame(x=runif(n))
#' SBC_test(~ reg(~ 1 + x, b0=c(0.25, 1), Q0=1, name="beta"),
#'   sigma.mod=pr_invchisq(df=1, scale=list(df=1, scale=1)), data=dat,
#'   pars=list(mu="beta[1]", beta_x="beta[2]", sigma="sigma_"),
#'   n.draws=25L, n.sim=100L*25L, thin=3L, burnin=50L
#' )
#' }
#'
#' @export
#' @keywords internal
#' @param ... passed to \code{\link{create_sampler}} (can be all parameters except \code{prior.only})
#' @param pars named list with univariate functions of the parameters to use in test. This list
#'   is passed to argument \code{pred} of \code{\link{MCMCsim}}.
#' @param n.draws number of posterior draws to retain in posterior simulations.
#' @param n.sim number of simulation iterations.
#' @param burnin burnin to use in posterior simulations, passed to \code{\link{MCMCsim}}.
#' @param thin thinning to use in posterior simulations, passed to \code{\link{MCMCsim}}.
#' @return A matrix with ranks.
#' @references
#'   S. Talts, M. Betancourt, D. Simpson, A. Vehtari and A. Gelman (2018).
#'     Validating Bayesian inference algorithms with simulation-based calibration.
#'     arXiv:1804.06788.
# TODO parallel option
SBC_test <- function(..., pars, n.draws=25L, n.sim=20L*n.draws, burnin=25L, thin=2L, show.progress=TRUE) {
  sampler_args <- list(...)
  sampler_args$prior.only <- TRUE
  sampler <- do.call("create_sampler", sampler_args)
  #sampler <- create_sampler(..., prior.only=TRUE)
  n.pars <- length(pars)
  ranks <- matrix(0, n.draws + 1L, n.pars, dimnames=list(NULL, names(pars)))
  # TODO: in some cases like TMVN prior draws may not be iid; in that case use burnin and thin here as well
  sim0 <- MCMCsim(sampler, burnin=0L, n.iter=n.sim, n.chain=1L, pred=pars, from.prior=TRUE, verbose=FALSE)
  pred <- predict(sim0, show.progress=FALSE)
  sampler_args$prior.only <- FALSE
  if (names(sampler_args)[1L] == "") names(sampler_args)[1L] <- "formula"
  sampler_args$formula <- update.formula(sampler_args$formula, y.tilde ~ .)
  environment(sampler_args$formula) <- environment()
  if (show.progress) pb <- txtProgressBar(min=1L, max=n.sim, style=3L)
  for (i in seq_len(n.sim)) {
    y.tilde <- pred[[1L]][i, ]
    # TODO function to update sampler to use y.tilde
    sampler <- do.call("create_sampler", sampler_args)
    sim <- MCMCsim(sampler, pred=pars, burnin=burnin, n.chain=1L, n.iter=thin*n.draws, thin=thin, verbose=FALSE)
    for (v in seq_along(pars)) {
      vname <- names(pars)[v]
      r <- sum(sim[[vname]][[1L]] < sim0[[vname]][[1L]][i, ]) + 1L
      ranks[r, v] <- ranks[r, v] + 1L
    }
    if (show.progress) setTxtProgressBar(pb, i)
  }
  if (show.progress) close(pb)
  ranks
}
