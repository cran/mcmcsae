
context("Poisson regression")

set.seed(1, kind="Mersenne-Twister", normal.kind="Inversion")

# generate Poisson data
n <- 1000
x <- rnorm(n)
eta <- 1 + 2*x
y <- rpois(n, lambda = exp(eta))

test_that("negative binomial approximation to Poisson regression works", {
  r <- 100  # should be large for negligible overdispersion, but too large values result in slow MCMC mixing
  o <- rep(-log(r), n)
  sampler <- create_sampler(y ~ 1 + x + offset(o), ry=r, r.mod=pr_fixed(1), family="negbinomial")
  sim <- MCMCsim(sampler, n.chain=2, burnin=150, n.iter=500, verbose=FALSE)
  summ <- summary(sim)
  expect_equal(summ$reg1[, "Mean"], c(`(Intercept)`=1, x=2), tolerance=0.25)
})

test_that("Poisson shortcut works", {
  sampler <- create_sampler(y ~ 1 + x, family="poisson")
  sim <- MCMCsim(sampler, n.chain=2, burnin=150, n.iter=500, verbose=FALSE)
  summ <- summary(sim)
  expect_equal(summ$reg1[, "Mean"], c(`(Intercept)`=1, x=2), tolerance=0.25)
})
