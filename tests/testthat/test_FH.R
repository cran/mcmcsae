
context("Fay-Herriot area-level model")

set.seed(1, kind="Mersenne-Twister", normal.kind="Inversion")

m <- 50L
df <- data.frame(area = seq_len(m), x = runif(m))
beta <- c(1, -0.5)
sigma_v <- 1
theta <- rnorm(m, beta[1] + beta[2]*df$x, sigma_v)
psi <- (1 + 2*runif(m))^2
df$y <- rnorm(m, theta, sqrt(psi))

test_that("MCMC for FH model runs", {
  sampler <- create_sampler(
    y ~ reg(~ 1 + x, name="beta") + gen(factor = ~ area, name="v"),
    Q0=1/psi, sigma.fixed=TRUE, data=df
  )
  sim <- MCMCsim(sampler, n.iter=500, burnin=100, n.chain=2, verbose=FALSE)
  summ <- summary(sim)
  compute_DIC(sim)
  expect_error(compute_WAIC(sim))  # for WAIC need store.all=TRUE
})

test_that("shortcut '_local' for area indicator runs and compute_WAIC runs", {
  sampler <- create_sampler(
    y ~ reg(~ 1 + x, name="beta") + gen(factor = ~ local_, name="v"),
    Q0=Diagonal(x=1/psi), sigma.fixed=TRUE, data=df
  )
  sim <- MCMCsim(sampler, n.iter=500, burnin=100, n.chain=2, verbose=FALSE, store.all=TRUE)
  expect_is(sim$v, "dc")
  summ <- summary(sim)
  expect_named(compute_DIC(sim))
  expect_named(compute_WAIC(sim))
})

test_that("single block Gibbs sampler runs", {
  sampler <- create_sampler(
    y ~ reg(~ 1 + x, name="beta") + gen(factor = ~ iid(area)),
    Q0=Cdiag(1/psi), sigma.fixed=TRUE,
    data=df, block=TRUE
  )
  sim <- MCMCsim(sampler, n.iter=500, burnin=100, n.chain=2, verbose=FALSE, store.all=TRUE)
  summ <- summary(sim)
  expect_is(summ$beta, "dc_summary")
  compute_DIC(sim)
  compute_WAIC(sim)
})

test_that("FH model with unit sampling variances runs", {
  # misspecified model with identity Q0
  sampler <- create_sampler(y ~ reg(~ 1 + x) + gen(factor = ~ area),
    sigma.fixed=TRUE, data=df
  )
  sim <- MCMCsim(sampler, n.iter=250, burnin=100, n.chain=2, verbose=FALSE)
  expect_is(sim, "draws")
})
