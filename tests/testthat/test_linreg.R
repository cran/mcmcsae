
context("Linear regression")

set.seed(1, kind="Mersenne-Twister", normal.kind="Inversion")

n <- 1000L
df <- data.frame(
  x1 = rnorm(n),
  x2 = runif(n),
  x3 = rbinom(n, 1, 0.1),
  x4 = rgamma(n, 1, 1)
)
df$y <- with(df, 1 + x1 + 2*x2 + 3*x3 + 4*x4 + rnorm(n))

test_that("linear regression works, as well as fitted and residuals methods", {
  sampler <- create_sampler(y ~ reg(~1+x1+x2+x3+x4, name="beta"), data=df)
  sim <- MCMCsim(sampler, n.iter=500, burnin=100, n.chain=2, verbose=FALSE)
  summ <- summary(sim)
  expect_true(all((0.5 * c(1,1,2,3,4) < summ$beta[, "Mean"]) & (summ$beta[, "Mean"] < 2 * c(1,1,2,3,4))))
  expect_true((0.5 < summ$sigma_[, "Mean"]) && (summ$sigma_[, "Mean"] < 2))
  expect_equal(fitted(sim, mean.only=TRUE), as.vector(summary(fitted(sim))[, "Mean"]))
  expect_equal(residuals(sim, mean.only=TRUE), as.vector(summary(residuals(sim))[, "Mean"]))
})

test_that("single-site Gibbs sampler works", {
  sampler <- create_sampler(
    y ~ reg(~1, name="mu") + reg(~x1-1, name="beta1") + reg(~x2-1, name="beta2") +
        reg(~x3-1, name="beta3") + reg(~x4-1, name="beta4"), data=df
  )
  sim <- MCMCsim(sampler, n.iter=500, burnin=100, n.chain=2, verbose=FALSE)
  summ <- summary(sim)
  expect_true((0.5 < summ$mu[, "Mean"]) && (summ$mu[, "Mean"] < 2))
  expect_true((0.5 < summ$beta1[, "Mean"]) && (summ$beta1[, "Mean"] < 2))
  expect_true((2*0.5 < summ$beta2[, "Mean"]) && (summ$beta2[, "Mean"] < 2*2))
  expect_true((3*0.5 < summ$beta3[, "Mean"]) && (summ$beta3[, "Mean"] < 3*2))
  expect_true((4*0.5 < summ$beta4[, "Mean"]) && (summ$beta4[, "Mean"] < 4*2))
  expect_true((0.5 < summ$sigma_[, "Mean"]) && (summ$sigma_[, "Mean"] < 2))
})

test_that("full blocking works", {
  sampler <- create_sampler(
    y ~ reg(~1, name="mu") + reg(~x1-1, name="beta1") + reg(~x2-1, name="beta2") +
      reg(~x3-1, name="beta3") + reg(~x4-1, name="beta4"),
    data=df, block=TRUE
  )
  sim <- MCMCsim(sampler, n.iter=500, burnin=100, n.chain=2, verbose=FALSE)
  summ <- summary(sim)
  expect_true((0.5 < summ$mu[, "Mean"]) && (summ$mu[, "Mean"] < 2))
  expect_true((0.5 < summ$beta1[, "Mean"]) && (summ$beta1[, "Mean"] < 2))
  expect_true((2*0.5 < summ$beta2[, "Mean"]) && (summ$beta2[, "Mean"] < 2*2))
  expect_true((3*0.5 < summ$beta3[, "Mean"]) && (summ$beta3[, "Mean"] < 3*2))
  expect_true((4*0.5 < summ$beta4[, "Mean"]) && (summ$beta4[, "Mean"] < 4*2))
  expect_true((0.5 < summ$sigma_[, "Mean"]) && (summ$sigma_[, "Mean"] < 2))
})

test_that("custom blocking works", {
  sampler <- create_sampler(
    y ~ reg(~1, name="mu") + reg(~x1-1, name="beta1") + reg(~x2-1, name="beta2") +
      reg(~x3-1, name="beta3") + reg(~x4-1, name="beta4"),
    data=df, block=list(c("beta4", "beta2"), c("mu", "beta3"))
  )
  sim <- MCMCsim(sampler, n.iter=500, burnin=100, n.chain=2, verbose=FALSE)
  summ <- summary(sim)
  expect_true((0.5 < summ$mu[, "Mean"]) && (summ$mu[, "Mean"] < 2))
  expect_true((0.5 < summ$beta1[, "Mean"]) && (summ$beta1[, "Mean"] < 2))
  expect_true((2*0.5 < summ$beta2[, "Mean"]) && (summ$beta2[, "Mean"] < 2*2))
  expect_true((3*0.5 < summ$beta3[, "Mean"]) && (summ$beta3[, "Mean"] < 3*2))
  expect_true((4*0.5 < summ$beta4[, "Mean"]) && (summ$beta4[, "Mean"] < 4*2))
  expect_true((0.5 < summ$sigma_[, "Mean"]) && (summ$sigma_[, "Mean"] < 2))
})

test_that("using an offset for linear regression works", {
  sampler <- create_sampler(
    formula = y ~ reg(~ x1 + x2 + x3, name="beta") + offset(4*x4), data=df
  )
  expect_equal(sampler$offset, 4*df$x4)
  sim <- MCMCsim(sampler, n.iter=500, burnin=100, n.chain=2, verbose=FALSE)
  summ <- summary(sim)
  expect_true(all((0.5 * c(1,1,2,3) < summ$beta[, "Mean"]) & (summ$beta[, "Mean"] < 2 * c(1,1,2,3))))
  expect_true((0.5 < summ$sigma_[, "Mean"]) && (summ$sigma_[, "Mean"] < 2))
})

test_that("simplified formula + offset works for linear regression", {
  sampler <- create_sampler(
    formula = y ~ x1 + x2 + x3 + offset(4*x4), data=df
  )
  expect_equal(sampler$offset, 4*df$x4)
  sim <- MCMCsim(sampler, n.iter=500, burnin=100, n.chain=2, verbose=FALSE)
  summ <- summary(sim)
  expect_true(all((0.5 * c(1,1,2,3) < summ$reg1[, "Mean"]) & (summ$reg1[, "Mean"] < 2 * c(1,1,2,3))))
  expect_true((0.5 < summ$sigma_[, "Mean"]) && (summ$sigma_[, "Mean"] < 2))
})

test_that("model specification in two steps works", {
  mod <- y ~ x1 + x2 + x3 + x4
  ml_mod <- y ~ reg(mod, Q0=1e-6, name="beta")
  sampler <- create_sampler(ml_mod, data=df, block=TRUE)
  sim <- MCMCsim(sampler, n.iter=500, burnin=100, n.chain=2, verbose=FALSE)
  summ <- summary(sim)
  expect_true(all((0.5 * c(1,1,2,3,4) < summ$beta[, "Mean"]) & (summ$beta[, "Mean"] < 2 * c(1,1,2,3,4))))
  expect_true((0.5 < summ$sigma_[, "Mean"]) && (summ$sigma_[, "Mean"] < 2))
})

n <- 1000
sd0 <- 0.41
y <- 1 + rnorm(n, sd=sd0)
test_that("different priors for residual variance work", {
  sampler <- create_sampler(y ~ 1, sigma.mod=pr_fixed(0.41^2))
  sim <- MCMCsim(sampler, n.iter=500, burnin=100, n.chain=2, verbose=FALSE)
  pm_sigma <- summary(sim$sigma_)[, "Mean"]
  expect_true((0.75*sd0 < pm_sigma) && (pm_sigma < 1.25*sd0))
  sampler <- create_sampler(y ~ 1, sigma.mod=pr_invchisq(df=1, scale="modeled"))
  sim <- MCMCsim(sampler, n.iter=500, burnin=100, n.chain=2, verbose=FALSE)
  pm_sigma <- summary(sim$sigma_)[, "Mean"]
  expect_true((0.75*sd0 < pm_sigma) && (pm_sigma < 1.25*sd0))
  sampler <- create_sampler(y ~ 1, sigma.mod=pr_exp(scale=1))
  sim <- MCMCsim(sampler, n.iter=500, burnin=100, n.chain=2, verbose=FALSE)
  pm_sigma <- summary(sim$sigma_)[, "Mean"]
  expect_true((0.75*sd0 < pm_sigma) && (pm_sigma < 1.25*sd0))
  sampler <- create_sampler(y ~ 1, sigma.mod=pr_gig(a=1, b=1, p=1))
  sim <- MCMCsim(sampler, n.iter=500, burnin=100, n.chain=2, verbose=FALSE)
  pm_sigma <- summary(sim$sigma_)[, "Mean"]
  expect_true((0.75*sd0 < pm_sigma) && (pm_sigma < 1.25*sd0))
})
