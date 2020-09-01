
context("Negative binomial regression")

set.seed(1, kind="Mersenne-Twister", normal.kind="Inversion")

# generate negative binomial data
n <- 1000L
r <- 0.17  # overdispersion parameter
df <- data.frame(x1=rnorm(n), x2=runif(n))
mod <- ~ reg(~ 1 + x1 + x2, b0=c(3, 1, 0.5), Q0=1e10*diag(3), name="beta")
dat <- generate_data(mod, family="negbinomial", ry=r, data=df)

test_that("fitting a negative binomial model works", {
  sampler <- create_sampler(dat$y ~ reg(~ 1 + x1 + x2, name="beta"), data=df, family="negbinomial", ry=r)
  sim <- MCMCsim(sampler, n.iter=500L, burnin=150L, n.chain=2L, verbose=FALSE)
  summ <- summary(sim)
  expect_equivalent(summ$beta[, "Mean"], dat$pars$beta, tol=0.5)
  DIC <- compute_DIC(sim)
  WAIC <- compute_WAIC(sim)
  expect_true(abs((DIC["DIC"] - WAIC["WAIC2"])/abs(DIC["DIC"])) < 0.05)
  pred <- as.matrix(predict(sim, type="response", iters=1:100, show.progress=FALSE))
  expect_equivalent(r*mean(pred/(1 - pred)), mean(dat$y), tol=2)
})

test_that("fitting a negative binomial model with unknown overdispersion parameter works", {
  sampler <- create_sampler(dat$y ~ x1 + x2, data=df, family="negbinomial", r.mod=pr_invchisq(df=1, scale="modeled"))
  sim <- MCMCsim(sampler, n.iter=500L, burnin=150L, n.chain=2L, verbose=FALSE)
  summ <- summary(sim)
  expect_true(abs(summ$negbin_r_[, "Mean"] - 0.17) < 0.05)
  compute_DIC(sim)
})
