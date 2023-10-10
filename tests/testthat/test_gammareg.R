
context("Gamma regression")

# for reproducibility, even across platforms:
set.seed(1, kind="Mersenne-Twister", normal.kind="Inversion")

n <- 1000L
df <- data.frame(
  x1 = rnorm(n),
  x2 = runif(n)
)
b <- c(0.8, 2, 1)
alpha <- 1
mu <- exp(b[1] + b[2]*df$x1 + b[3]*df$x2)
df$y <- rgamma(n, shape=alpha, rate=alpha/mu)

test_that("gamma regression works", {
  sampler <- create_sampler(y ~ reg(~ x1+x2), family="gamma", data=df)
  sim <- MCMCsim(sampler, n.iter=600, burnin=250, n.chain=2, verbose=FALSE)
  summ <- summary(sim)
  expect_equal(unname(summ$reg1[, "Mean"]), b, tolerance=1)
  compute_DIC(sim)
  compute_WAIC(sim, show.progress=FALSE)
})

test_that("gamma regression prediction works", {
  sampler <- create_sampler(y ~ reg(~ x1+x2), family="gamma", data=df[1:900, ])
  sim <- MCMCsim(sampler, n.iter=600, burnin=250, n.chain=2, verbose=FALSE)
  summ <- summary(sim)
  expect_equal(unname(summ$reg1[, "Mean"]), b, tolerance=1)
  pred <- predict(sim, newdata = df[901:1000, ], show.progress = FALSE)
  summpred <- summary(pred)
  #plot(summpred[, "Mean"], df$y[901:1000]); abline(0, 1)
})

alpha <- 10
df$y <- rgamma(n, shape=alpha, rate=alpha/mu)
test_that("gamma regression works when shape has to be inferred", {
  sampler <- create_sampler(y ~ reg(~ x1+x2), family=f_gamma(shape.prior = pr_gamma(1, 1)), data=df)
  sim <- MCMCsim(sampler, n.iter=600, burnin=250, n.chain=2L, verbose=FALSE)
  summ <- summary(sim)
  #plot(sim, "gamma_shape_")
  expect_length(acceptance_rates(sim)[["gamma_shape_"]], 2L)
  expect_equal(unname(summ$gamma_shape_[, "Mean"]), alpha, tolerance=1)
  expect_equal(unname(summ$reg1[, "Mean"]), b, tolerance=1)
  compute_DIC(sim)
  compute_WAIC(sim, show.progress=FALSE)
})

test_that("exponential prior on shape works", {
  sampler <- create_sampler(y ~ reg(~ x1+x2), family=f_gamma(shape.prior = pr_exp(1)), data=df)
  expect_equal(sampler$family$shape.prior$type, "gamma")
  sim <- MCMCsim(sampler, n.iter=600, burnin=250, n.chain=2, verbose=FALSE)
  summ <- summary(sim)
  expect_equal(unname(summ$gamma_shape_[, "Mean"]), alpha, tolerance=1)
  expect_equal(unname(summ$reg1[, "Mean"]), b, tolerance=1)
})
