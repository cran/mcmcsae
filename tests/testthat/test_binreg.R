
context("Binomial regression")

# for reproducibility, even over platforms:
set.seed(1, kind="Mersenne-Twister", normal.kind="Inversion")

n <- 1000L
df <- data.frame(
  x1 = rnorm(n),
  x2 = runif(n)
)
b <- c(0.8, 2, 1)
df$y <- rbinom(n, 1, prob=1/(1 + exp(-(b[1] + b[2]*df$x1 + b[3]*df$x2))))

test_that("logistic regression works", {
  sampler <- create_sampler(y ~ reg(~ x1+x2, Q0=0.1), family="binomial", data=df)
  sim <- MCMCsim(sampler, n.iter=600, burnin=250, n.chain=2, verbose=FALSE)
  summ <- summary(sim)
  expect_true(all.equal(unname(summ$reg1[, "Mean"]), b, tolerance=1))
  compute_DIC(sim)
  compute_WAIC(sim, show.progress=FALSE)
})

df$y <- rbinom(n, 1, prob=pnorm(b[1] + b[2]*df$x1 + b[3]*df$x2))
test_that("probit regression works", {
  sampler <- create_sampler(y ~ reg(~  x1 + x2), ny=1, family=f_binomial(link="probit"), data=df)
  sim <- MCMCsim(sampler, n.iter=600, burnin=250, n.chain=2, verbose=FALSE)
  summ <- summary(sim)
  expect_true(all.equal(unname(summ$reg1[, "Mean"]), b, tolerance=1))
  compute_DIC(sim)
  compute_WAIC(sim, show.progress=FALSE)
})

y <- 30; ny <- 100
test_that("minimal logistic binomial regression example works", {
  sampler <- create_sampler(y ~ 1, ny=ny, family="binomial")
  sim <- MCMCsim(sampler, n.iter=500, burnin=200, n.chain=2, verbose=FALSE)
  summ <- summary(transform_dc(sim$reg1, fun=function(x) 1/(1+exp(-x))))
  expect_true(0.2 < summ[, "Mean"] && summ[, "Mean"] < 0.4)
  summ <- summary(predict(sim, show.progress=FALSE))  # in-sample
  expect_true(20 < summ[, "Mean"] && summ[, "Mean"] < 40)
  summ <- summary(predict(sim, type="response", show.progress=FALSE))  # in-sample probability
  expect_true(0.2 < summ[, "Mean"] && summ[, "Mean"] < 0.4)
  summ <- summary(predict(sim, newdata=data.frame(1), ny=200, show.progress=FALSE))  # out-of-sample
  expect_true(40 < summ[, "Mean"] && summ[, "Mean"] < 80)
  summ <- summary(predict(sim, newdata=data.frame(1), ny=200, type="response", show.progress=FALSE))  # out-of-sample probability
  expect_true(0.2 < summ[, "Mean"] && summ[, "Mean"] < 0.4)
})

y <- rbinom(1000L, 1L, prob=0.3)
test_that("minimal probit binary regression example works", {
  sampler <- create_sampler(y ~ 1, family=f_binomial(link="probit"))
  sim <- MCMCsim(sampler, n.iter=500, burnin=200, n.chain=2, verbose=FALSE)
  summ <- summary(transform_dc(sim$reg1, fun=function(x) pnorm(x)))
  expect_true(0.2 < summ[, "Mean"] && summ[, "Mean"] < 0.4)
  summ <- summary(predict(sim, type="response", iters=sample(1:500, 250), show.progress=FALSE))  # in-sample probability
  expect_true(all(0.2 < summ[, "Mean"] & summ[, "Mean"] < 0.4))
  summ <- summary(predict(sim, newdata=data.frame(1), show.progress=FALSE))  # out-of-sample
  expect_true(0.2 < summ[, "Mean"] && summ[, "Mean"] < 0.4)
  summ <- summary(predict(sim, newdata=data.frame(1), type="response", show.progress=FALSE))  # out-of-sample probability
  expect_true(0.2 < summ[, "Mean"] && summ[, "Mean"] < 0.4)
})
