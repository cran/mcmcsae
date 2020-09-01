
context("Model with spline component")

set.seed(1, kind="Mersenne-Twister", normal.kind="Inversion")

x <- seq(0, 1, 0.01)
dat <- data.frame(x, y=x*sin(5*x))

test_that("spline model works", {
  sampler <- create_sampler(y ~ gen(factor = ~ spline(x, knots=15)),
                            linpred="fitted", data=dat)
  sim <- MCMCsim(sampler, n.iter=500, verbose=FALSE)
  (summ <- summary(sim))
  expect_true(max(abs(dat$y - summ$linpred_[, "Mean"])) < 0.02)
})
