
context("Generate data")

set.seed(1, kind="Mersenne-Twister", normal.kind="Inversion")

test_that("prior predictive sampling works", {
  df <- data.frame(int=rep(1, 1000))
  dat <- generate_data(~ reg(~ 1, b0=2, Q0=1e10), sigma.mod=pr_fixed(value=4), data=df)
  expect_true((1.5 < mean(dat$y) && mean(dat$y) < 2.5))
  expect_true((1.5 < sd(dat$y) && sd(dat$y) < 2.5))
})

test_that("prior specification of regression component works", {
  library(survey)
  data(api)
  mod <- api00 ~
    reg(~ (api99 + stype) * sch.wide, b0=1:8, Q0=c(1e3,0.01,1,1e6,100,1,1,1e6)) +
    gen(~ api99, factor= ~ cname)
  samplr <- create_sampler(mod, sigma.mod=pr_invchisq(df=10, scale=1), data=apisrs)
  sim <- MCMCsim(samplr, verbose=FALSE)
  summ <- summary(sim)
  expect_equal(as.numeric(summ$reg1[c(4,8), "Mean"]), c(4,8), tol=1e-2)
  # prior sampling
  sim <- MCMCsim(samplr, from.prior=TRUE, verbose=FALSE)
  summ <- summary(sim)
  expect_equal(as.numeric(summ$reg1[c(4,8), "Mean"]), c(4,8), tol=1e-2)
})
