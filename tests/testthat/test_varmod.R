
context("Variance modeling")

set.seed(1, kind="Mersenne-Twister", normal.kind="Inversion")

n <- 100L
df <- data.frame(x=runif(n), f=factor(sample(1:8, n, replace=TRUE)))
Vmodel <- ~ vfac(factor="f", prior=pr_invchisq(df=5))
dat <- generate_data(~ reg(~ x + f, prior=pr_normal(precision=1)),
                     sigma.fixed=TRUE, formula.V=Vmodel, data=df)

test_that("generated data based on vfac model is OK", {
  expect_true(length(dat$y) == n)
  expect_false(any(is.na(dat$y)))
})

test_that("modeling vfac variance structure works", {
  df$y <- dat$y
  Vmodel <- ~ vfac(factor="f", prior=pr_invchisq(df=2), name="varf")
  sampler <- create_sampler(y ~ x + f, sigma.fixed=TRUE, formula.V=Vmodel, data=df)
  sim <- MCMCsim(sampler, burnin=100, n.iter=400, n.chain=2, store.all=TRUE, verbose=FALSE)
  expect_is(sim$varf, "dc")
  summ <- summary(sim)
  compute_DIC(sim)
  compute_WAIC(sim)
})

# non-diagonal sampling covariance matrix (NB very contrived example)
subdiv <- c(25,10,25,25,8,7)
df$f <- as.factor(rep(1:6, subdiv))
Q0 <- bdiag(mapply(Q_RW2, subdiv)) + Diagonal(n)
dat <- generate_data(~ reg(~ x + f, prior=pr_normal(precision=1)),
         sigma.fixed=TRUE, Q0=Q0, formula.V=Vmodel, data=df)

test_that("vfac variance structure works for compatible non-diagonal sampling variance matrix", {
  df$y <- dat$y
  Vmodel <- ~ vfac(factor="f", prior=pr_invchisq(df=2), name="varf")
  sampler <- create_sampler(y ~ x + f, sigma.fixed=TRUE, Q0=Q0, formula.V=Vmodel, data=df)
  sim <- MCMCsim(sampler, burnin=100, n.iter=400, n.chain=2, store.all=TRUE, verbose=FALSE)
  expect_is(sim$varf, "dc")
  summ <- summary(sim)
  compute_DIC(sim)
  compute_WAIC(sim)
})


Vmodel <- ~ vreg(~ x, prior=pr_normal(precision=1))
dat <- generate_data(~ reg(~ x + f, prior=pr_normal(precision=1)),
         sigma.fixed=TRUE, formula.V=Vmodel, data=df)

test_that("generated data based on vreg model is OK", {
  expect_true(length(dat$y) == n)
  expect_false(any(is.na(dat$y)))
})

test_that("modeling vreg variance structure works", {
  df$y <- dat$y
  Vmodel <- ~ vreg(formula = ~ x, name="varf")
  sampler <- create_sampler(y ~ x + f, sigma.fixed=TRUE, formula.V=Vmodel, data=df)
  sim <- MCMCsim(sampler, burnin=100, n.iter=400, n.chain=2, store.all=TRUE, verbose=FALSE)
  expect_is(sim$varf, "dc")
  summ <- summary(sim)
  compute_DIC(sim)
  compute_WAIC(sim)
})
