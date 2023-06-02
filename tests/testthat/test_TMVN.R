
context("Truncated multivariate normal sampling")

set.seed(1, kind="Mersenne-Twister", normal.kind="Inversion")

# 5d example, 1 equality, one inequality
R <- cbind(c(1,-1,-1,-1,1))  # equalities R'x = r
r <- 0
S <- cbind(c(0,0,0,1,-1))    # inequalities S'x >= s
s <- 1.7
mu0 <- c(6.9, 4.2, 1.0, 4.9, 4.3)
sd0 <- c(0.3, 0.2, 0.1, 0.3, 0.3)
Q0 <- diag(1/sd0^2)  # precision

# 'exact' answer
answer <- c(6.974, 4.167, 0.9918, 5.507, 3.693)
answer.sd <- c(0.191, 0.172, 0.097, 0.218, 0.218)


test_that("HMC TMVN method works", {
  sampler <- create_TMVN_sampler(Q=Q0, mu=mu0, R=R, r=r, S=S, s=s)
  sim <- MCMCsim(sampler, burnin=1000, n.iter=2500, verbose=FALSE)
  summ <- summary(sim)
  expect_equal(crossprod_mv(R, summ$x[, "Mean"]), r)
  expect_true(crossprod_mv(S, summ$x[, "Mean"]) >= s)
  expect_equal(unname(summ$x[, "Mean"]), answer, tolerance=0.04)
  expect_equal(unname(summ$x[, "SD"]), answer.sd, tolerance=0.02)
})

test_that("HMC TMVN method works after projection on equality constraint surface", {
  sampler <- create_TMVN_sampler(Q=Q0, mu=mu0, R=R, r=r, S=S, s=s, reduce=TRUE)
  sim <- MCMCsim(sampler, burnin=1000, n.iter=2500, verbose=FALSE)
  summ <- summary(sim)
  expect_equal(crossprod_mv(R, summ$x[, "Mean"]), r)
  expect_true(crossprod_mv(S, summ$x[, "Mean"]) >= s)
  expect_equal(unname(summ$x[, "Mean"]), answer, tolerance=0.04)
  expect_equal(unname(summ$x[, "SD"]), answer.sd, tolerance=0.02)
})

test_that("HMCZigZag TMVN method works", {
  sampler <- create_TMVN_sampler(Q=Q0, mu=mu0, R=R, r=r, S=S, s=s,
                                 method=m_HMCZigZag(rate=sqrt(diag(Q0))/5))
  sim <- MCMCsim(sampler, burnin=50, n.iter=100, verbose=FALSE)
  summ <- summary(sim)
  expect_equal(crossprod_mv(R, summ$x[, "Mean"]), r, tolerance=0.2)
  expect_true(crossprod_mv(S, summ$x[, "Mean"]) >= s)
  expect_equal(unname(summ$x[, "Mean"]), answer, tolerance=0.2)
  expect_equal(unname(summ$x[, "SD"]), answer.sd, tolerance=0.2)
})

test_that("Gibbs TMVN method works", {
  sampler <- create_TMVN_sampler(Q=Q0, mu=mu0, R=R, r=r, S=S, s=s, method="Gibbs")
  sim <- MCMCsim(sampler, burnin=500, n.iter=3000, verbose=FALSE)
  summ <- summary(sim)
  expect_equal(crossprod_mv(R, summ$x[, "Mean"]), r)
  expect_true(crossprod_mv(S, summ$x[, "Mean"]) >= s)
  expect_equal(unname(summ$x[, "Mean"]), answer, tolerance=0.04)
  expect_equal(unname(summ$x[, "SD"]), answer.sd, tolerance=0.02)
})

test_that("slice-Gibbs TMVN method works", {
  sampler <- create_TMVN_sampler(Q=Q0, mu=mu0, R=R, r=r, S=S, s=s, method=m_Gibbs(slice=TRUE))
  sim <- MCMCsim(sampler, burnin=500, n.iter=3000, verbose=FALSE)
  summ <- summary(sim)
  expect_equal(crossprod_mv(R, summ$x[, "Mean"]), r)
  expect_true(crossprod_mv(S, summ$x[, "Mean"]) >= s)
  expect_equal(unname(summ$x[, "Mean"]), answer, tolerance=0.05)
  expect_equal(unname(summ$x[, "SD"]), answer.sd, tolerance=0.05)
})

test_that("Soft TMVN method works", {
  sampler <- create_TMVN_sampler(Q=Q0, mu=mu0, R=R, r=r, S=S, s=s, method="softTMVN")
  sim <- MCMCsim(sampler, burnin=500, n.iter=2500, verbose=FALSE)
  summ <- summary(sim)
  expect_equal(crossprod_mv(R, summ$x[, "Mean"]), r)
  expect_true(crossprod_mv(S, summ$x[, "Mean"]) >= s)
  expect_equal(unname(summ$x[, "Mean"]), answer, tolerance=0.04)
  expect_equal(unname(summ$x[, "SD"]), answer.sd, tolerance=0.02)
})
