
context("Truncated multivariate normal sampling")

set.seed(1, kind="Mersenne-Twister", normal.kind="Inversion")

# 5d example, 1 equality, one inequality
R <- cbind(c(1,-1,-1,-1,1))
S <- cbind(c(0,0,0,1,-1))
mu <- c(6.9, 4.2, 1.0, 4.9, 4.3)
sd0 <- c(0.3, 0.2, 0.1, 0.3, 0.3)

# 'exact' answer
answer <- c(6.591, 4.337, 1.034, 5.209, 3.991)
answer.sd <- c(0.254, 0.187, 0.098, 0.254, 0.254)


test_that("HMC TMVN method works", {
  sampler <- create_TMVN_sampler(Q=diag(1/sd0^2), mu=mu, R=R, r=0, S=S, s=0)
  sim <- MCMCsim(sampler, n.iter=2500, verbose=FALSE)
  summ <- summary(sim)
  expect_equal(crossprod_mv(R, summ$x[, "Mean"]), 0)
  expect_true(crossprod_mv(S, summ$x[, "Mean"]) > 0)
  expect_true(all.equal(unname(summ$x[, "Mean"]), answer, tolerance=0.04))
  expect_true(all.equal(unname(summ$x[, "SD"]), answer.sd, tolerance=0.02))
})

test_that("HMC TMVN method works after projection on equality constraint surface", {
  sampler <- create_TMVN_sampler(Q=diag(1/sd0^2), mu=mu, R=R, r=0, S=S, s=0, reduce=TRUE)
  sim <- MCMCsim(sampler, n.iter=2500, verbose=FALSE)
  summ <- summary(sim)
  expect_equal(crossprod_mv(R, summ$x[, "Mean"]), 0)
  expect_true(crossprod_mv(S, summ$x[, "Mean"]) > 0)
  expect_true(all.equal(unname(summ$x[, "Mean"]), answer, tolerance=0.04))
  expect_true(all.equal(unname(summ$x[, "SD"]), answer.sd, tolerance=0.02))
})

test_that("Gibbs TMVN method works", {
  sampler <- create_TMVN_sampler(Q=diag(1/sd0^2), mu=mu, R=R, r=0, S=S, s=0, method="Gibbs")
  sim <- MCMCsim(sampler, burnin=100, n.iter=2500, verbose=FALSE)
  summ <- summary(sim)
  expect_equal(crossprod_mv(R, summ$x[, "Mean"]), 0)
  expect_true(crossprod_mv(S, summ$x[, "Mean"]) > 0)
  expect_true(all.equal(unname(summ$x[, "Mean"]), answer, tolerance=0.04))
  expect_true(all.equal(unname(summ$x[, "SD"]), answer.sd, tolerance=0.02))
})

test_that("Soft TMVN method works", {
  sampler <- create_TMVN_sampler(Q=diag(1/sd0^2), mu=mu, R=R, r=0, S=S, s=0, method="softTMVN")
  sim <- MCMCsim(sampler, burnin=100, n.iter=2500, verbose=FALSE)
  summ <- summary(sim)
  expect_equal(crossprod_mv(R, summ$x[, "Mean"]), 0)
  expect_true(crossprod_mv(S, summ$x[, "Mean"]) > 0)
  expect_true(all.equal(unname(summ$x[, "Mean"]), answer, tolerance=0.04))
  expect_true(all.equal(unname(summ$x[, "SD"]), answer.sd, tolerance=0.02))
})
