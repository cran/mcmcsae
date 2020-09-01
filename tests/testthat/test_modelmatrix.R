
context("Model matrix")

set.seed(1, kind="Mersenne-Twister", normal.kind="Inversion")

library(survey)
data(api)
model <- ~ stype + api99 + sch.wide*meals + ell + pct.resp

test_that("model.matrix and model_matrix yield same result", {
  M1 <- model.matrix(model, apipop)
  M2 <- model_matrix(model, apipop, sparse=FALSE)
  expect_equivalent(M1, M2)
  expect_identical(order_interactions(colnames(M1)), order_interactions(colnames(M2)))
})

test_that("contrasts as a list work", {
  expect_equal(
    model_matrix(model, apipop, contrasts.arg=list(sch.wide=levels(apipop$sch.wide)[1L], stype=levels(apipop$stype)[1L])),
    model_matrix(model, apipop)
  )
  expect_equal(
    model_matrix(model, apipop, contrasts.arg=list(stype="E")),
    model_matrix(model, apipop)
  )
})

test_that("contrasts work as intended", {
  n <- 100
  dat <- data.frame(
    f = as.factor(sample(1:3, n, replace=TRUE))
  )
  expect_equivalent(
    model.matrix(~ f, dat),
    model_matrix(~ f, dat, sparse=FALSE)
  )
  expect_equivalent(
    model.matrix(~ 0 + f, dat),
    model_matrix(~ 0 + f, dat, sparse=FALSE)
  )
  expect_equivalent(
    model_matrix(~ f, dat, contrasts.arg="contr.none", sparse=FALSE),
    cbind(1, model.matrix(~ 0 + f, dat))
  )
  expect_equal(sum(abs(
    model_matrix(~ f, dat, contrasts.arg="contr.none", sparse=FALSE) -
      model_matrix(~ f, dat, contrasts.arg="contr.none", sparse=TRUE)
  )), 0
  )
})

test_that("aggregation with model_matrix works", {
  expect_equivalent(
    model_matrix(model, apipop, by=apipop$sch.wide, sparse=FALSE),
    rowsum(model.matrix(model, apipop), apipop$sch.wide)
  )
  expect_equivalent(
    as.matrix(model_matrix(model, apipop, by=apipop$sch.wide, sparse=TRUE)),
    rowsum(model.matrix(model, apipop), apipop$sch.wide)
  )
})

test_that("matrix term works", {
  x <- rnorm(3)
  M <- matrix(1:6, 3, 2)
  a <- model.matrix(~ M)
  b <- model_matrix(~ M)
  expect_equivalent(a, b)
  expect_identical(colnames(a), colnames(b))
  a <- model.matrix(~ x + M + -1)
  b <- model_matrix(~ x + M + -1)
  expect_equivalent(a, b)
  expect_identical(colnames(a), colnames(b))
})

test_that("mix of various term types works", {
  n <- 100
  dat <- data.frame(
    x = rnorm(n),
    f = factor(sample(1:10, n, replace=TRUE)),
    g = factor(sample(1:25, n, replace=TRUE)),
    st = sample(c("a", "b", "c"), n, replace=TRUE),
    i = 1:n,
    d = as.Date(sample(c("21/5/69", "22/7/67"), n, replace=TRUE), "%d/%m/%y"),
    stringsAsFactors = FALSE
  )
  mod <- ~ d + i + g + st + g:f + x + f + I(g) + I(2*x) + g*x + sin(x) + interaction(g,f)
  a <- model.matrix(mod, dat)
  b <- model_matrix(mod, dat, sparse=FALSE)
  expect_equivalent(a, b)
  expect_identical(order_interactions(colnames(a)), order_interactions(colnames(b)))
  a <- model_matrix(mod, dat, sparse=TRUE)
  expect_equal(sum(abs(a - b)), 0)
  expect_identical(order_interactions(colnames(a)), order_interactions(colnames(b)))
})


n <- 100
data <- data.frame(
  t = 1:n,
  x = rnorm(n),
  f1 = as.factor(sample(1:10, n, replace=TRUE)),
  f2 = as.factor(sample(1:16, n, replace=TRUE)),
  f3 = as.factor(sample(1:7, n, replace=TRUE)),
  f4 = as.factor(1)
)
formulas <- list(
  ~ 0, ~ 1, ~ x,  ~ 0 + x, ~ 0 + exp(x), ~ 1 + offset(x), ~ poly(x, 3),
  ~ 0 + f1 + f2, ~ 0 + f2 + f1,
  ~ 0 + f1*f2 + f2*f3 + f1*f3, ~ 0 + x:f2 + f1*f2 + f2*f3 + f1*f3,
  ~ 0 + f2:f1, ~ 0 + f1:f2 + f3:f2,
  ~ x*t, ~ I(x*t), ~ x + t*f1*f2 + f1*f2*f3,
  ~ 0 + t:x:f1:f2
)
test_that("many special case formulas work", {
  for (mod in formulas) {
    expect_equivalent(
      model.matrix(mod, data),
      model_matrix(mod, data, sparse=FALSE)
    )
  }
  expect_error(model.matrix(~ f4, data))
  expect_error(model_matrix(~ f4, data))
})

test_that("special characters in variable names pose no problem", {
  dat <- data.frame(runif(10)); names(dat) <- "a[b]"  # also for ;,! etc etc
  expect_equivalent(
    model.matrix(~ ., dat),
    model_matrix(~ ., dat, sparse=FALSE)
  )
})
