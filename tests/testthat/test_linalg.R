
context("Matrix algebra")

set.seed(1, kind="Mersenne-Twister", normal.kind="Inversion")

test_that("checking for redundant columns works", {
  M <- matrix(rnorm(5*2), 5, 2)
  expect_null(detect_redundancy(crossprod(M)))
  expect_null(detect_redundancy(M, method="qr"))
  expect_identical(remove_redundancy(M), M)
  Mext <- M[, c(1, 1, 1, 2)]
  expect_identical(remove_redundancy(Mext), M)
})

test_that("is_zero_matrix works", {
  expect_false(is_zero_matrix(matrix(0.1, 2, 1)))
  expect_true(is_zero_matrix(matrix(0, 2, 1)))
  expect_false(is_zero_matrix(Diagonal(3)))
  expect_true(is_zero_matrix(0*Diagonal(3)))
  expect_false(is_zero_matrix(as(Diagonal(1), "tabMatrix")))
  expect_true(is_zero_matrix(0*as(Diagonal(1), "tabMatrix")))
})

test_that("inverseSPD works", {
  n <- 4
  M <- crossprod(matrix(rnorm(n*n), n, n)) + diag(n)
  expect_equal(solve(M), inverseSPD(M))
})

test_that("dotprodC works", {
  x <- rnorm(10)
  y <- runif(10)
  expect_equal(sum(x*y), dotprodC(x, y))
})

test_that("add_diagC works", {
  n <- 7
  M <- matrix(rnorm(n*n), n, n)
  d <- rnorm(n)
  Md <- add_diagC(M, d)
  expect_equal(Md, M + diag(d))
  expect_equal(diag(M), diag(Md) - d)
})

test_that("matrix-vector products work", {
  n <- 4
  M <- matrix(rnorm(n*n), n)
  x <- rnorm(n)
  expect_equal(as.numeric(M %*% x), M %m*v% x)
  M <- as(M, "CsparseMatrix")
  expect_equal(as.vector(M %*% x), M %m*v% x)
  M <- crossprod(M)
  expect_equal(as.vector(M %*% x), M %m*v% x)
  M <- Diagonal(n)
  expect_equal(as.vector(M %*% x), M %m*v% x)
  M <- Diagonal(x=rnorm(n))
  expect_equal(as.vector(M %*% x), M %m*v% x)
  M <- as(matrix(c(1,0,1,0,1,0,0,0,0), 3, 3), "tabMatrix")
  x <- rnorm(3)
  expect_equal(as.vector(M %*% x), M %m*v% x)
  M <- as(rnorm(3)*matrix(c(1,0,1,0,1,0,0,0,0), 3, 3), "tabMatrix")
  x <- rnorm(3)
  expect_equal(Ctab2dgC(M) %m*v% x, M %m*v% x)
  M <- as(matrix(c(1,0,0), 3, 1), "tabMatrix")
  x <- 2
  expect_equal(Ctab2dgC(M) %m*v% x, M %m*v% x)
  M <- as(matrix(c(1,0,0,0,2,0), 3, 2), "tabMatrix")
  x <- rnorm(2)
  expect_equal(Ctab2dgC(M) %m*v% x, M %m*v% x)
})

test_that("matrix-vector crossproducts work", {
  n <- 4
  M <- matrix(rnorm(n*n), n)
  x <- rnorm(n)
  expect_equal(as.numeric(crossprod(M, x)), crossprod_mv(M, x))
  M <- as(M, "CsparseMatrix")
  expect_equal(as.vector(crossprod(M, x)), crossprod_mv(M, x))
  M <- crossprod(M)
  expect_equal(as.vector(crossprod(M, x)), crossprod_mv(M, x))
  M <- Diagonal(n)
  expect_equal(x, crossprod_mv(M, x))
  M <- Diagonal(x=rnorm(n))
  expect_equal(M@x * x, crossprod_mv(M, x))
  M <- as(matrix(c(1,0,1,0,1,0,0,0,0), 3, 3), "tabMatrix")
  x <- rnorm(3)
  expect_equal(as.vector(t(M) %*% x), crossprod_mv(M, x))
  M <- as(rnorm(3)*matrix(c(1,0,1,0,1,0,0,0,0), 3, 3), "tabMatrix")
  expect_equal(crossprod_mv(Ctab2dgC(M), x), crossprod_mv(M, x))
  x <- rnorm(2)
  expect_error(crossprod_mv(M, x))  # incompatible dimensions
  M <- as(matrix(c(1,0,0), 3, 1), "tabMatrix")
  x <- rnorm(3)
  expect_equal(crossprod_mv(Ctab2dgC(M), x), crossprod_mv(M, x))
  M <- as(matrix(c(1,0,0,0,2,0), 3, 2), "tabMatrix")
  expect_equal(crossprod_mv(Ctab2dgC(M), x), crossprod_mv(M, x))
})

test_that("diagonal-dgC product works", {
  nr <- 10; nc <- 25
  Q <- Diagonal(x=rnorm(nr))
  X <- rsparsematrix(nr, nc, 0.01)
  expect_equal(Q %*% X, Cdiag_sparse_prod(Q@x, X))
})

test_that("tab to dgC conversion works", {
  m1 <- matrix(c(1,0,1,0,1,0,0,0,0), 3, 3)
  M1 <- as(m1, "tabMatrix")
  expect_false(M1@num)
  expect_equal(Ctab2dgC(M1), as(as(m1, "CsparseMatrix"), "generalMatrix"))
  m1 <- matrix(c(2,0,0), 3, 1)
  M1 <- as(m1, "tabMatrix")
  expect_true(M1@num)
  expect_equal(Ctab2dgC(M1), as(as(m1, "CsparseMatrix"), "generalMatrix"))
  M2 <- aggrMatrix(sample(1:7, 11, replace=TRUE))
  expect_equal(Ctab2dgC(M2), as(as(M2, "CsparseMatrix"), "generalMatrix"))
  M2 <- aggrMatrix(sample(1:7, 11, replace=TRUE), w=runif(11))
  expect_equivalent(as(Ctab2dgC(M2), "matrix"), Ctab2mat(M2))
  expect_error(as(matrix(c(1,1,1,0), 2, 2), "tabMatrix"))
})

test_that("tabMatrix column selection works", {
  m <- matrix(c(1,0,1,0,1,0,0,0,0), 3, 3)
  M <- as(m, "tabMatrix")
  expect_identical(M[, 1], m[, 1])
  expect_identical(Ctab2mat(M[, 1, drop=FALSE]), m[, 1, drop=FALSE])
  expect_identical(Ctab2mat(M[, c(1, 3)]), m[, c(1, 3)])
  expect_identical(Ctab2mat(M[, -2]), m[, -2])
  expect_identical(M[, c(-2, -3)], c(1, 0, 1))
  expect_identical(Ctab2mat(M[, c(-2, -3), drop=FALSE]), m[, c(-2, -3), drop=FALSE])
})

test_that("crossprod_sym works", {
  n <- 25L
  q <- runif(n)
  Q <- Diagonal(x=q)
  Qsym <- crossprod(rsparsematrix(n, n, density=0.1))
  Qmat <- as(Qsym, "matrix")
  M <- matrix(rnorm(n^2), n, n)
  expect_equal(crossprod(M, Q %*% M), crossprod_sym(M, q))
  expect_equal(crossprod_sym(M, q), crossprod_sym(M, Q))
  expect_equal(crossprod(M, Qsym %*% M), crossprod_sym(M, Qsym))
  expect_equal(crossprod(M, Qmat %*% M), crossprod_sym(M, Qmat))
  M <- Diagonal(n)
  expect_equal(crossprod(M, Q %*% M), crossprod_sym(M, q))
  expect_equal(crossprod_sym(M, q), crossprod_sym(M, Q))
  expect_equal(crossprod(M, Qsym %*% M), crossprod_sym(M, Qsym))
  expect_equivalent(as(crossprod(M, Qmat %*% M), "matrix"), crossprod_sym(M, Qmat))
  M <- Diagonal(x=rnorm(n))
  expect_equal(crossprod(M, Q %*% M), crossprod_sym(M, q))
  expect_equal(crossprod_sym(M, q), crossprod_sym(M, Q))
  expect_equal(as(crossprod(M, Qsym %*% M), "symmetricMatrix"), crossprod_sym(M, Qsym))
  expect_equivalent(as(crossprod(M, Qmat %*% M), "matrix"), crossprod_sym(M, Qmat))
  M <- aggrMatrix(sample(1:n, n))
  expect_true(tab_isPermutation(M))
  expect_equal(as(crossprod(M, Q %*% M), "diagonalMatrix"), crossprod_sym(M, q))
  expect_equal(crossprod_sym(M, q), crossprod_sym(M, Q))
  expect_equal(as(crossprod(M, Qsym %*% M), "symmetricMatrix"), crossprod_sym(M, Qsym))
  expect_equivalent(as(crossprod(M, Qmat %*% M), "matrix"), crossprod_sym(M, Qmat))
  M <- aggrMatrix(sample(1:n, n), w=runif(n))
  expect_false(tab_isPermutation(M))
  expect_equal(as(crossprod(M, Q %*% M), "diagonalMatrix"), crossprod_sym(M, q))
  expect_equal(crossprod_sym(M, q), crossprod_sym(M, Q))
  attr(M, "isPermutation") <- FALSE  # TODO add isPermutation slot to tabMatrix class definition
  expect_equal(as(crossprod(M, Qsym %*% M), "symmetricMatrix"), crossprod_sym(M, Qsym))
  expect_equivalent(as(crossprod(M, Qmat %*% M), "matrix"), crossprod_sym(M, Qmat))
  m <- 12L
  M <- rsparsematrix(n, m, density=0.1)
  expect_equal(as(crossprod(M, Q %*% M), "symmetricMatrix"), crossprod_sym(M, q))
  expect_equal(crossprod_sym(M, q), crossprod_sym(M, Q))
  expect_equal(as(crossprod(M, Qsym %*% M), "symmetricMatrix"), crossprod_sym(M, Qsym))
  expect_equivalent(as(crossprod(M, Qmat %*% M), "matrix"), crossprod_sym(M, Qmat))
})

test_that("crossprod_sym2 works", {
  n <- 11L; m <- 5L
  M1 <- matrix(runif(n*m), n, m)
  expect_equal(crossprod_sym2(M1), crossprod(M1))
  M2 <- diag(runif(n)) %*% M1
  expect_equal(crossprod_sym2(M1, M2), crossprod(M1, M2))
  M1 <- as(as(M1, "CsparseMatrix"), "generalMatrix")
  M2 <- as(as(M2, "CsparseMatrix"), "generalMatrix")
  expect_equal(crossprod_sym2(M1), crossprod(M1))
  expect_equal(crossprod_sym2(M1, M2), as(crossprod(M1, M2), "symmetricMatrix"))
})

test_that("(re)defined S4 methods work", {
  n <- 10L; m <- 5L
  M <- matrix(rnorm(n*m), n, m)
  expect_equal(Diagonal(n) %*% M, M)
  expect_equal(Diagonal(x=1:n) %*% M, (1:n) * M)
  Ms <- rsparsematrix(n, n, density=0.2)
  expect_equal(Ms %*% M, as.matrix(Ms) %*% M)
  Ms <- crossprod(Ms)
  expect_equal(Ms %*% M, as.matrix(Ms) %*% M)
  Ms <- aggrMatrix(sample(1:n, n))
  expect_equal(Ms %*% M, as.matrix(Ms) %*% M)
})
