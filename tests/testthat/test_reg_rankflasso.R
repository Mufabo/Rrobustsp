library(Rrobustsp)
library(MASS)
library(compiler)
library(Rcpp)

test_that('enetpath test 1', {
  skip('skip')
  # Arguments ----
  load('~/Rrobustsp/data/images.RData')
  load(path_test('blas'))

  y20n <- unlist(unname(images['y20n']))

  X <- diag(1, 400, 400)
  L <- 20

  lambda2 <- 340
  lambda1 <- 124

  B1 <- rankflasso(y20n, X, lambda1, lambda2, Blas)


  #MSE_rank1 <- colSums((scaledata(B1) - y20n)^2)
  #expect_equal(MSE_rank1, 1.058938435615084e+02)


})

test_that('enetpath test 1 long', {
  #skip('takes too long')
  # Args ----

  X <- diag(1, 400, 400) # Matrix::sparseMatrix(1:400, 1:400, x = 1)#
  load('~/Rrobustsp/data/images.RData')
  load(path_test('blas'))

  y <- unlist(unname(images['y20n']))
  lambda2 <- 340
  lambda1 <- 124
  b0 <- Blas

  printitn <- 0

  # ----
  n <- nrow(X)
  p <- ncol(X)

  intcpt <- F

  if(is.null(b0)) {
    b0 <- MASS::ginv(cbind(rep(1, n), X)) %*% y #qr.solve(cbind(rep(1, n), X), y)
    b0 <- bo[2:length(b0)]
  }

  B <- repmat(1:n, n)

  A <- t(B)

  a <- A[A < B]
  b <- B[A < B]

  D <- diag(-1, p-1, p-1)
  D[seq(p, (p-1)^2, p)] <- 1

  D <- cbind(D, c(rep(0, p-2), 1))

  ytilde <- c(y[a] - y[b], rep(0, p-1))

  Xtilde <- rbind( X[a,] - X[b,], lambda2 * D)

  if(printitn > 0) sprintf('rankflasso: starting iterations\n')

  r <- ladlasso(ytilde, Xtilde, lambda1, intcpt, b0, reltol = 1)
  iter <- r[[2]]

  r <- r[[1]]

  r[abs(r) < 1e-7] <- 0
  expect_equal(iter, 11)

  load(path_test('reg_rankflasso_1'))
  expect_equal(round(r, digits = 4), sol)
})

test_that('short', {
  skip('short')
  # Args ----
  X <- diag(1, 6,3)

  y <- c(0.01739562 ,-1.28630053, -1.64060553 , 0.45018710, -0.01855983 ,-0.31806837)

  lambda2 <- 3
  lambda1 <- 1
  b0 <- NULL

  printitn <- 0
  # ----
  n <- nrow(X)
  p <- ncol(X)

  if(is.null(b0)) {
    b0 <- qr.solve(cbind(rep(1, n), X), y)
    b0 <- b0[2:length(b0)]
  }

  B <- repmat(1:n, n)

  A <- t(B)

  a <- A[A < B]
  b <- B[A < B]

  D <- diag(-1, p-1, p-1)
  D[seq(p, (p-1)^2, p)] <- 1

  D <- cbind(D, c(rep(0, p-2), 1))

  ytilde <- c(y[a] - y[b], rep(0, p-1)) # stimmt

  Xtilde <- rbind( X[a,] - X[b,], lambda2 * D) # stimmt

  if(printitn > 0) sprintf('rankflasso: starting iterations\n')

  r <- ladlasso(ytilde, Xtilde, lambda1, F, b0, printitn = printitn)
  iter <- r[[2]]
  r <- r[[1]]
  r[abs(r) < 1e-7] <- 0

  expect_equal(round(r, digits = 4), c(0, -1.2255, -1.3225))
  expect_equal(iter, 2000)
})