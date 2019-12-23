library(R.matlab)
library(MASS)
library(Matrix)

test_that('rankflasso test 1', {

  # Arguments ----
  load('~/Rrobustsp/data/images.RData')
  load(path_test('blas'))

  y20n <- unlist(unname(images['y20n']))

  X <- diag(1, 400, 400)
  L <- 20

  lambda2 <- 340
  lambda1 <- 124

  B1 <- rankflasso(y20n, X, lambda1, lambda2, Blas)


  MSE_rank1 <- colSums((scaledata(B1) - y20n)^2)
  expect_equal(MSE_rank1, 1.058938435615084e+02)


})

test_that('rankflasso test 1 long', {
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

  # ladlasso ----
  # r <- ladlasso(ytilde, Xtilde, lambda1, intcpt, b0, reltol = 1)

  # ladlasso end ----
  iter <- r[[2]]

  r <- r[[1]]

  r[abs(r) < 1e-7] <- 0
  expect_equal(iter, 11)

  load(path_test('reg_rankflasso_1'))
  expect_equal(round(r, digits = 3), sol)
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

test_that('rankflasso with Ladlasso 1',{
  # Args ----

  X <- diag(1, 400, 400) # Matrix::sparseMatrix(1:400, 1:400, x = 1)#
  load('~/Rrobustsp/data/images.RData')
  load(path_test('blas'))

  y <- unlist(unname(images['y20n']))
  lambda2 <- 340
  lambda1 <- 124
  b0 <- Blas

  printitn <- 1

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

  D <- sparseMatrix(row(D)[which(!D == 0)], row(D)[which(!D == 0)])

  ytilde <- c(y[a] - y[b], rep(0, p-1))

  Xtilde <- rbind( X[a,] - X[b,], lambda2 * D)

  if(printitn > 0) sprintf('rankflasso: starting iterations\n')

  # args ladlasso ----

  # r <- ladlasso(ytilde, Xtilde, lambda1, intcpt, b0, reltol = 1)
  y <- ytilde
  X <- Xtilde
  lambda <- lambda1

  # ladlasso ----

  N <- nrow(X)
  p <- ncol(X)

  if(intcpt) X <- cbind(matrix(1, N, 1), X)

  if(is.null(b0)) b0 <- qr.solve(X, y) # ginv(X) %*% y

  iter <- NULL
  if(printitn > 0) sprintf('Computing the solution for lambda = %.3f\n',lambda)

  # The case of only one predictor
  if(p == 1){
    if(!intcpt) b1 <- wmed( rbind(y / X, 0), rbind(abs(X), lambda))
    if(!is.complex(y) & N < 200 & intcpt){
      if(lambda == 0){
        b <- elemfits(X[,2], y) # b is a matrix
        b <- b[[1]]
      }else{
        b <- elemfits(c(X[,2], 0), c(y, lambda))
        b <- b[[1]]
      }}
    res <- colSums(abs(repmat(y, ncol(b)) - X %*% b))
    indx <- which.min(res)
    b1 <- b[,indx]
  }
  else {
    # use IRWLS always when p > 1
    if(printitn > 0) print('Starting the IRWLS algorithm..\n')
    if(lambda > 0){
      y <- c(y, rep(0, p))

      # slow
      if(intcpt) X <- rbind(X, cbind(rep(0, p), diag(lambda, p, p))) else X <- rbind(X, diag(lambda, p, p))
    }

    for(iter in 1:iter_max){
      resid <- abs(y - X %*% b0)
      resid[resid < 1e-6] <- 1e-6

      # slow
      Xstar <- sweep(X, 1, resid, FUN = "/")

      # slow
      # b1 <- fRcpp(Xstar, X, matrix(y))

      # b1 <- qr.solve(Xstar %*% X, Xstar %*% y)

      # b1 <- c(MASS::ginv(t(Xstar) %*% X) %*% (t(Xstar) %*% y))

      crit <- norm(b1-b0, type = "2") / norm(b0, type = "2")

      if(printitn > 0 & iter %% printitn) sprintf('ladlasso: crit(%4d) = %.9f\n',iter,crit)
      if(crit < reltol && iter > 10) break
      b0 <- b1
    }}
  # ladlasso end ----

  b1[abs(b1) < 1e-7] <- 0
  expect_equal(iter, 11)

  load(path_test('reg_rankflasso_1'))
  expect_equal(round(b1, digits = 3), sol)
})
