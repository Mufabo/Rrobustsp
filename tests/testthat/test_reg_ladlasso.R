library(Rrobustsp)
library(MASS)
library(Matrix)
test_that('ladlasso1 p>1, no intcpt', {

  data('images')
  load(path_test('Blas'))

  load(path_test('sol1'))
  y20n <- unlist(unname(images['y20n']))
  lambda1 <- 124
  intcpt <- F

  Rsol <- ladlasso(y20n, diag(1, 400, 400), lambda1, intcpt, Blas)

  R_sol <- Rsol[[1]]
  R_it  <- Rsol[[2]]


  expect_equal(R_it, 11)
  expect_equal(round(R_sol - sol, digits = 3), rep(0, length(sol)))
})

test_that('ladlasso2 intcpt, p>1', {

  data('images')

  y20n <- unlist(unname(images['y20n']))
  lambda1 <- 124
  intcpt <- T

  load(path_test('sol2'))

  Rsol <- ladlasso(y20n, diag(1, 400, 400), lambda1, intcpt)

  R_it  <- Rsol[[2]]
  R_sol <- Rsol[[1]] # N x 1 matrix, whereas sol is N


  expect_equal(round(R_sol - sol, digits = 3), rep(0,length(sol))) # one of 401 mismatches with a difference of 0.003
  expect_equal(R_it, 20) # returns 41 instead of 20 ...
})

#test_that('ladlasso 3, intcpt, p=1'{
# eye-check, is correct
#})

test_that('rankflasso test 1 long', {
  skip('takes forever')
  # Args ----

  X <- diag(1, 400, 400)
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

  r <- ladlasso(ytilde, Xtilde, lambda1, F, b0, printitn = printitn)
  iter <- r[[2]]

  r <- r[[1]]

  r[abs(r) < 1e-7] <- 0
  expect_equal(iter, 11)

  load(path_test('reg_rankflasso_1'))
  expect_equal(r, reg_rankflasso_1)
})

test_that('rankflasso short', {
  skip('takes forever')
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

  r <- ladlasso(ytilde, Xtilde, lambda1, F, b0, printitn)
  iter <- r[[2]]
  r <- r[[1]]
  r[abs(r) < 1e-7] <- 0

  expect_equal(round(r, digits = 4), c(0, -1.2255, -1.3225))
  expect_equal(iter, 49)
})
