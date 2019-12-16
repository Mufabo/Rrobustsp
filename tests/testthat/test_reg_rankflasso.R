library(Rrobustsp)

test_that('enetpath test 1', {
  # Arguments ----
  load('~/Rrobustsp/data/images.RData')
  load(path_test('enetpath_test_1_stats'))
  load(path_test('enetpath_test_1_B'))

  y20n <- unlist(unname(images['y20n']))
  y20 <- unlist(unname(images['y20']))
  X <- diag(1, 400, 400)
  L <- 20

  lasso_sol <- enetpath(y20n, X, 1, L, eps = 1e-3, intcpt = F)

  Blas20n <- lasso_sol[[1]]
  statsn <- lasso_sol[[2]]

  Blas20n <- Blas20n[,2:ncol(Blas20n)]

  ero <- scaledata(Blas20n) - repmat(y20, ncol(Blas20n))

  indx <- which.min(colSums(ero^2))

  Blas <- Blas20n[,indx]
  lam_las <- statsn[['Lambda']][1+indx]



  lambda2 <- 340
  lambda1 <- 124

  B1 <- rankflasso(y20n, diag(1, n, n), lambda1, lambda2, Blas, 1)[[1]]


  MSE_rank1 <- colSums((scaledata(B1) - y20n)^2)
  expect_equal(MSE_rank1, )


})
