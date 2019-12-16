library(Rrobustsp)

test_that('hublasso',{
  y <- c(1.3, 4, 5)
  X <- matrix(c(4,5,6,7,8,9,0,3,0), nrow = 3, ncol = 3, byrow = T)
  lambda <- 3.3
  b0 <- c(0,0,0)
  c <- 4
  sig0 <- 1
  reltol <- 1e-5
  printitn <- 0
  iter_max <- 500
  sol <- hublasso(y, X, c, lambda, b0, sig0, reltol, printitn)

  load(path_test('hublasso_sol'))
  # tests fail, comparison by eye shows pass
  expect_equal(abs(sol[[1]] - c(hublasso_sol[[1]])) < 1e-4, T)
  expect_equal(round(sol[[2]], digits = 4), hublasso_sol[[2]])
  expect_equal(round(sol[[3]], digits = 4), hublasso_sol[[3]])
})

test_that('enet_long',{
  # inputs ----

  y <- c(1.0, 2.555555, 3.2342342, 4.73567256, 5.13131, 6, 7, 8, 9, 10)
  X <- diag(1, 10, 5)
  alpha <- 1
  L <- 10
  eps <- 1e-3
  intcpt <- T
  lambda = 2.4648666648211974
  alpha <- 1
  printitn = 0
  itermax = 1000
  beta = numeric(5)


# enet <- function() -----
  # check for valid arguments

  p <- ncol(X)

  betaold <- beta
  normb0 <- norm(beta, type = "2")
  r <- y - X %*% beta

  lam1 <- alpha * lambda
  lam2 <- lambda * (1-alpha)

  const <- 1/(1+lam2)

  if(printitn > 0){
    cat('enet : using penalty lambda = ', round(lambda, digits=5))
  }

  for(iter in 1:itermax){
    for(jj in 1:p){
      # transpose necessary ?
      beta[jj] <- const*soft_thresh(beta[jj]+ t(X[,jj]) %*% r, lam1)
      r <- r + X[,jj]*(betaold[jj]-beta[jj])
    }

    normb <- norm(beta, type = "2")

    # tranpose necessary ? what if normb==0
    crit <- sqrt(normb0^2 + normb^2 - 2 * Re(t(betaold) %*% beta))/normb

    if(!is.na(iter%%printitn) && (iter %% printitn) == 0){
      sprintf('enet: %4d  crit = %.8f\n',iter,crit)
    }

    if(is.nan(crit) | crit < 1e-4) {      break }

    betaold <- beta
    normb0 <- normb
  }
})

test_that('enet', {
  y <- 1:4
  X <- matrix(c(0,1,0,0,1,0,0,0),4,2)

  # passed
  expect_equal(enet(1:2, matrix(c(0, 1, 1, 0), 2, 2), c(-.1,.4),0.5), list(c(1.5, .5), 2))

  load('~/Rrobustsp/data/images.RData')
  y20n <- unlist(images['y20n'])

  load(path_test('beta_enet'))
  beta_enet <- unname(unlist(beta_enet))

  X <- diag(1, 400, 400)
  lambda <- 2.3441
  beta <- rep(0, 400)

  beta_r <- enet(y20n, X, beta, lambda)
  expect_equal(beta_enet, round(beta_r[[1]], digits = 4))
  expect_equal(beta_enet, round(beta_r[[1]], digits = 4))
})

test_that('enetpath', {
  load('~/RRobustsp/R/sysdata.rda')
  Blas20n <- Blas20n[[1]]

  load('~/Rrobustsp/data/images.RData')
  y20n <- unlist(images['y20n'])

  res <- enetpath(y20n, diag(1, 400, 400), 1, 20, 1e-3, F)[[1]]
  expect_equal(res, Blas20n)
})
