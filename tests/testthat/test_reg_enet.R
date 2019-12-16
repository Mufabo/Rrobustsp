test_that('enet test 1', {
  # Arguments ----
  load('~/Rrobustsp/data/images.RData')
  y <- unlist(unname(images['y20n']))
  alpha <- 1
  printitn <- 0
  itermax <- 500
  load(path_test('enet_test 1_beta'))


  X <- diag(1, 400, 400)
  lambda <- 2.3441
  beta <- rep(0, 400)

  # Function ----
  p <- ncol(X)

  betaold <- beta
  normb0 <- norm(beta, type = "2")

  # test
  r <- y - X %*% beta
  expect_equal(c(r), c(y))

  lam1 <- alpha * lambda
  lam2 <- lambda * (1-alpha)

  const <- 1/(1+lam2)

  if(printitn > 0){
    sprintf('enet : using penalty lambda = %f.5', lambda)
  }

  for(iter in 1:itermax){
    for(jj in 1:p){
      beta[jj] <- const*soft_thresh(beta[jj]+ t(X[,jj]) %*% r, lam1)
      r <- r + X[,jj]*(betaold[jj]-beta[jj])
    }

    normb <- norm(beta, type = "2")

    crit <- sqrt(normb0^2 + normb^2 - 2 * Re(t(betaold) %*% beta))/normb

    if(!is.na(iter%%printitn) && (iter %% printitn) == 0){
      sprintf('enet: %4d  crit = %.8f\n',iter,crit)
    }

    if(is.nan(crit) | crit < 1e-4) {      break }

    betaold <- beta
    normb0 <- normb
  }

  expect_equal(iter, 2)
  expect_equal(beta, enet_test_1_beta)

  beta_sol <- enet(y, X, beta, lambda)
  expect_equal(beta_sol[[1]], enet_test_1_beta)
  expect_equal(beta_sol[[2]], iter)
  expect_equal(beta_sol[[2]], 2)
  # paste above return
})
