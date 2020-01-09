test_that('', {
  beta <- c(1.3667  , -1.5199  ,  1.1383  , -0.8206  , -0.0255  , -0.8495 ,  -0.0965)
  load(path_test('arma43'))

  p <- 4
  q <- 3

  res <- arma_tau_resid_sc(arma43, beta, p, q)

})

