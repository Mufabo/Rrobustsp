test_that('arma_resid AR', {
  x <- c(-0.001527736, 0.876177007, 0.233397402, -0.409891350, -0.183410372)
  p <- 1
  q <- 0
  beta_hat <- c(1)

  res <- arma_resid(x, beta_hat, p, q)

  target <- c(0   , 0.8777 ,  -0.6428 ,  -0.6433,    0.2265)

  expect_equal(round(res, digits = 4), target)
})

test_that('arma_resid MA', {
  x <- c(-0.001527736, 0.876177007, 0.233397402, -0.409891350, -0.183410372)
  p <- 0
  q <- 1
  beta_hat <- c(1)

  res <- arma_resid(x, beta_hat, p, q)

  target <- c(0,0.8762,1.1096,0.6997 ,0.5163)

  expect_equal(round(res, digits = 4), target)
})

test_that('arma_resid ARMA', {
  x <- c(-0.001527736, 0.876177007, 0.233397402, -0.409891350, -0.183410372)
  p <- 1
  q <- 1
  beta_hat <- c(.7, 0.5)

  res <- arma_resid(x, beta_hat, p, q)

  target <- c(0,0.8772,0.0587, -0.5439, -0.1684)

  expect_equal(round(res, digits = 4), target)
})

test_that('arma_s_resid AR', {})

test_that('arma_s_resid MA', {
  x <- c(-0.001527736, 0.876177007, 0.233397402, -0.409891350, -0.183410372)
  p <- 0
  q <- 1
  beta_hat <- c(1)


  })

test_that('arma_s_resid ARMA', {})

test_that('arma_s_resid_sc AR', {})

test_that('arma_s_resid_sc MA', {})

test_that('arma_s_resid_sc ARMA', {})

test_that('arma_tau_resid_sc AR', {})

test_that('arma_tau_resid_sc MA', {})

test_that('arma_tau_resid_sc ARMA', {})

test_that('bip_resid AR', {

})

test_that('bip_resid MA', {})

test_that('bip_resid ARMA', {
  library(pracma)
  x <- c(-0.001527736, 0.876177007, 0.233397402, -0.409891350, -0.183410372)
  p <- 1
  q <- 1
  beta_hat <- c(0.7, 0.5)

  res <- bip_resid(x, beta_hat, p, q)

  target <- c(0.8772  ,  0.0955  , -0.5255,   -0.1593)

  expect_equal(round(res, digits = 4), target)
})
