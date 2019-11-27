library(Rrobustsp)

test_that('reg stuff', {
  y <- 1:4
  X <- matrix(c(0,1,0,0,1,0,0,0),4,2)
  expect_equal(enet(1:2, matrix(c(0, 1, 1, 0), 2, 2), c(-.1,.4),0.5), list(c(1.5, .5), 2))

})

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
  expect_equal(hublasso(y, X, lambda, c, sig0, reltol, printitn),
               list(

               ))
})
