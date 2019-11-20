library(Rrobustsp)

test_that('reg stuff', {
  y <- 1:4
  X <- matrix(c(0,1,0,0,1,0,0,0),4,2)
  expect_equal(enet(1:2, matrix(c(0, 1, 1, 0), 2, 2), c(-.1,.4),0.5), list(c(1.5, .5), 2))
})
