library(Rrobustsp)

test_that('spatmed', {
  X <- matrix(c(0,1,0,0,1,0,0,0),4,2)
  expect_equal(spatmed(X), c(0.5, 0.5))
  expect_equal(round(spatmed(X+1i), digits = 4), c(0.0000+1i, 0.0000+1i))
})
