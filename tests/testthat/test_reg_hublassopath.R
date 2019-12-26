library(Matrix)
library(MASS)
test_that('hublassopath 1', {

  data('prostate')

  X <- prostate$X
  y <- c(prostate$y)

  tmp <- hublassopath(y, X)

  B <- tmp[[1]]
  stats <- tmp[[3]]

  load(path_test('Bhub'))
  load(path_test('statshub'))

  rd3 <- function(x) return(round(x, digits = 3))

  expect_equal(rd3(B), Bhub)
  expect_equal(rd3(stats[[1]]), statshub[[1]])
  expect_equal(rd3(stats[[2]]), statshub[[2]])
  expect_equal(rd3(stats[[3]]), statshub[[3]])
})
