library(Matrix)

test_that('ranklassopath 1', {

  data('prostate')

  X <- prostate$X
  y <- c(prostate$y)

  tmp <- ranklassopath(y, X)

  B <- tmp[[1]]
  stats <- tmp[[3]]

  load(path_test('Brlad'))
  load(path_test('statsrlad'))

  expect_equal(B, Brlad)
  expect_equal(stats[[1]], statsrlad[[1]])
  expect_equal(stats[[2]], statsrlad[[2]])
  expect_equal(stats[[3]], statsrlad[[3]])
})
