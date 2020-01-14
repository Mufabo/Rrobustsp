
test_that('arma_est_bip_s ar1', {
  library(pracma)
  library(signal)
  library(zeallot)

  load(path_test('x_ar1'))
  load(path_test('result_s'))

  p <- 1
  q <- 0

  result <- arma_est_bip_s(x_ar1, p, q)
  expect_equal(result$cleaned_signal, result_s[[4]])
})
