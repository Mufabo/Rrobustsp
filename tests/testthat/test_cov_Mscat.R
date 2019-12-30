
test_that('mscat Tyler',{
  load(path_test('mscat_x')) # x is the input
  load(path_test('mscat'))
  result <- Mscat(Conj(t(x)), 't_loss', 0)
  flag <- T

  expect_equal(result$C, mscat[[1]])
  expect_equal(result$invC, mscat[[2]])
  expect_equal(result$flag, flag)
  expect_equal(result$iter, 9)
})
