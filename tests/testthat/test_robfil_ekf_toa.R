test_that('test 1', {
  library(R.matlab)
  library(tensorA)
  library(pracma)

  mat <- readMat("C:/users/computer/desktop/testdaten_Rrobustsp/ekf_toa_test1.mat")

  res <- ekf_toa(mat$res.ges, mat$theta.init, mat$BS)

  expect_equal(res$th_hat , mat$a)
  expect_equal(res$P_min , mat$b)
  expect_equal(res$P , mat$c)
  expect_equal(res$parameter , mat$d)
})
