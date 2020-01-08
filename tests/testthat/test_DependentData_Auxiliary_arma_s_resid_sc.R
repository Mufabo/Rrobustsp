test_that('p=1, q=0', {
  data('x_ao_mat')
  p <- 1
  q <- 0

  beta_initial <- c(0.7275)
  res <- arma_s_resid_sc(x_ao_mat, beta_initial, p, q)
  expect_equal(round(res, digits = 4), 3.9014)


})
