test_that('AR 1', {
  data("x_ao_mat")
  p <- 1
  q <- 0
  beta_initial <- c(0.7275)

  # 3.1797
  res <- bip_tau_resid_sc(x_ao_mat, beta_initial, p, q)$scale
})

test_that('tau-scale', {
  data("x_ao_mat")

  # 4.2918
  res <- tau_scale(x_ao_mat)
  expect_equal(round(res, digits = 4), 4.2918)

})
