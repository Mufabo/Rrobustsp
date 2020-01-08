test_that('p = 0', {
  set.seed(1)
  x <- rnorm(5)
  p <- 0

  tmp <- ar_est_bip_tau(x, p)
})

test_that('p = 1', {
  # phi_hat =
  #
  #   -0.2290
  #
  # a_scale_final =
  #
  #   0.9736    0.9837
  #
  # x_filt =
  #        0    0.1836   -0.8356    1.5953    0.3295


})

test_that('p > 1', {
#   phi_hat =
#
#     0.0441    0.0481   -0.9900
#
#
#   x_filt =
#
# 0 0  0  1.5953  0.3295
#
#   a_scale_final =
#
#     0.9736    0.9837    1.1303    0.8811

  set.seed(1)
  x <- rnorm(5)
  p <- 3

  tmp <- ar_est_bip_tau(x, p)
})
