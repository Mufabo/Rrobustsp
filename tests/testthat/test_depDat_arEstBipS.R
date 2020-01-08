test_that('p = 0', {
  set.seed(1)
  x <- rnorm(5)
  p <- 0

  tmp <- ar_est_bip_s(x, p)
})

test_that('p = 1', {
  set.seed(1)
  x <- rnorm(5)
  p <- 1

  tmp <- ar_est_bip_s(x, 1)
})

test_that('p > 1', {
  # >> phi_hat
  #
  # phi_hat =
  #
  #   0.3829    0.3849   -0.9900
  #
  # >> x_filt
  #
  # x_filt =
  #
  #   0
  # 0
  # 0
  # 1.5953
  # 0.3295
  #
  # >> a_scale_final
  #
  # a_scale_final =
  #
  #   0.8721    0.8607    1.0843    1.1368

  set.seed(1)
  x <- rnorm(5)
  p <- 3

  tmp <- ar_est_bip_s(x, p)
})
