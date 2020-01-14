test_that('check that bip_s_resid_sc returns the correct cleaned signal', {

  '
  rng(0)
  a = randn(500, 1);
  %% AR(1)
  p = 1;
  q = 0;
  x_ar1 = filter(1, [1 -.8], a);
  '
  load(path_test('x_ar1'))

  load(path_test('bipsresidsc2_ar1'))

  p <- 1
  q <- 0

  # scale 1.8860
  result <- bip_s_resid_sc(x_ar1, c(1, -.8), p, q)
  expect_equal(result$filtered_signal[,1], unname(bipsresidsc2_ar1))
})
