library(Rrobustsp)

test_that("All Auxi functions for x=1.0,2.1,4.9 and evtl. c=2.2", {
  x <- c(1.0,2.1,4.9)
  expect_equal(round(madn(x),digits=4), c(1.6297))
  expect_equal(round(psihub(x,2.2),digits=4), c(1.0,2.1,2.2))
  expect_equal(round(psituk(x,2.2),digits=4), c(0.6295, 0.0166,0.0))
  expect_equal(round(rhohub(x,2.2),digits=4), c(1.0,4.41,16.72))
  expect_equal(round(rhotuk(x,2.2),digits=4), c(0.8076, 1.6122, 1.6133))
  expect_equal(round(soft_thresh(x,2.2),digits=4), c(0.0,0.0,2.7))
  expect_equal(round(whub(x,2.2),digits=4), c(1.0, 1.0, 0.449))
  expect_equal(round(wtuk(x,2.2),digits=4), c(0.6295, 0.0079, 0.0))
})

test_that("rho, soft_thresh and weights for x=1.0i,2.1i,4.9i and evtl. c=2.2", {
  x <- c(1.0+1i,2.1+1i,4.9+1i)
  expect_equal(round(rhohub(x,2.2),digits=4), c(2.0000, 5.3941, 17.1644))
  expect_equal(round(rhotuk(x,2.2),digits=4), c(1.2874, 1.6133, 1.6133))
  expect_equal(round(soft_thresh(x,2.2),digits=4), c(0.0000 + 0.0000i, 0.1137 + 0.0541i, 2.7444 + 0.5601i))
  expect_equal(round(whub(x,2.2),digits=4), c(1.0000 + 0.0000i, 1.0000 + 0.0000i, 0.4310 - 0.0880i))
  expect_equal(round(wtuk(x,2.2),digits=4), c(0.8292 - 0.8264i, -0.6657 - 0.5128i, 0.0000 + 0.0000i))
})


