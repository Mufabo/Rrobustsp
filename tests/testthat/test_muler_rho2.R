
test_that('muler_rho2', {
  library(Rrobustsp)

  data("PPG_Data") # from Rrobustsp

  ibi_ppg <- diff(PPG_Data$ppg.pos[1,])

  x <- ibi_ppg[length(ibi_ppg):1] - median(ibi_ppg)

  x2 <- x / m_scale(x)

  res <- muler_rho2(x2)

  load(path_test('res_muler_rho2'))

  expect_equal(round(res, digits = 4), unname(res_muler_rho2))
})

test_that('muler_rho2', {
  library(Rrobustsp)

  data("PPG_Data") # from Rrobustsp

  ibi_ppg <- diff(PPG_Data$ppg.pos[1,])

  x <- ibi_ppg[length(ibi_ppg):1] - median(ibi_ppg)

  x <- x / m_scale(x)

  #res <- muler_rho2(x2) ----

  rho <- rep(3.25, length(x))

  intv <- (abs(x) > 2) & (abs(x) <= 3)
  rho[intv] <-0.002 * x[intv]^8- 0.052 * x[intv]^6+ 0.432 * x[intv]^4- 0.972 * x[intv]^2+ 1.792

  rho[abs(x) <= 2] <- 0.5 * x[abs(x) <= 2]^2
  # ---------

  load(path_test('res_muler_rho2'))

  expect_equal(round(rho, digits = 4), unname(res_muler_rho2))
})

