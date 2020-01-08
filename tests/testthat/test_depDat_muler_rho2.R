test_that('muler_rho2', {
  # Rrobustsp/R/dependentData_Auxiliary.R
  data('x_ao_mat')
  sigma_m <- 1.943347

  # 347.3194
  res <- sum(muler_rho2(x_ao_mat/sigma_m))
})
