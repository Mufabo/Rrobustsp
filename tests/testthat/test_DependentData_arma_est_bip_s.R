test_that('long' , {
  data("x_ao_mat")
  p <- 1
  q <- 0
  x <- x_ao_mat
  # ----

  # Robust starting point by BIP AR-S approximation
  beta_initial <- robust_starting_point(x, p, q)$beta_initial
  beta_initial <- 0.7275 # from matlab
  # ----

  # objective function for ARMA model and BIP-ARMA model
  F <- function(beta) arma_s_resid_sc(x, beta, p, q)
  F_bip <- function(beta) bip_s_resid_sc(x, beta, p, q)[[1]]

  # 0.7709
  beta_arma <- lsqnonlin(F, beta_initial, options = list(
    'tolx' = 5e-7))$x

  # .81
  beta_bip <- lsqnonlin(F_bip, beta_initial, options = list(
    'tolx' = 5e-7))$x
  # ----

  # innovations m-scale for ARMA model
  # 3.8799
  a_sc <- arma_s_resid_sc(x, beta_arma, p, q)

  # innovations m-scale for BIP-ARMA model
  # 1.8736
  c(a_bip_sc, x_filt) %<-% bip_s_resid_sc(x, beta_bip, p, q)
})
