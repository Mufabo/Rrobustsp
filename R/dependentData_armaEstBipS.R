#' arma_est_bip_s
#'
#' @export
arma_est_bip_s <- function(x, p, q, tolX = 1e-2){
  if(p == 0 && q == 0){
    result$inno_scale < m_scale(x)
    warning('Please choose a nonzero value for p or q')
    return(result)}
  if(length(x) <= (p + q)){
    warning('There are too many parameters to estimate for chosen data size. Reduce model order or use a larger data set.')
    return(result)
  }

  # Robust starting point by BIP AR-S approximation
  beta_initial <- robust_starting_point(x, p, q)$beta_initial
  beta_initial <- head(beta_initial, -1) # remove intercept
  result <- list()

  # objective function for ARMA model and BIP-ARMA model
  F <- function(beta) arma_s_resid_sc(x, beta, p, q)
  F_bip <- function(beta) bip_s_resid_sc(x, beta, p, q)[[1]]

  beta_arma <- lsqnonlin(F, -beta_initial, options = list(
    'tolx' = tolX))$x

  beta_bip <- lsqnonlin(F_bip, -beta_initial, options = list(
    'tolx' = tolX))$x

  # innovations m-scale for ARMA model
  a_sc <- arma_s_resid_sc(x, beta_arma, p, q)

  # innovations m-scale for BIP-ARMA model
  c(a_bip_sc, x_filt) %<-% bip_s_resid_sc(x, beta_bip, p, q)

  if( a_sc < a_bip_sc) beta_hat <- beta_arma else beta_hat <- beta_bip

  # final m-scale
  a_m_sc <- min(a_sc, a_bip_sc)

  # Output the results
  phi_bip_s <- c()
  phi_bip_s_init <- c()

  if(0 < p){
    phi_bip_s <- -beta_hat[1:p]
    phi_bip_s_init <- -beta_initial[1:p]
  }

  theta_bip_s <- c()
  theta_bip_s_init <- c()
  if(0 < q){
    theta_bip_s <- -beta_hat[(p+1):(p+q)]
    theta_bip_s_init <- -beta_hat[(p+1):(p+q)]
  }

  result$ar_coeffs <- phi_bip_s
  result$ma_coeffs <- theta_bip_s
  result$inno_scale <- a_m_sc
  result$cleaned_signal <- x_filt
  result$ar_coeffs_init <- phi_bip_s_init
  result$ma_coeffs_init <- theta_bip_s_init

  return(result)
} # end function
