#' arma_est_bip_tau
#'
#' The function arma_est_bip_tau(x,p,q) comuptes BIP tau-estimates of the
#' ARMA model parameters. It also computes an outlier
#' cleaned signal using BIP-ARMA(p,q) predictions
#'
#' @param x: numeric vector. The signal
#' @param p: AR order
#' @param q: MA order
#' @param tolx: threshold value that is passed to pracma::lsqnonlin.
#'              Default = 5e-7
#'
#' @return result: named list with following fields
#'                 ar_coeffs : numeric vector of length p. BIP-AR(p) tau estimates
#'                 ma_coeffs : numeric vector of length q. BIP-AR(q) tau estimates
#'                 inno_scale : numeric, BIP s-estimate of the innovations scale
#'                 ar_coeffs_init : numeric vector of length p. Robust starting point for estimation
#'                 ma_coeffs_init : numeric vector of length q. Robust starting point for estimation
#'
#'
#' @examples
#'
#' N <- 500
#' a <- rnorm(N)
#' p <- 1
#' q <- 0
#' x <- signal::filter(1, c(1, -0.8), a)
#'
#' arma_est_bip_tau(x, p, q)
#'
#' @note
#' file is in dependentData_arma_Est_BipTau.R
#'
#' @export
arma_est_bip_tau <- function(x, p, q, tolx = 5e-7){
  result <- list()

  if(p == 0 && q == 0){
    result$inno_scale < tau_scale(x)
    warning('Please choose a nonzero value for p or q')
    return(result)}
  if(length(x) <= (p + q)){
    warning('There are too many parameters to estimate for chosen data size. Reduce model order or use a larger data set.')
    return(result)
  }

  # Robust starting point by BIP AR-S approximation
  beta_initial <- robust_starting_point(x, p, q)[[1]]
  # remove intercept from coefficients
  beta_initial <- head(beta_initial, -1)

  # objective function for ARMA model and BIP-ARMA model
  F <- function(beta) arma_tau_resid_sc(x, beta, p, q)
  F_bip <- function(beta) bip_tau_resid_sc(x, beta, p, q)[[1]]

  beta_arma <- lsqnonlin(F, -beta_initial, options = list(tolx = tolx))$x

  beta_bip <- lsqnonlin(F_bip, -beta_initial)$x

  # innovations tau-scale for ARMA model
  a_sc <- arma_tau_resid_sc(x, beta_arma, p, q)[[1]]

  # innovations m-scale for BIP-ARMA model
  c(a_bip_sc, x_filt) %<-% bip_tau_resid_sc(x, beta_bip, p, q)

  if( a_sc < a_bip_sc) beta_hat <- beta_arma else beta_hat <- beta_bip

  # final m-scale
  a_tau_sc <- min(a_sc, a_bip_sc)

  # Output the results
  phi_bip_tau <- c()
  phi_bip_tau_init <- c()

  if(0 < p){
    phi_bip_tau <- -beta_hat[1:p]
    phi_bip_tau_init <- -beta_initial[1:p]
  }

  theta_bip_tau <- c()
  theta_bip_tau_init <- c()
  if(0 < q){
    theta_bip_tau <- -beta_hat[(p+1):(p+q)]
    theta_bip_tau_init <- -beta_hat[(p+1):(p+q)]
  }

  result$ar_coeffs <- phi_bip_tau
  result$ma_coeffs <- theta_bip_tau
  result$inno_scale <- a_tau_sc
  result$cleaned_signal <- x_filt
  result$ar_coeffs_init <- phi_bip_tau_init
  result$ma_coeffs_init <- theta_bip_tau_init
  return(result)
} # end function
