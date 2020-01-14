#' arma_est_bip_s
#' The function  arma_est_bip_mm(x,p,q) comuptes BIP S-estimates of the
#' ARMA model parameters.
#' It also computes an outlier cleaned signal using BIP-ARMA(p,q) predictions
#'
#' @param x: data (observations/measurements/signal)
#' @param p: autoregressive order
#' @param q: moving-average order
#' @param tolX: numeric. Threshold passed to pracma::lsqnonlin. Default = 1e-2
#'
#' @return result: named list with following fields
#'                 \item{ar_coeffs}{numeric vector of length p. BIP-AR(p) S estimates}
#'                 \item{ma_coeffs}{numeric vector of length q. BIP-AR(q) S estimates}
#'                 \item{inno_scale}{numeric, BIP s-estimate of the innovations scale}
#'                 \item{ar_coeffs_init}{numeric vector of length p. Robust starting point for estimation}
#'                 \item{ma_coeffs_init}{numeric vector of length q. Robust starting point for estimation}
#'
#'
#'
#' @references
#'
#'   "Robust Statistics for Signal Processing"
#'   Zoubir, A.M. and Koivunen, V. and Ollila, E. and Muma, M.
#'   Cambridge University Press, 2018.
#'
#'  "Bounded Influence Propagation $\tau$-Estimation: A New Robust Method for ARMA Model Estimation."
#'   Muma, M. and Zoubir, A.M.
#'   IEEE Transactions on Signal Processing, 65(7), 1712-1727, 2017.
#'
#' @examples
#' library(signal)
#'
#' N <- 500
#' a <- rnorm(N)
#' p <- 1
#' q <- 0
#' x <- signal::filter(1, c(1, -0.8), a)
#'
#' arma_est_bip_s(x, p, q)
#' @note
#'
#' file is in dependentData_armaEstBipS.R
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
