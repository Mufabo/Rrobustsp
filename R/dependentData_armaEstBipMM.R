#' arma_est_bip_m
#'
#'  The function  arma_est_bip_m(x,p,q) comuptes the BIP M-estimation step for BIP MM estimates of the
#'  ARMA model parameters. It can also be used as a stand-alone
#'  M-estimator.
#'
#'
#' @param x: data (observations/measurements/signal)
#' @param p: autoregressive order
#' @param q: moving-average order
#' @param beta_hat_s: BIP S-estimate
#' @param a_sc_final: M scale estimate of residuals of BIP S-estimate
#'
#'
#' @return ar_coeffs: vector of BIP-AR(p) MM-estimates
#' @return ma_coeffs: vector of BIP-MA(q) MM-estimates
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
#'
#' @export
arma_est_bip_mm <- function(x, p, q, tolX = 5e-7){
  # Get starting point and residual scale from S-estimtor
  bip_s_est = arma_est_bip_s(x, p, q, tolX)
  beta_hat_s = c(bip_s_est$ar_coeffs, bip_s_est$ma_coeffs)

  N <- length(x)

  F_mm <- function(beta){
    (1/(N-p) * muler_rho2(arma_s_resid(x, beta, p, q)) / bip_s_est$inno_scale)
  }

  F_bip_mm <- function(beta){
    (1/(N-p) *
       muler_rho2(bip_s_resid(x, beta, p, q)[[1]]) / bip_s_est$inno_scale)
  }

  beta_arma_mm <- lsqnonlin(F_mm, -beta_hat_s, options = list('tolx' = tolX))$x
  beta_bip_mm <- lsqnonlin(F_bip_mm, -beta_hat_s, options = list('tolx' = tolX))$x

  a <- arma_s_resid(x, beta_arma_mm, p, q)
  a_sc <- m_scale(a)

  a_bip <- bip_s_resid(x, beta_bip_mm, p, q)[[1]]
  a_bip_sc <- m_scale(a_bip)

  if(a_sc < a_bip_sc) beta_hat <- beta_arma_mm
  else beta_hat <- beta_bip_mm

  a_m_sc <- min(a_sc, a_bip_sc)

  if(0 < p) phi_bip_mm <- -beta_hat[1:p] else phi_bip_mm <- numeric(0)

  if(0 < q) theta_bip_mm <- -beta_hat[(p+1):(p+q)] else theta_bip_mm <- numeric(0)

  result <- list()
  result$ar_coeffs <- phi_bip_mm
  result$ma_coeffs <- theta_bip_mm
  result$inno_scale <- a_m_sc
  result$cleaned_signal <- bip_s_est$cleaned_signal
  result$ar_coeffs_init <- bip_s_est$ar_coeffs
  result$ma_coeffs_init <- bip_s_est$ma_coeffs

  return(result)
}
