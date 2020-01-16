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
#'  "Bounded Influence Propagation \eqn{\tau}-Estimation: A New Robust Method for ARMA Model Estimation."
#'   Muma, M. and Zoubir, A.M.
#'   IEEE Transactions on Signal Processing, 65(7), 1712-1727, 2017.
#'
#' @examples
#' library(signal)
#' library(zeallot)
#' library(pracma)
#'
#' N <- 500
#' a <- rnorm(N)
#' p <- 1
#' q <- 0
#' x <- signal::filter(1, c(1, -0.8), a)
#'
#' beta_s <- arma_est_bip_s(x, p, q, tolX = 1e-8)
#' beta <- c(beta_s$ar_coeffs, beta_s$ma_coeffs)
#'
#' arma_est_bip_m(x, p, q, beta, beta_s$inno_scale)
#'@note
#'
#'File location: dependentData_armaEstBipM.R
#'
#' @export
arma_est_bip_m <- function(x, p, q, beta_hat_s, a_sc_final){
  N <- length(x)

  F_mm <- function(beta) (sqrt(1/(N-p)*sum(muler_rho2(arma_resid(
    x/a_sc_final, beta, p, q)))))

  F_bip_mm <- function(beta) (sqrt(1/(N-p)*sum(muler_rho2(bip_resid(
    x/a_sc_final, beta, p, q))/a_sc_final)))

  beta_arma_mm <- lsqnonlin(F_mm, beta_hat_s)$x
  beta_bip_mm <- lsqnonlin(F_bip_mm, beta_hat_s)$x

  a_rho2_mm <- 1 / (N - p) * sum(muler_rho2(arma_resid(x/a_sc_final,
                                                       beta_arma_mm,
                                                       p, q)))

  a_bip_rho2_mm <- 1 / (N - p) * sum(muler_rho2(bip_resid(x/a_sc_final,
                                                           beta_bip_mm,
                                                           p, q)))

  if(a_rho2_mm < a_bip_rho2_mm) beta_hat <- beta_arma_mm
  else beta_hat <- beta_bip_mm

  if(0 < p) phi_bip_mm <- - beta_hat[1:p] else phi_bip_mm <- numeric(0)

  if(0 < q) theta_bip_mm <- - beta_hat[(p+1):(p+q)] else theta_bip_mm <- numeric(0)

  return(list('ar_coeffs' = phi_bip_mm,
              'ma_coeffs' = theta_bip_mm))
  }
