#' bip_tau_arma_PSD_book
#'
#' @note
#'
#' File located in SpectrumEstimation_spec_arma_est_bip.R
#'
#' @export
bip_tau_arma_PSD_book <- function(x, p, q){
  x <- x - median(x)
  N <- length(x)
  w <- linspace(0, pi, N/2)
  s <- exp(1i * w)

  result <- arma_est_bip_tau(x, p, q)

  beta_hat <- c(result$ar_coeffs, result$ma_coeffs)

  Xx <- polyval(c(1, beta_hat[(p+1):(p+q)]) , s) / polyval(c(1, beta_hat[1:p]), s)
  sigma_hat <- result$inno_scale
  Pxx <- 1 / (2 * pi) * abs(Xx)^2
  PxxdB <- 10 * log10(Pxx)

  return(list('PxxdB' = PxxdB, 'Pxx' = Pxx, 'w' = w, 'sigma_hat' = sigma_hat))
}



#' spec_arma_est_bip_mm
#'
#' The function spec_arma_est_bip_mm(x,p,q) comuptes
#' spectral estimates using the BIP mm-estimates of the
#' ARMA model parameters
#'
#' @param x: numeric vector. The data
#' @param p: AR order
#' @param q: MA order
#'
#' @return PxxdB : Spectral estimate in dB
#' @return Pxx : spectral estimate
#' @return q : frequency in (0, pi)
#' @return sigma_hat : BIO mm-scale estemate of the innovations
#'
#' @references
#'
#' "Robust Statistics for Signal Processing"
#' Zoubir, A.M. and Koivunen, V. and Ollila, E. and Muma, M.
#' Cambridge University Press, 2018.
#'
#' "Bounded Influence Propagation \eqn{\tau}-Estimation: A New Robust Method for ARMA Model Estimation."
#' Muma, M. and Zoubir, A.M.
#' IEEE Transactions on Signal Processing, 65(7), 1712-1727, 2017.
#' @note
#'
#' File located in SpectrumEstimation_spec_arma_est_bip.R
#'
#' @export
spec_arma_est_bip_mm <- function(x, p, q, tolx = 5e-7){
  x <- x - median(x)
  N <- length(x)
  w <- linspace(0, pi, N/2)
  s <- exp(1i * w)

  result <- arma_est_bip_mm(x, p, q, tolX = tolx)

  beta_hat <- c(result$ar_coeffs, result$ma_coeffs)

  Xx <- polyval(c(1, beta_hat[(p+1):(p+q)]) , s) / polyval(c(1, beta_hat[1:p]), s)
  sigma_hat <- result$inno_scale
  Pxx <- sigma_hat^2 / (2 * pi) * abs(Xx)^2
  PxxdB <- 10 * log10(Pxx)

  return(list('PxxdB' = PxxdB, 'Pxx' = Pxx, 'w' = w, 'sigma_hat' = sigma_hat))
}


#' spec_arma_est_bip_tau
#'
#' The function spec_arma_est_bip_s(x,p,q) comuptes
#' spectral estimates using the BIP s-estimates of the
#' ARMA model parameters
#'
#' @param x: numeric vector. The data
#' @param p: AR order
#' @param q: MA order
#'
#' @return PxxdB : Spectral estimate in dB
#' @return Pxx : spectral estimate
#' @return q : frequency in (0, pi)
#' @return sigma_hat : BIO s-scale estemate of the innovations
#'
#' @references
#'
#' "Robust Statistics for Signal Processing"
#' Zoubir, A.M. and Koivunen, V. and Ollila, E. and Muma, M.
#' Cambridge University Press, 2018.
#'
#' "Bounded Influence Propagation \eqn{\tau}-Estimation: A New Robust Method for ARMA Model Estimation."
#' Muma, M. and Zoubir, A.M.
#' IEEE Transactions on Signal Processing, 65(7), 1712-1727, 2017.
#' @note
#'
#' File located in SpectrumEstimation_spec_arma_est_bip.R
#'
#' @export
spec_arma_est_bip_s <- function(x, p, q, tolx = 5e-7){
  x <- x - median(x)
  N <- length(x)
  w <- linspace(0, pi, N/2)
  s <- exp(1i * w)

  result <- arma_est_bip_s(x, p, q, tolX = tolx)

  beta_hat <- c(result$ar_coeffs, result$ma_coeffs)

  Xx <- polyval(c(1, beta_hat[(p+1):(p+q)]) , s) / polyval(c(1, beta_hat[1:p]), s)
  sigma_hat <- result$inno_scale
  Pxx <- sigma_hat^2 / (2 * pi) * abs(Xx)^2
  PxxdB <- 10 * log10(Pxx)

  return(list('PxxdB' = PxxdB, 'Pxx' = Pxx, 'w' = w, 'sigma_hat' = sigma_hat))
}


#' spec_arma_est_bip_tau
#'
#' The function spec_arma_est_bip_tau(x,p,q) comuptes
#' spectral estimates using the BIP tau-estimates of the
#' ARMA model parameters
#'
#' @param x: numeric vector. The data
#' @param p: AR order
#' @param q: MA order
#'
#' @return PxxdB : Spectral estimate in dB
#' @return Pxx : spectral estimate
#' @return q : frequency in (0, pi)
#' @return sigma_hat : BIO tau_scale estemate of the innovations
#'
#' @references
#'
#' "Robust Statistics for Signal Processing"
#' Zoubir, A.M. and Koivunen, V. and Ollila, E. and Muma, M.
#' Cambridge University Press, 2018.
#'
#' "Bounded Influence Propagation \eqn{\tau}-Estimation: A New Robust Method for ARMA Model Estimation."
#' Muma, M. and Zoubir, A.M.
#' IEEE Transactions on Signal Processing, 65(7), 1712-1727, 2017.
#' @note
#'
#' File located in SpectrumEstimation_spec_arma_est_bip.R
#'
#' @export
spec_arma_est_bip_tau <- function(x, p, q, tolx = 5e-7){
  x <- x - median(x)
  N <- length(x)
  w <- linspace(0, pi, N/2)
  s <- exp(1i * w)

  result <- arma_est_bip_tau(x, p, q, tolx = tolx)

  Xx <- polyval(c(1, result$ma_coeffs) , s) / polyval(c(1, result$ar_coeffs), s)
  sigma_hat <- result$inno_scale
  Pxx <- sigma_hat^2 / (2 * pi) * abs(Xx)^2
  PxxdB <- 10 * log10(Pxx)

  return(list('PxxdB' = PxxdB, 'Pxx' = Pxx, 'w' = w, 'sigma_hat' = sigma_hat))
}
