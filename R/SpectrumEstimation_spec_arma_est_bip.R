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
#' @note
#'
#' File located in SpectrumEstimation_spec_arma_est_bip.R
#'
#' @export
spec_arma_est_bip_mm <- function(x, p, q){
  x <- x - median(x)
  N <- length(x)
  w <- linspace(0, pi, N/2)
  s <- exp(1i * w)

  result <- arma_est_bip_mm(x, p, q)

  beta_hat <- c(result$ar_coeffs, result$ma_coeffs)

  Xx <- polyval(c(1, beta_hat[(p+1):(p+q)]) , s) / polyval(c(1, beta_hat[1:p]), s)
  sigma_hat <- result$inno_scale
  Pxx <- sigma_hat^2 / (2 * pi) * abs(Xx)^2
  PxxdB <- 10 * log10(Pxx)

  return(list('PxxdB' = PxxdB, 'Pxx' = Pxx, 'w' = w, 'sigma_hat' = sigma_hat))
}


#' spec_arma_est_bip_s
#'
#' @note
#'
#' File located in SpectrumEstimation_spec_arma_est_bip.R
#'
#' @export
spec_arma_est_bip_s <- function(x, p, q){
  x <- x - median(x)
  N <- length(x)
  w <- linspace(0, pi, N/2)
  s <- exp(1i * w)

  result <- arma_est_bip_s(x, p, q)

  beta_hat <- c(result$ar_coeffs, result$ma_coeffs)

  Xx <- polyval(c(1, beta_hat[(p+1):(p+q)]) , s) / polyval(c(1, beta_hat[1:p]), s)
  sigma_hat <- result$inno_scale
  Pxx <- sigma_hat^2 / (2 * pi) * abs(Xx)^2
  PxxdB <- 10 * log10(Pxx)

  return(list('PxxdB' = PxxdB, 'Pxx' = Pxx, 'w' = w, 'sigma_hat' = sigma_hat))
}


#' spec_arma_est_bip_tau
#'
#' @note
#'
#' File located in SpectrumEstimation_spec_arma_est_bip.R
#'
#' @export
spec_arma_est_bip_tau <- function(x, p, q){
  x <- x - median(x)
  N <- length(x)
  w <- linspace(0, pi, N/2)
  s <- exp(1i * w)

  result <- arma_est_bip_tau(x, p, q)

  beta_hat <- c(result$ar_coeffs, result$ma_coeffs)

  Xx <- polyval(c(1, beta_hat[(p+1):(p+q)]) , s) / polyval(c(1, beta_hat[1:p]), s)
  sigma_hat <- result$inno_scale
  Pxx <- sigma_hat^2 / (2 * pi) * abs(Xx)^2
  PxxdB <- 10 * log10(Pxx)

  return(list('PxxdB' = PxxdB, 'Pxx' = Pxx, 'w' = w, 'sigma_hat' = sigma_hat))
}
