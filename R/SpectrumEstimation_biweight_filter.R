#' biweight_filter
#'
#' The biweight_filter(x) is our implementation of the
#' method described in
#'
#' "High breakdown methods of time series analysis.
#' Tatum, L.G., and Hurvich, C. M.
#' Journal of the Royal Statistical Society. Series B (Methodological),
#' pp. 881-896, 1993.
#'
#' The code is based on an implementation by Falco Strasser, Signal Processing
#' Group, TU Darmstadt, October 2010.
#'
#'
#' @param x: data (observations/measurements/signal), real-valued vector
#'
#' @return xFBi: Biweight filtered (outlier cleaned) signal
#' @return ABi:  Fourier coefficients for cosine
#' @return BBi:  Fourier coefficients for sine
#'
#'
#' @references
#' "Robust Statistics for Signal Processing"
#' Zoubir, A.M. and Koivunen, V. and Ollila, E. and Muma, M.
#' Cambridge University Press, 2018.
#'
#' @note
#' ATM takes forever
#' File location: SpectrumEstimation_biweight_filter.R
#'
#' @export
biweight_filter <- function(x){
  # to avoid no visible binding for global variable note
  xFRM <- ARM <- BRM <- nls <- N_Prime <- xFb2
  N <- length(x)
  xFB <- numeric(N) # filter cleaned signal

  # Works on signals of prime length, thus signal is split is into
  # 2 overlapping signals of prime length
  x_split <- split_into_prime(x)

  # Biweight filter each prime segment
  for(ii in 1:2){
    x_part <- x_split[,ii]
    wr <- order_wk(x_part) # Fourier frequencies in descending order
    c(xFRM, ARM, BRM) %<-% repeated_median_filter(x_part)
    N_prime <- length(x_part) # length of the prime segment
    t <- 0:(N_prime - 1)

    K <- (N_prime - 1) / 2 # number of Fourier coefficients
    ABi <- numeric(k) # biweight estimate of cosine coefficients at w(k)
    BBi <- numeric(k) # biweight estimate of sine coefficients at w(k)

    k <- 4 # tuning constant as recommended in the paper by Tatum and Hurvich
    xb <- MlocTUK(x_part, k) # Tukeys location M-estimate
    xc <- (x_part - xb) # robustly centered time series

    for(k in 1:K){
      c <- nls(xc ~ a*cos(wr[k]*t)+b*sin(wr[k]*t), start = c(a = ARM[k], b= BRM[k]))
      ABi[k] <- c$a
      BBi[k] <- c$b
      xc <- xc - ABi[k]*cos(wr[k]*t)+BBi[k]*sin(wr[k]*t)
    }

    # recover the core process by regression of the Biweight estimates onto the
    # independent parameters
    if(ii == 1){
      xFB1 <- xb
      for(k in 1:K){
        sumAB <- ABi[k]*cos(wr[k]*t) + BBi[k]*sin(wr[k]*t)
        xFB1 <- xFB1 + sumAB
      }
    }
    if(ii==2){
      xFB2 <- xb
      for(k in 1:K){
        sumAB <- ABi[k]*cos(wr[k]*t) + BBi[k]*sin(wr[k]*t)
        xFB2 <- xFB2 + sumAB
      }
    }


  }

  # fuse the cleaned segments of prime length
  if(ii==1) xFB <- xFB1
  else{
    xFB[1:(N - N_prime)] <- xFB1[1:(N - N_prime)]
    xFB[(N-N_Prime+1):N_prime] <- 0.5 * (xFB1[(N-N_Prime+1):N_prime] + xFb2[(N-N_Prime+1):N_prime])
    xFB[(N_prime+1):N] <- xFB2[(N_prime+1):N]
  }

  return(list('xFB' = xFB, 'ABi' = ABi, 'BBi' = BBi))


}
