# arma_s_resid ----

#' arma_s_resid
#'
#' Computes the residuals of S-estimates for an ARMA process
#'
#' @param x: numeric vector. The signal
#' @param beta_hat: numeric vector, the AR and MA coefficients
#' @param p: AR order
#' @param q: MA order
#'
#' @return a_sc: m-scale based on residuals.
#'
#' @note
#' File is in dependentData_Auxiliary_armaResid.R
#'
#' @examples
#'
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
#' arma_s_resid(x, beta, p, q)
#'
#'@export
arma_s_resid <- function(x, beta_hat, p, q){
  # phi_hat := AR coefficients
  if(0 < p){phi_hat <- beta_hat[1:p]} else phi_hat <- numeric(0)
  if(0 < q){theta_hat <- beta_hat[(p+1):(p+q)]} else theta_hat <- numeric(0)

  N <- length(x)
  r <- max(p, q)

  a <- numeric(N)
  x_sc <- m_scale(x)

  if(r == 0) {
    a <- x
    a_sc <- x_sc}
  else{
    if(or(sum(abs(roots(c(1,-phi_hat))) > 1)
          , sum(abs(roots(c(1,-theta_hat))) > 1))) a_sc <- 10^10
    else{
    #ARMA
    if(p >= 1 && q >= 1){
      for(ii in (r+1):N){
        # ARMA residuals
        a[ii] <- x[ii] - phi_hat %*% x[seq.int(ii-1, ii-p, by = -1)] + theta_hat %*% a[seq.int(ii-1, ii-q, by = -1)]
      }
    }
    # MA model
    if(p == 0 && q >= 1){
      for(ii in (r+1):N){
        a[ii] <- x[ii] + theta_hat %*% a[seq.int(ii-1, ii-q, by = -1)]
      }
    }
    # AR model
    if(p >= 1 && q == 0){
      for(ii in (r+1):N){

        a[ii] <- x[ii] - phi_hat %*% x[seq.int(ii-1, ii-p, by = -1)]
      }
    }

      }
    }
  a_sc <- m_scale(a[(p+1):N])
  return(a)
}


# arma_tau_resid_sc ----
#' arma_tau_resid_sc
#'
#' @param x: numeric vector. The signal
#' @param beta_hat: numeric vector, the AR and MA coefficients
#' @param p: AR order
#' @param q: MA order
#'
#'
#' @return a: numeric vector of residuals
#'
#' @examples
#'
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
#' arma_tau_resid_sc(x, beta, p, q)
#'
#' @note
#' File is in dependentData_Auxiliary_armaResid.R
#'
#' @export
arma_tau_resid_sc <- function(x, beta_hat, p, q){
  # phi_hat := AR coefficients
  if(0 < p){phi_hat <- beta_hat[1:p]} else phi_hat <- numeric(0)
  if(0 < q){theta_hat <- beta_hat[(p+1):(p+q)]} else theta_hat <- numeric(0)

  N <- length(x)
  r <- max(p, q)

  a <- numeric(N)
  x_sc <- tau_scale(x)

  if(r == 0) {
    a <- x
    return(x_sc)}
  else{
    if(or(sum(abs(roots(c(1,-phi_hat))) > 1)
          , sum(abs(roots(c(1,-theta_hat))) > 1))) return(10^10)
    else{
      #ARMA
      if(p >= 1 && q >= 1){
        for(ii in (r+1):N){
          # ARMA residuals
          a[ii] <- x[ii] - phi_hat %*% x[seq.int(ii-1, ii-p, by = -1)] + theta_hat %*% a[seq.int(ii-1, ii-q, by = -1)]
        }
      }
      # MA model
      if(p == 0 && q >= 1){
        for(ii in (r+1):N){
          a[ii] <- x[ii] + theta_hat %*% a[seq.int(ii-1, ii-q, by = -1)]
        }
      }
      # AR model
      if(p >= 1 && q == 0){
        for(ii in (r+1):N){
          a[ii] <- x[ii] - phi_hat %*% x[seq.int(ii-1, ii-p, by = -1)]
        }
      }
    }
  }
  a_sc <- tau_scale(a[(p+1):N])
  return(a_sc)
}


# bip_resid ----
#' bounded influence propagation residuals
#'
#' @note
#' File is in dependentData_Auxiliary_armaResid.R
#'
#' @note
#' File in dependentData_Auxiliary_armaResid
#'
#'
#' @examples
#'
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
#' bip_resid(x, beta, p, q)
#' @export
bip_resid <- function(x, beta_hat, p, q){
  # phi_hat := AR coefficients
  if(0 < p){phi_hat <- beta_hat[1:p]} else phi_hat <- numeric(0)
  if(0 < q){theta_hat <- beta_hat[(p+1):(p+q)]} else theta_hat <- numeric(0)

  N <- length(x)
  r <- max(p, q)

  a_bip <- numeric(N)
  x_sc <- m_scale(x)

  kap2 <- 0.8724286 # kap=var(eta(randn(10000,1)));

  if(r == 0) {
    a_bip <- x
  }
  else{
    if(or(sum(abs(roots(c(1,-phi_hat))) > 1)
          , sum(abs(roots(c(1,theta_hat))) > 1))) {
      a_bip <- x
      sigma_hat <- x_sc
    }
    else{
      lamb <- ma_infinity(phi_hat, -theta_hat, 100);
      sigma_hat <- sqrt(x_sc^2/(1+kap2*sum(lamb^2)));

      #ARMA
      if(p >= 1 && q >= 1){
        for(ii in (r+1):N){
          # ARMA residuals
          a_bip[ii] <- x[ii] - phi_hat %*% (x[seq.int(ii-1, ii-p, by = -1)] - a_bip[seq.int(ii-1, ii-p, by = -1)] +
                                             sigma_hat *Rrobustsp::eta(a_bip[seq.int(ii-1, ii-p, by = -1)]/sigma_hat)) +
           (theta_hat * sigma_hat) %*% Rrobustsp::eta(a_bip[seq.int(ii-1, ii-q, by = -1)] / sigma_hat)
        }
      }
      # MA model
      if(p == 0 && q >= 1){
        for(ii in (r+1):N){
          a_bip[ii] <- x[ii] + (theta_hat * sigma_hat) %*% Rrobustsp::eta(a_bip[seq.int(ii-1, ii-q, by = -1)] / sigma_hat)
        }
      }
      # AR model
      if(p >= 1 && q == 0){
        for(ii in (r+1):N){
        a_bip[ii] <- x[ii] - phi_hat %*% (x[seq.int(ii-1, ii-p, by = -1)] - a_bip[seq.int(ii-1, ii-p, by = -1)] +
                                           sigma_hat * Rrobustsp::eta(a_bip[seq.int(ii-1, ii-p, by = -1)]/sigma_hat))
      }}

    }
  }

  a_bip <- a_bip[(p + 1): length(a_bip)]
  return(a_bip)
}


# bip_s_resid ----
#' bip_s_resid
#'
#' @note
#' file in dependentData_Auxiliary_armaResid.R
#'
#' @examples
#'
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
#' bip_s_resid(x, beta, p, q)
#' @export
bip_s_resid <- function(x, beta_hat, p, q){
  # phi_hat := AR coefficients
  if(0 < p){phi_hat <- beta_hat[1:p]} else phi_hat <- numeric(0)
  if(0 < q){theta_hat <- beta_hat[(p+1):(p+q)]} else theta_hat <- numeric(0)

  N <- length(x)
  r <- max(p, q)

  a_bip <- numeric(N)
  x_sc <- m_scale(x)

  kap2 <- 0.8724286 # kap=var(eta(randn(10000,1)));


  if(sum(abs(roots(c(1,-phi_hat))) > 1) | sum(abs(roots(c(1,-theta_hat))) > 1)) {
    sigma_hat <- x_sc
  }
  else{
    lamb <- ma_infinity(phi_hat, -theta_hat, 100) # MA infinity approximation to compute scale used in eta function
    sigma_hat <- sqrt(x_sc^2/(1+kap2*sum(lamb^2))); # scale used in eta function
  }


  if(r == 0) {
    return(list(x, numeric(0), x_sc))
    # a_bip <- x
    # a_sc_bip <- x_sc
  }
  else{
    if(sum(abs(roots(c(1,-phi_hat))) > 1) | sum(abs(roots(c(1,-theta_hat))) > 1)) {
      a_bip_sc <- 10^10
    }
    else{
      #ARMA
      if(p >= 1 && q >= 1){
        for(ii in (r+1):N){
          # ARMA residuals
          a_bip[ii] <- x[ii] - phi_hat %*% (x[seq.int(ii-1, ii-p, by = -1)] - a_bip[seq.int(ii-1, ii-p, by = -1)] +
                                              sigma_hat * Rrobustsp::eta(a_bip[seq.int(ii-1, ii-p, by = -1)]/sigma_hat)) +
            (theta_hat * sigma_hat) %*% Rrobustsp::eta(a_bip[seq.int(ii-1, ii-q, by = -1)] / sigma_hat)
        }
      }
      # MA model
      if(p == 0 && q >= 1){
        for(ii in (r+1):N){
          a_bip[ii] <- x[ii] + (theta_hat * sigma_hat) %*% Rrobustsp::eta(a_bip[seq.int(ii-1, ii-q, by = -1)] / sigma_hat)
        }
      }
      # AR model
      if(p >= 1 && q == 0){
        for(ii in (r+1):N){
          a_bip[ii] <- x[ii] - phi_hat %*% (x[seq.int(ii-1, ii-p, by = -1)] - a_bip[seq.int(ii-1, ii-p, by = -1)] +
                                            sigma_hat * Rrobustsp::eta(a_bip[seq.int(ii-1, ii-p, by = -1)]/sigma_hat))
        }
      }

      a_bip_sc <- m_scale(a_bip[(p + 1): length(a_bip)])

      x_filt <- x

      for(ii in (p+1):N){
        x_filt[ii] <- x[ii] - a_bip[ii] + sigma_hat * Rrobustsp::eta(a_bip[ii]/sigma_hat)
      }
      return(list(a_bip, 'filtered_signal' = x_filt, 'scale' = a_bip_sc))

    }
  }
}


# bip_s_resid_sc ----
#' bip_s_resid_sc
#'
#' @note
#' file in dependentData_Auxiliary_armaResid.R
#'
#' @examples
#'
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
#' bip_s_resid_sc(x, beta, p, q)
#' @export
bip_s_resid_sc <- function(x, beta_hat, p, q){
  # phi_hat := AR coefficients
  if(0 < p){phi_hat <- beta_hat[1:p]} else phi_hat <- numeric(0)
  if(0 < q){theta_hat <- beta_hat[(p+1):(p+q)]} else theta_hat <- numeric(0)

  N <- length(x)
  r <- max(p, q)

  a_bip <- numeric(N)
  x_sc <- m_scale(x)

  kap2 <- 0.8724286 # kap=var(eta(randn(10000,1)));


  if(sum(abs(roots(c(1,-phi_hat))) > 1) | sum(abs(roots(c(1,-theta_hat))) > 1)) {
    sigma_hat <- x_sc
  }
  else{
    lamb <- ma_infinity(phi_hat, -theta_hat, 100) # MA infinity approximation to compute scale used in eta function
    sigma_hat <- sqrt(x_sc^2/(1+kap2*sum(lamb^2))) # scale used in eta function
  }


  if(r == 0) {
    return(list(x_sc, numeric(0)))
    # a_bip <- x
    # a_sc_bip <- x_sc
  }
  else{
    if(sum(abs(roots(c(1,-phi_hat))) > 1) | sum(abs(roots(c(1,-theta_hat))) > 1)) {
      a_bip_sc <- 10^10
      return(list(a_bip_sc, numeric(0)))
    }
    else{
      #ARMA
      if(p >= 1 && q >= 1){
        for(ii in (r+1):N){
          # ARMA residuals
          a_bip[ii] <- x[ii] - phi_hat %*% (x[seq.int(ii-1, ii-p, by = -1)] - a_bip[seq.int(ii-1, ii-p, by = -1)] +
                                              sigma_hat * Rrobustsp::eta(a_bip[seq.int(ii-1, ii-p, by = -1)]/sigma_hat)) +
            (theta_hat * sigma_hat) %*% Rrobustsp::eta(a_bip[seq.int(ii-1, ii-q, by = -1)] / sigma_hat)
        }
      }
      # MA model
      if(p == 0 && q >= 1){
        for(ii in (r+1):N){
          a_bip[ii] <- x[ii] + (theta_hat * sigma_hat) %*% Rrobustsp::eta(a_bip[seq.int(ii-1, ii-q, by = -1)] / sigma_hat)
        }
      }
      # AR model
      if(p >= 1 && q == 0){
        for(ii in (r+1):N){
          a_bip[ii] <- x[ii] - phi_hat %*% (x[seq.int(ii-1, ii-p, by = -1)] - a_bip[seq.int(ii-1, ii-p, by = -1)] +
                                              sigma_hat * Rrobustsp::eta(a_bip[seq.int(ii-1, ii-p, by = -1)]/sigma_hat))
        }
      }

      a_bip_sc <- m_scale(a_bip[(p + 1): length(a_bip)])

      x_filt <- x

      for(ii in (p+1):N){
        x_filt[ii] <- x[ii] - a_bip[ii] + sigma_hat * Rrobustsp::eta(a_bip[ii]/sigma_hat)
      }

      return(list('scale' = a_bip_sc, 'filtered_signal' = x_filt))

    }
  }
}


# bip_tau_resid_sc ----

#' bip_tau_resid_sc
#'
#' Computes the tau-estimate residual scale
#'
#' @param x : numeric vector. The signal
#' @param beta_hat : ARMA parameter estimates
#' @param p : AR order
#' @param q : MA order
#'
#' @return scale : the tau estimate residual scale
#' @return filtered_signal : Filtered signal
#'
#' @note
#' In file dependentData_Auxiliary_armaResid
#'
#' @examples
#' data("x_ao_mat")
#'
#' # 4.2918
#' res <- tau_scale(x_ao_mat)
#'
#' ## Example 2
#'
#' data("PPG_Data")
#'
#' library(zeallot)
#'
#' ibi_ppg <- diff(PPG_Data$ppg.pos[1,])
#'
#' x <- ibi_ppg[length(ibi_ppg):1] - median(ibi_ppg)
#'
#' beta_mat <- c(-0.6640 ,  -0.0815 ,  -0.2180 ,  -0.4591 ,  -0.0830 ,   0.2936 ,  -0.0434   , 0.0190 ,   0.0757 ,  -0.0938,   -0.1758)
#'
#' bip_tau_resid_sc(x, beta_mat, 0, 11)
#'
#' @export
bip_tau_resid_sc <- function(x, beta_hat, p, q){
  # phi_hat := AR coefficients
  if(0 < p){phi_hat <- beta_hat[1:p]} else phi_hat <- 0
  if(0 < q){theta_hat <- beta_hat[(p+1):(p+q)]} else theta_hat <- 0

  N <- length(x)
  r <- max(p, q)

  a_bip <- numeric(N)
  x_sc <- tau_scale(x)

  kap2 <- 0.8724286 # kap=var(eta(randn(10000,1)));


  if(sum(abs(roots(c(1,-phi_hat))) > 1) || sum(abs(roots(c(1,-theta_hat))) > 1)) {
    sigma_hat <- x_sc
  }
  else{
    lamb <- ma_infinity(phi_hat, -theta_hat, 100) # MA infinity approximation to compute scale used in eta function
    sigma_hat <- sqrt(x_sc^2/(1+kap2*sum(lamb^2))); # scale used in eta function
  }


  if(r == 0) {
    return(list(x_sc, x))
    # a_bip <- x
    # a_sc_bip <- x_sc
  }
  else{
    if(sum(abs(roots(c(1,-phi_hat))) > 1) || sum(abs(roots(c(1,-theta_hat))) > 1)) {
      a_bip_sc <- 10^10
      return(list(a_bip_sc, numeric(0)))
    }
    else{
      #ARMA
      if(p >= 1 && q >= 1){
        for(ii in (r+1):N){
          # ARMA residuals
          a_bip[ii] <- x[ii] - phi_hat %*% (x[seq.int(ii-1, ii-p, by = -1)] - a_bip[seq.int(ii-1, ii-p, by = -1)] +
                                              sigma_hat * Rrobustsp::eta(a_bip[seq.int(ii-1, ii-p, by = -1)]/sigma_hat)) +
            (theta_hat * sigma_hat) %*% Rrobustsp::eta(a_bip[seq.int(ii-1, ii-q, by = -1)] / sigma_hat)
        }
      }
      # MA model
      if(p == 0 && q >= 1){
        for(ii in (r+1):N){
          a_bip[ii] <- x[ii] + (theta_hat * sigma_hat) %*% Rrobustsp::eta(a_bip[seq.int(ii-1, ii-q, by = -1)] / sigma_hat)
        }
      }
      # AR model
      if(p >= 1 && q == 0){
        for(ii in (r+1):N){
          a_bip[ii] <- x[ii] - phi_hat %*% (x[seq.int(ii-1, ii-p, by = -1)] - a_bip[seq.int(ii-1, ii-p, by = -1)] +
                                              sigma_hat * Rrobustsp::eta(a_bip[seq.int(ii-1, ii-p, by = -1)]/sigma_hat))
        }
      }

      a_bip_sc <- tau_scale(a_bip[(p + 1): length(a_bip)])

      x_filt <- x

      for(ii in (p+1):N){
        x_filt[ii] <- x[ii] - a_bip[ii] + sigma_hat * Rrobustsp::eta(a_bip[ii]/sigma_hat)
      }
      return(list('scale' = a_bip_sc, 'filtered_signal' = x_filt))

    }
  }
}

