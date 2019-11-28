#'  Mloc computes the M-estimates of location using an auxiliary scale
#'  estimate. It uses the iterative reweighted least squares (IRWLS) algorithm
#'
#'
#'@param         y : (numeric) data vector of size N x 1
#'@param   lossfun : (string) either 'huber' or 'tukey' to identify the desired
#'             loss function
#'
#'@return             mu_hat: M-estimate of location
#' @export
Mloc <- function(y, lossfun){
  return (Mreg(y, numeric(length(y)), lossfun, median(y)))
}

#'   Mloc_HUB computes Huber's M-estimate of
#'   location, i.e.,
#'
#'   \eqn{mu_hat = arg_{min_mu} \sum{i} rho_{HUB}(y_{i} - \mu)}
#'
#'
#'
#'@param            y: real valued data vector of size N x 1
#'@param             c: tuning constant c>=0
#'
#'
#'@return             mu_hat: Huber's M-estimate of location
#' @export
MlocHUB <- function(y,c=1.345, max_iters=1000, tol_err=1e-5){
  # previously computated scale estimate
  sigma_0 <- madn(y)

  # initial robust location estimate
  mu_n <- median(y)

  for (n in 0:max_iters){
    w_n <- whub(abs(y-mu_n)/sigma_0, c) # compute weights

    mu_n_plus1 <- sum(w_n*y)/sum(w_n) # compute weighted average

    if (abs(mu_n_plus1-mu_n)/sigma_0 > tol_err) mu_n <- mu_n_plus1 else break
  }

  return (mu_n)
}

#'   Mlocscale computes Huber's M-estimates of location and scale.
#'
#'
#'@param             y: real valued data vector of size N x 1
#'@param             c: tuning constant c>=0
#'
#'
#'@return             mu_hat: Huber's M-estimate of location
#'@return             sigma_hat: Huber's M-estimate of scale
#' @export
Mlocscale <- function(y, c=NULL){
  # approx 95 efficiency for Gaussian errors
  if(is.null(c)){
    if (is.complex(y)) c<-1.215 else c<-1.3415
  }
  return (hubreg(y, rep(1, length(y)), c, madn(y), median(y)))
}

#'  Mloc_TUK computes Tukey's M-estimate of
#'   location, i.e.,
#'
#'   mu_hat = arg min_mu SUM_i rho_TUK(y_i - mu)
#'
#'
#'
#'@param             y: real valued data vector of size N x 1
#'@param             c: tuning constant c>=0
#'
#'
#'@return             mu_hat: Tukey's M-estimate of location
#'   @export
MlocTUK <- function(y, c=4.685, max_iters=1000, tol_err=1e-5){
  # previously computed scale estimate
  sigma_0 <- madn(y)

  # initial robust location estimate
  mu_n <- median(y)

  for (n in 0:max_iters){
    w_n <- wtuk(abs(y-mu_n)/sigma_0, c) # compute weights

    mu_n_plus1 <- sum(w_n*y)/sum(w_n) # compute weighted average

    if (abs(mu_n_plus1-mu_n)/sigma_0 > tol_err) mu_n <- mu_n_plus1 else break
  }

  return (mu_n)
}

#'   Mscale_HUB computes Huber's M-estimate of
#'   scale.
#'
#'
#'
#'@param             y: real valued data vector of size N x 1
#'@param             c: tuning constant c>=0
#'
#'
#'@return             sigma_hat: Huber's M-estimate of scale
#' @export
MscaleHUB <- function(y, c=1.345, max_iters=1000, tol_err=1e-5){
  # initial scale estimate
  sigma_n <- madn(y)

  # subtract previously computed location
  y <- y-median(y)


  N <- length(y)

  # consistency with the standard deviation at the Gaussian
  # set.seed(1)
  # delta <- mean(rhohub(rnorm(10000), c))
  delta <- 0.9379 # delta from matlab
  for (n in 0:max_iters){
    w_n <- whub(abs(y)/sigma_n, c)

    sigma_n_plus1 <- sqrt(1/(N*delta)*sum(w_n*y^2))
    if (abs(sigma_n_plus1/sigma_n-1)>tol_err)
      {sigma_n <- sigma_n_plus1}
    else break
  }

  return (sigma_n)
}


#'   Mscale_TUK computes Tukey's M-estimate of
#'   scale.
#'
#'
#'
#'@param             y: real valued data vector of size N x 1
#'@param             c: tuning constant c>=0
#'
#'
#'@return             sigma_hat: Tukey's M-estimate of scale
#' @export
MscaleTUK <- function(y, c=4.685, max_iters=1000, tol_err=1e-5){
  # initial scale estimate
  sigma_n <- madn(y)

  # subtract previously computed location
  y <- y-median(y)


  N <- length(y)

  # consistency with the standard deviation at the Gaussian
  # set.seed(1)
  # delta <- mean(rhotuk(rnorm(10000), c))
  delta <- 0.8788 # delta from matlab
  for (n in 0:max_iters){
    w_n <- wtuk(abs(y)/sigma_n, c)

    sigma_n_plus1 <- sqrt(1/(N*delta)*sum(w_n*y^2))
    if (abs(sigma_n_plus1/sigma_n-1)>tol_err) sigma_n <- sigma_n_plus1 else break
  }

  return (sigma_n)
}
