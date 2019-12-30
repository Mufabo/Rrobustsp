#' Mscat
#'
#' computes M-estimator of scatter matrix for the n x p data matrix X
#' using the loss function 'Huber' or 't-loss' and for a given parameter of
#' the loss function (i.e., q for Huber's or degrees of freedom v for
#' the t-distribution).
#'
#' @param x : n x p matrix
#' @param loss : 'Huber', 't_loss' or 'Tyler'
#' @param losspar: parameter of the loss function: q in [0,1) for Huber and
#'          d.o.f. v >= 0 for t-loss. For Tyler you do not need to specify
#'          this value. Parameter q determines the treshold
#'          c^2 as the qth quantile of chi-squared distribution with p
#'          degrees of freedom distribution (Default q = 0.8). Parameter v
#'          is the def.freedom of t-distribution (Default v = 3)
#'          if v = 0, then one computes Tyler's M-estimator
#' @param    invC: initial estimate is the inverse scatter matrix (default =
#'           inverse of the sample covariance matrix)
#' @param  printitn: integer, print iteration number (default = 0, no printing)
#'
#' @return       C: the M-estimate of scatter using Huber's weights
#' @return    invC: the inverse of C
#' @return      iter: nr of iterations
#' @return     flag: flag (true/false) for convergence
#'
#' @examples
#'
#' @export
Mscat <- function(x, loss, losspar = NULL, invC = NULL, printitn = 0){

  if(! loss %in% c('Huber', 't_loss', 'Tyler')){
    sprintf('invalid loss. Choose one of Huber, t_loss, Tyler. Your loss is %s', loss)
  }
  if(is.complex(losspar) | is.infinite(losspar) | (0>losspar) ){
    if(loss == 'Huber' && (losspar>1)){
      sprintf('Argument losspar needs to be a scalar in [0, 1[')
      return(list(NULL, NULL, NULL, NULL))
    }
    if(loss == 't_loss'){
      sprintf('Argument losspar needs to be a scalar greater zero')
      return(list(NULL, NULL, NULL, NULL))}
  }

  n <- nrow(x)
  p <- ncol(x)

  if(is.complex(x)){
    realdata <- F
    if(is.null(invC)){
      C <- Conj(t(x)) %*% x / n
      invC <- solve(C, diag(1, p, p), tol = 1e-30)
    }
  } else {
    realdata <- T
    if(is.null(invC)){
      C <- t(x) %*% x / n
      invC <- solve(C, diag(1, p, p), tol = 1e-30)
    }
  }

  # loss functions below
  switch(loss,
         Huber = tmp <- loss_Huber(n, p, losspar, realdata),
         Tyler = tmp <- loss_t_loss(n, p, losspar, realdata),
         t_loss= tmp <- loss_t_loss(n, p, losspar, realdata)
         )

  MAX_ITER <- 1000
  EPS <- 1e-5

  flag <- F
  for(iter in 1:MAX_ITER){
    t <- Re(rowSums(x %*% invC * Conj(x)))

    C <- tmp$const * Conj(t(x)) %*% (x * tmp$ufun(t, tmp$upar))
    d <- inf_norm(diag(1, p, p) - invC %*% C)

    if(printitn && iter %% printitn == 0) sprintf('At iter = %4d, dis=%.6f\n',iter,d)

    invC <- solve(C, diag(1, p, p), tol = 1e-30)

    if(d <= EPS){
      flag <- T
      break
    }
  }
  if(iter == MAX_ITER) print('Warning: slow convergence')
  return(list('C' = C, 'invC' = invC, 'iter' = iter, 'flag' = flag))
}

loss_Huber <- function(n, p, losspar, realdata){
  if(is.null(losspar)) q <- 0.9 else q <- losspar

  ufun <- function(t, c) (t <= c) + (c/t) * (t > c)

  if(realdata){
    upar <- qchisq(q, p)
    b <- pchisq(upar, p+2) + (upar / p) * (1 - q)
  } else {
    upar <- qchisq(q, 2 * p) / 2
    b <- pchisq(2 * upar, 2 * (p + 1)) + (upar / p) * (1 - q)
  }

  const <- 1 / (b * n)

  return(list('ufun' = ufun,'upar' = upar, 'b' = b, 'const' = const))
}

loss_t_loss <- function(n, p, losspar, realdata){
  if(is.null(losspar)) upar <- 3 else upar <- losspar

  if(realdata && upar != 0){
    ufun <- function(t, v) 1 / (v + t)
    b <- tloss_consistency_factor(p, upar)
    const <- (upar + p) / (b * n)
  } else {
    if(!realdata && upar != 0){
      # complex case
      ufun <- function(t, v) 1 / (v + 2 * t)
      b <- tloss_consistency_factor(2 * p, upar)
      const <- (upar + 2 * p) / (b * n)
      } else {
        # Tylers M-Estimator
        ufun <- function(t, v) 1 / t
        const <- p / n
        b <- NULL
      }
    }
  return(list('ufun' = ufun, 'b' = b, 'const' = const, 'upar' = upar))
}

tloss_consistency_factor <- function(p, v){
  sfun <- function(x) x^(p / 2) / (v + x) * exp(-x / 2)

  c <- 2^(p / 2) * gamma(p / 2)
  q <- (1 / c) * integrate(sfun, 0, Inf)$value
  b <- ((v + p) / p) * q # consistency factor
  return(b)
}
