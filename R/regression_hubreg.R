
#'   hubreg
#'
#'   regression and scale using Huber's criterion
#'
#'   hubreg computes the joint M-estimates of regression and scale using
#'   Huber's criterion. Function works for both real- and complex-valued data.
#'
#'
#'@param   y: Numeric data vector of size N x 1 (output, respones)
#'@param   X: Numeric data matrix of size N x p. Each row represents one
#'             observation, and each column represents one predictor (feature).
#'             If the model has an intercept, then first column needs to be a
#'             vector of ones.
#'@param    c: numeric threshold constant of Huber's function
#'@param    sig0: (numeric) initial estimator of scale \cr
#'          default = SQRT(1/(n-p)*RSS)
#'@param    b0: initial estimator of regression (default: LSE)
#'@param    printitn: print iteration number (default = 0, no printing)
#'@param    iter_max: maximum number of iterations. \cr default = 2000
#'@param    errortol: ERROR TOLERANCE FOR HALTING CRITERION. \cr default = 1e-5
#'
#'@return    b1: the regression coefficient vector estimate
#'@return    sig1: the estimate of scale
#'@return    iter: the # of iterations
#'
#'@references
#'uses \code{\link[MASS]{ginv}} from the MASS package
#'
#'@examples
#'y <- c(1.0347, 0.7269, -0.3034, 0.2939, -0.7873)
#'X <- matrix(c(0.884, -1.1471, -1.0689, -0.8095, -2.9443, 1.4384, 0.3252, -0.7549, 1.3703, -1.7115), 5, 2)
#'
#'hubreg(y, X)
#'hubreg(y+1i, X)
#'hubreg(y+1i, X+1i)
#'hubreg(y, X+1i)
#'
#'@export
hubreg <- function(y, X, c = NULL, sig0 = NULL, b0 = NULL, printitn = 0, iter_max = 2000, errortol = 1e-5){
  n <- nrow(X)
  p <- ncol(X)

  if(is.complex(y)) real_data <- F else real_data <- T
  # Default: approx 95 efficiency for Gaussian errors
  if(is.null(c)){if(real_data) c <- 1.3415 else c <- 1.215}

  if(is.null(b0)) b0 <- qr.solve(X, y) # b0 <- ginv(X) %*% y

  # use unbiased residual sum of squares as initial estimate of scale
  if(is.null(sig0)){sig0 <- norm(y - X %*% b0, type = "2") /
    sqrt(n - p)}

  csq <- c^2

  # consistency factor for scale
  if(real_data){
    qn <- pchisq(csq, 1)
    alpha <- pchisq(csq, 3) + csq * (1 - qn)
  } else{
    qn <- pchisq(2 * csq, 2)
    alpha <- pchisq(2 * csq, 4) + csq * (1 - qn)
  }

  Z <- ginv(X)

  con <- sqrt((n - p) * alpha)

  for(iter in 1:iter_max){
    # Step 1: Update residual
    r <- y - X %*% b0
    psires <- psihub(r / sig0, c) * sig0

    # Step 2: Update the scale
    sig1 <- norm(psires, type ="2") / con

    # Step 3: Update the pseudo-residual
    psires <- psihub(r / sig1, c) * sig1

    # Step 4: regress X on pseudo-residual
    update <- Z %*% psires

    # Step 5: update the beta
    b1 <- b0 + update

    # Step 6 check convergence
    crit2 <- norm(update, type = "2") / norm(b0, type = "2")

    if(printitn > 0 & iter %% printitn) sprintf('hubreg: crit(%4d) = %.9f\n',iter,crit2)

    if(is.na(crit2) | crit2 < errortol) break

    b0 <- b1
    sig0 <- sig1

    if(iter == iter_max) sprintf('error!!! MAXiter = %d crit2 = %.7f\n',iter,crit2)
  }

  return( list(b1, sig1, iter))
}
