# enet ----

#'   enet
#'
#'   enet computes the elastic net estimator using the cyclic co-ordinate
#'   descent (CCD) algorithm.
#'
#'@param y      : (numeric) data vector of size N x 1 (output, respones)
#'                if the intercept is in the model, then y needs to be centered.
#'@param X      : (numeric) data matrix of size N x p (input, features)
#'                Columns are assumed to be standardized (i.e., norm(X(:,j))=1)
#'                as well as centered (if intercept is in the model).
#'@param beta   : (numeric) regression vector for initial start for CCD algorithm
#'@param lambda : (numeric) a postive penalty parameter value
#'@param alpha  : (numeric) elastic net tuning parameter in the range [0,1]. If
#'                 not given then use alpha = 1 (Lasso)
#'@param printitn: print iteration number (default = 0, no printing)
#'
#'@return beta  : (numeric) the regression coefficient vector
#'@return iter  : (numeric) # of iterations
#'@export
enet <- function(y, X, beta, lambda, alpha=1, printitn=0, itermax = 1000){
  # check for valid arguments

  p <- ncol(X)

  betaold <- beta
  normb0 <- norm(beta, type = "2")
  r <- y - X %*% beta

  lam1 <- alpha * lambda
  lam2 <- lambda * (1-alpha)

  const <- 1/(1+lam2)

  if(printitn > 0){
    cat('enet : using penalty lambda = ', round(lambda, digits=5))
  }

  for(iter in 1:itermax){
    for(jj in 1:p){
      beta[jj] <- const*soft_thresh(beta[jj]+ t(X[,jj]) %*% r, lam1)
      r <- r + X[,jj]*(betaold[jj]-beta[jj])
    }

    normb <- norm(beta, type = "2")

    crit <- sqrt(normb0^2 + normb^2 - 2 * Re(t(betaold) %*% beta))/normb

    if(!is.na(iter%%printitn) && (iter %% printitn) == 0){
      sprintf('enet: %4d  crit = %.8f\n',iter,crit)
    }

    if(is.nan(crit) | crit < 1e-4) {return (list(beta, iter)) }

    betaold <- beta
    normb0 <- normb
  }
  return(list(beta, iter))
}
