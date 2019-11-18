# % [beta, w] = elemfits(x,y)
# % elemfits compute the Nx(N-1)/2 elemental fits, i.e., intercepts b_{0,ij}
# % and slopes b_{1,ij}, that define a line y = b_0+b_1 x that passes through
# % the data points (x_i,y_i) and (x_j,y_j), i<j, where i, j in {1, ..., N}.
# % and the respective weights | x_i - x_j |
#   % INPUT:
#   %    y : (numeric) N x 1 vector of real-valued outputs (response vector)
# %    x : (numeric) N x 1 vector of inputs (feature vector)
# % OUTPUT:
#   %  beta: (numeric) N*(N-1)/2 matrix of elemental fits
# %     w: (numeric) N*(N-1)/2 matrix of weights
#' @export
elemfits <- function(y, x){
  N <- length(x)
  B <- matrix(1:N, N, N)
  a <- A[A<B]
  b <- B[A<B]
  w <- x[a] - x[b]

  beta <- matrix(0, 2, length(a))

  beta[2,] <- (y[a] - y[b])/w
  beta[1,] <- (x[a] * y[b] - x[b]*y[a])/w
  return (list(beta, abs(w)))
}

#' % enet computes the elastic net estimator using the cyclic co-ordinate
#' % descent (CCD) algorithm.
#' %
#' % INPUT:
#'   %       y : (numeric) data vector of size N x 1 (output, respones)
#' %         if the intercept is in the model, then y needs to be centered.
#' %       X : (numeric) data matrix of size N x p (input, features)
#' %           Columns are assumed to be standardized (i.e., norm(X(:,j))=1)
#' %           as well as centered (if intercept is in the model).
#' %    beta : (numeric) regression vector for initial start for CCD algorithm
#' %  lambda : (numeric) a postive penalty parameter value
#' %  alpha  : (numeric) elastic net tuning parameter in the range [0,1]. If
#' %           not given then use alpha = 1 (Lasso)
#' % printitn: print iteration number (default = 0, no printing)
#' % OUTPUT:
#'   %   beta    : (numberic) the regression coefficient vector
#' %   iter  : (numeric) # of iterations
#' @export
enet <- function(y, X, beta, lambda, alpha=1, printitn=0, itermax = 1000){
  # check for valid arguments

  p <- ncol(X)

  betaold <- beta
  normb0 <- norm(beta)
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

    normb <- norm(beta)
    crit <- sqrt(normb0^2 + normb^2 - 2 * Re(t(betaold) %*% beta))/normb

    if(!is.na(iter%%printitn) && (iter %% printitn) == 0){
      sprintf('enet: %4d  crit = %.8f\n',iter,crit)
    }

    if(crit<1e-4) break

    betaold <- beta
    normb0 <- normb

    return (list(beta, iter))
  }
}
